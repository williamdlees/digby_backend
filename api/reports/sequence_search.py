# Sequence search for AIRR-seq and genomic samples

from werkzeug.exceptions import BadRequest
from api.reports.report_utils import make_output_file, collate_samples, chunk_list, collate_gen_samples, splitlines, chunks
from api.reports.reports import send_report

from app import app, vdjbase_dbs, genomic_dbs
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype, AllelesSample, Gene, Allele, AllelesPattern
from db.genomic_db import Sequence as GenomicSequence, Subject as GenomicSubject, SubjectSequence as GenomicSubjectSequence, Gene as GenomicGene

from db.vdjbase_airr_model import Patient, Sample
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
from Bio import pairwise2
from receptor_utils import simple_bio_seq as simple

SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format not in ["html", "xls"]:
        raise BadRequest('Invalid format requested')

    r_chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)

    seqs_to_search = {}

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        allele_recs = []

        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            query = session.query(Allele) \
                .join(Gene) \
                .join(AllelesSample) \
                .join(Sample) \
                .join(Patient, Patient.id == Sample.patient_id) \
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Allele.name.notlike('%Del%')) \
                .filter(Allele.name.notlike('%OR%')) \
                .filter(Sample.sample_name.in_(sample_list))

            if params['novel_alleles'] == 'Exclude':
                query = query.filter(Allele.novel == 0)

            if params['ambiguous_alleles'] == 'Exclude':
                query = query.filter(Allele.is_single_allele == 1)

            for allele in query.all():
                if allele not in allele_recs:
                    allele_recs.append(allele)

        for allele in allele_recs:
            seq = allele.seq.replace('.', '')
            if seq not in seqs_to_search:
                seqs_to_search[seq] = {}
                seqs_to_search[seq]['datasets'] = []
                seqs_to_search[seq]['sequence'] = seq
            seqs_to_search[seq]['datasets'].append({
                'type': 'AIRR-seq',
                'dataset': dataset,
                'appearances': allele.appears,
                'gene': allele.gene.name,
                'allele_name': allele.name,
            })

    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].session
        sequence_recs = []

        for sample_chunk in chunk_list(gen_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(GenomicSubject.identifier).filter(GenomicSubject.identifier.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            query = session.query(GenomicSequence)\
                .join(GenomicGene)\
                .join(GenomicSubjectSequence)\
                .join(GenomicSubject) \
                .filter(GenomicSubject.id == GenomicSubjectSequence.subject_id) \
                .filter(GenomicSequence.id == GenomicSubjectSequence.sequence_id) \
                .filter(GenomicGene.id == GenomicSequence.gene_id)\
                .filter(GenomicSequence.type.in_(['V-REGION', 'D-REGION', 'J-REGION'])) \
                .filter(GenomicSubject.identifier.in_(sample_list))\
                .filter(GenomicGene.name.in_(wanted_genes))

            if params['novel_alleles'] == 'Exclude':
                query = query.filter(GenomicSequence.novel == 0)

            for sequence in query.all():
                if sequence not in sequence_recs:
                    sequence_recs.append(sequence)

        for sequence in sequence_recs:
            seq = sequence.sequence.lower()
            if seq not in seqs_to_search:
                seqs_to_search[seq] = {}
                seqs_to_search[seq]['datasets'] = []
                seqs_to_search[seq]['sequence'] = seq
            seqs_to_search[seq]['datasets'].append({
                'type': 'Genomic',
                'dataset': dataset,
                'appearances': sequence.appearances,
                'gene': sequence.gene.name,
                'allele_name': sequence.name,
            })

    target_seq = params['sequence'].lower().replace('.', '').replace('\n', '')
    target_score = len(target_seq) * 0.9

    for seq in seqs_to_search.keys():
        seqs_to_search[seq]['score'] = pairwise2.align.globalms(target_seq, seq, 1, -0.5, -0.5, -0.1, score_only=True)

    results = [x for x in seqs_to_search.values() if x['score'] >= target_score]

    if results:
        results = sorted(results, key=lambda x: x['score'], reverse=True)


    if format == 'html':
        target = params['sequence'].replace('\n', '').replace('\r', '').lower()

        while target[-1] == '':
            target = target[:-1]

        doc = '<html><head><title>Sequence Search Results</title></head><body>'
        doc += f'Target sequence: <pre><code>'
        for line in chunks(target, 80):
            doc += line + '\n'

        doc += '</code></pre><br>'

        for result in results:
            alignment = pairwise2.align.globalms(target_seq, result['sequence'], 1, -0.5, -0.5, -0.1, one_alignment_only=True)
            pretty = pairwise2.format_alignment(*alignment[0])

            pretty = pretty.replace('\r', '')
            pretty = pretty.split('\n')
            rep = pretty[:-2]

            maxlen = max([len(x) for x in rep])
            for i in range(len(rep)):
                rep[i] = rep[i].rjust(maxlen)

            rep = splitlines('\n'.join(rep), 80, 0)

            headers = []
            for ds in result['datasets']:
                headers.append(f"{ds['type']} {ds['dataset']} {ds['allele_name']} ({ds['appearances']} appearances)")

            headers = '; '.join(headers)
            doc += f'<pre><code>\n\n{headers}\n{pretty[-2]}\n{rep}<code></pre>'

        if not results:
            doc += '<div>No results found with >= 90% sequence identity.</div>'

        doc += '</body></html>'

        output_path = make_output_file('html')
        with open(output_path, 'w') as fo:
            fo.write(doc)

    if format == 'xls':
        doc = []

        for result in results:
            for ds in result['datasets']:
                doc.append({
                    'Score': result['score'],
                    'Dataset': ds['dataset'],
                    'Type': ds['type'],
                    'Gene': ds['gene'],
                    'Allele': ds['allele_name'],
                    'Appearances': ds['appearances'],
                    'Allele Sequence': result['sequence'],
                    'Target Sequence': target_seq,
                })

        if not results:
            doc.append({
                'Score': '',
                'Dataset': '',
                'Type': '',
                'Gene': '',
                'Allele': '',
                'Appearances': '',
                'Allele Sequence': '',
                'Target Sequence': target_seq,
            })

        output_path = make_output_file('csv')
        simple.write_csv(output_path, doc)

    return send_report(output_path, format, 'sequence_search.csv')
