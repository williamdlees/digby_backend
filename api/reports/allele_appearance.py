# Allele Appearances report

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, chunk_list
from app import app, vdjbase_dbs, genomic_dbs
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype, AllelesSample, Gene, Allele
from db.vdjbase_airr_model import Sample, Patient
from db.genomic_db import Sequence as GenomicSequence, Subject as GenomicSubject, SubjectSequence as GenomicSubjectSequence, Gene as GenomicGene
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
import xlwt

APPEARANCE_SCRIPT = 'allele_appeareance2.R'
SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format != 'pdf' and format != 'xls':
        raise BadRequest('Invalid format requested')

    rep_samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in rep_samples_by_dataset:
            rep_samples_by_dataset[rep_sample['dataset']] = []
        rep_samples_by_dataset[rep_sample['dataset']].append(rep_sample['sample_name'])

    # Format we need to produce is [gene_name, [allele names], [allele appearances], gene appearances]
    # Start with a dict indexed by gene, then convert to appropriately sorted list
    counts = {}

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        appearances = []

        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            app_query = session.query(AllelesSample.patient_id, Patient.patient_name, Gene.name, Allele.name, Sample.sample_name, Gene.locus_order, Gene.alpha_order)\
                                .filter(Sample.id == AllelesSample.sample_id)\
                                .filter(Allele.id == AllelesSample.allele_id)\
                                .filter(Gene.id == Allele.gene_id)\
                                .filter(Patient.id == AllelesSample.patient_id)\
                                .filter(Sample.sample_name.in_(sample_list))\
                                .filter(Gene.name.in_(wanted_genes))

            if params['novel_alleles'] == 'Exclude':
                app_query = app_query.filter(Allele.novel == 0)

            if params['ambiguous_alleles'] == 'Exclude':
                app_query = app_query.filter(Allele.is_single_allele == 1)

            appearances.extend(app_query.all())

        for app in appearances:
            _, patient_name, gene, allele, sample, locus_order, alpha_order = app
            allele = allele.split('*', 1)[1].upper()
            if gene not in counts:
                if params['sort_order'] == 'Alphabetic':
                    counts[gene] = [{}, [], alpha_order]
                else:
                    counts[gene] = [{}, [], locus_order]
            if allele not in counts[gene][0]:
                counts[gene][0][allele] = []
            if patient_name not in counts[gene][0][allele]:
                counts[gene][0][allele].append(patient_name)
            if patient_name not in counts[gene][1]:
                counts[gene][1].append(patient_name)

    gen_samples_by_dataset = {}
    for gen_sample in genomic_samples:
        if gen_sample['dataset'] not in gen_samples_by_dataset:
            gen_samples_by_dataset[gen_sample['dataset']] = []
        gen_samples_by_dataset[gen_sample['dataset']].append(gen_sample['identifier'])

    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].session
        appearances = []

        for sample_chunk in chunk_list(gen_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(GenomicSubject.identifier).filter(GenomicSubject.identifier.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            app_query = session.query(GenomicSubject.identifier,
                                      GenomicGene.name,
                                      GenomicSequence.name,
                                      GenomicGene.locus_order,
                                      GenomicGene.alpha_order)\
                .filter(GenomicSubject.id == GenomicSubjectSequence.subject_id) \
                .filter(GenomicSequence.id == GenomicSubjectSequence.sequence_id) \
                .filter(GenomicGene.id == GenomicSequence.gene_id) \
                .filter(GenomicSequence.type.in_(['V-REGION', 'D-REGION', 'J-REGION'])) \
                .filter(GenomicSubject.identifier.in_(sample_list))\
                .filter(GenomicGene.name.in_(wanted_genes))

            if params['novel_alleles'] == 'Exclude':
                app_query = app_query.filter(GenomicSequence.novel == 0)
            appearances.extend(app_query.all())

        for app in appearances:
            patient_name, gene, allele, locus_order, alpha_order = app
            allele = allele.split('*', 1)[1].upper()
            if gene not in counts:
                if params['sort_order'] == 'Alphabetic':
                    counts[gene] = [{}, [], alpha_order]
                else:
                    counts[gene] = [{}, [], locus_order]
            if allele not in counts[gene][0]:
                counts[gene][0][allele] = []
            if patient_name not in counts[gene][0][allele]:
                counts[gene][0][allele].append(patient_name)
            if patient_name not in counts[gene][1]:
                counts[gene][1].append(patient_name)

    single_alleles = []
    multi_alleles = []

    for gene, (alleles, total, order) in counts.items():
        row = [gene, sorted(list(alleles.keys())), [len(alleles[a]) for a in sorted(alleles.keys())], len(total), order]
        if len(alleles) > 1:
            multi_alleles.append(row)
        else:
            single_alleles.append(row)

    multi_alleles.sort(key=lambda row: row[4])
    multi_alleles = [m[:4] for m in multi_alleles]
    single_alleles.sort(key=lambda row: row[4])

    s = ['Single allele genes', [], []]
    for (gene, alleles, counts, _, _) in single_alleles:
        s[1].append(gene + '\n' + alleles[0])
        s[2].append(counts[0])

    multi_alleles.append(s)

    input_path = make_output_file('xls')
    output_path = make_output_file('pdf')
    book = xlwt.Workbook()

    for row in multi_alleles:
        if len(row[1]) > 0:
            write_gene(book, row)

    book.save(input_path)

    if format == 'xls':
        return send_report(input_path, format, '%s_allele_appearance.xls' % species)

    cmd_line = ["-i", input_path,
                "-o", output_path]

    if run_rscript(APPEARANCE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, '%s_allele_appearance.pdf' % species)
    else:
        raise BadRequest('No output from report')


def write_gene(book, gene):
    sheet = book.add_sheet(gene[0].replace("/","-"))

    sheet.write(0,0, "Allele")
    sheet.write(0,1, "Appeared")

    for i in range (0, len(gene[1])):
        sheet.write(i + 1, 0, gene[1][i])
        sheet.write(i + 1, 1, gene[2][i])

    if len(gene) == 4:
        sheet.write(0,2, "Total")
        sheet.write(1,2, gene[3])

