# Allele usage plot for AIRR-seq samples

from werkzeug.exceptions import BadRequest
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file, collate_samples, chunk_list, collate_gen_samples

from app import vdjbase_dbs, genomic_dbs
from db.vdjbase_model import AllelesSample, Gene, Allele, AllelesPattern
from db.genomic_db import Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence, Gene as GenomicGene
from db.genomic_airr_model import Sample as GenomicSample

from db.vdjbase_airr_model import Patient, Sample
import os
from api.vdjbase.vdjbase import apply_rep_filter_params
import pandas as pd

ALLELE_USAGE_SCRIPT = 'Alleles_Usage.R'
SAMPLE_CHUNKS = 400



def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format != 'html':
        raise BadRequest('Invalid format requested')

    kdiff = float(params['f_kdiff']) if 'f_kdiff' in params and params['f_kdiff'] != '' else 0
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)

    if (chain and g_chain) and chain != g_chain:
        raise BadRequest('This report requires all samples to be selected from the same chain (IGH, IGK, ...')

    if not chain:
        chain = g_chain

    if len(rep_samples_by_dataset) + len(gen_samples_by_dataset) > 1 and params['ambiguous_alleles'] != 'Exclude':
        raise BadRequest('Ambiguous alleles cannot be processed across multiple datasets')

    # Format we need to produce is [(gene_name, hetero count, homo count),...]

    gene_allele_counts = {}

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].get_session()
        allele_recs = []

        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            query = session.query(Gene.name, Allele.name, Gene.type) \
                .join(Allele) \
                .join(AllelesSample) \
                .join(Sample) \
                .join(Patient, Patient.id == Sample.patient_id) \
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Allele.name.notlike('%Del%')) \
                .filter(Allele.name.notlike('%OR%')) \
                .filter(Sample.sample_name.in_(sample_list)) \
                .filter(AllelesSample.kdiff >= kdiff)

            if 'sort_order' in params and params['sort_order'] == 'Locus':
                query = query.order_by(Gene.locus_order, Patient.id, Allele.name)
            else:
                query = query.order_by(Gene.alpha_order, Patient.id, Allele.name)

            if params['novel_alleles'] == 'Exclude':
                query = query.filter(Allele.novel == 0)

            if params['ambiguous_alleles'] == 'Exclude':
                query = query.filter(Allele.is_single_allele == 1)

            allele_recs.extend(query.all())

        # Collect up ids for each gene. Note this assumes a sorted list of genes (the order_by above)

        i = 0
        while i < len(allele_recs):
            (gene_name, allele_name, gene_type) = allele_recs[i]
            gene_allele_names = []

            while i < len(allele_recs):
                if allele_recs[i][0] != gene_name:
                    break

                allele_name = allele_recs[i][1]
                gene_allele_names.append(allele_name)
                i += 1

            gene_allele_names = set(gene_allele_names)

            # If we have both an unambiguous allele and an ambiguous allele containing that unambiguous one,
            # drop the unambiguous one because it is already counted

            if(params['ambiguous_alleles'] != 'Exclude'):
                patterns = session.query(AllelesPattern.pattern_id)\
                    .filter(AllelesPattern.allele_in_p_id.in_(gene_allele_names))\
                    .filter(AllelesPattern.pattern_id.in_(gene_allele_names))\
                    .all()

                if patterns is not None and len(patterns) >0:
                    patterns = set([pattern[0] for pattern in patterns])
                    gene_allele_names = gene_allele_names - patterns

            if gene_name not in gene_allele_counts:
                gene_allele_counts[gene_name] = gene_allele_names
            else:
                gene_allele_counts[gene_name] |= gene_allele_names

    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].get_session()
        allele_recs = []

        for sample_chunk in chunk_list(gen_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(GenomicSample.sample_name).filter(GenomicSample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            query = session.query(GenomicGene.name, GenomicSequence.name, GenomicGene.type) \
                .filter(GenomicSample.id == GenomicSampleSequence.sample_id) \
                .filter(GenomicSequence.id == GenomicSampleSequence.sequence_id) \
                .filter(GenomicGene.id == GenomicSequence.gene_id)\
                .filter(GenomicSequence.type.in_(['V-REGION', 'D-REGION', 'J-REGION'])) \
                .filter(GenomicSample.sample_name.in_(sample_list))\
                .filter(GenomicGene.name.in_(wanted_genes))

            if 'sort_order' in params and params['sort_order'] == 'Locus':
                query = query.order_by(GenomicGene.locus_order, GenomicSample.sample_name, GenomicSequence.name)
            else:
                query = query.order_by(GenomicGene.alpha_order, GenomicSample.sample_name, GenomicSequence.name)

            if params['novel_alleles'] == 'Exclude':
                query = query.filter(GenomicSequence.novel == 0)

            if not params['f_pseudo_genes']:
                query = query.filter(GenomicSequence.functional == 'Functional')

            allele_recs.extend(query.all())

        # Collect up ids for each gene. Note this assumes a sorted list of genes (the order_by above)

        i = 0
        while i < len(allele_recs):
            (gene_name, allele_name, gene_type) = allele_recs[i]
            gene_allele_ids = []

            while i < len(allele_recs):
                if allele_recs[i][0] != gene_name:
                    break

                allele_name = allele_recs[i][1]
                gene_allele_ids.append(allele_name)
                i += 1

            gene_allele_ids = set(gene_allele_ids)

            if gene_name not in gene_allele_counts:
                gene_allele_counts[gene_name] = gene_allele_ids
            else:
                gene_allele_counts[gene_name] |= gene_allele_ids

    listed_allele_count = []
    for gene, alleles in gene_allele_counts.items():
        listed_allele_count.append((gene, len(alleles)))

    labels = ['GENE', 'COUNT']
    input_path = make_output_file('tab')
    df = pd.DataFrame(listed_allele_count, columns=labels)
    df.to_csv(input_path, sep='\t', index=False)
    output_path = make_output_file('html')

    cmd_line = ["-i", input_path,
                "-o", output_path,
                "-c", chain
                ]

    if run_rscript(ALLELE_USAGE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format)
    else:
        raise BadRequest('No output from report')

