# Haplotype heatmap for AIRR-seq and Genomic samples

from werkzeug.exceptions import BadRequest


from api.reports.rep_genotype import fake_gene, process_igenotyper_genotype, process_genomic_genotype, process_repseq_genotype
from api.reports.report_utils import trans_df, collate_samples, chunk_list, collate_gen_samples
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file, find_primer_translations, translate_primer_alleles, translate_primer_genes
from app import vdjbase_dbs, genomic_dbs
from db.genomic_db import Subject as GenomicSubject

from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, get_multiple_order_file
import pandas as pd


HEATMAP_GENOTYPE_SCRIPT = "genotype_heatmap.R"
SAMPLE_CHUNKS = 400



def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    html = (format == 'html')
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)
    genotypes = {}
    all_wanted_genes = set()

    if not chain:
        chain = g_chain

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session

        sample_list = []
        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list.extend(session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all())

        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            all_wanted_genes |= set(wanted_genes)
            for (name, genotype, patient_id) in sample_list:
                session = vdjbase_dbs[species][dataset].session
                genotype = process_repseq_genotype(name, all_wanted_genes, session, False)
                genotypes[name] = genotype

    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].session
        subjects = session.query(GenomicSubject).filter(GenomicSubject.identifier.in_(gen_samples_by_dataset[dataset])).all()

        if not subjects:
            continue

        sample_list = [(subject.identifier, subject.annotation_path, subject.identifier) for subject in subjects]
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
        all_wanted_genes |= set(wanted_genes)

        # now that we have the filtered samples, build a richer list for genomic processing
        filtered_samples = [s[0] for s in sample_list]
        sample_list = [(subject.identifier, subject.name_in_study, subject.study.name, subject.annotation_path) for subject in subjects if subject.identifier in filtered_samples]

        if len(wanted_genes) > 0:
            for (subject_name, name_in_study, study_name, annotation_path) in sample_list:
                genotype = process_genomic_genotype(subject_name, all_wanted_genes, session, not params['f_pseudo_genes'])
                genotypes[subject_name] = genotype

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    # add fakes to each genotype for missing genes

    for subject_name, genotype in genotypes.items():
        contained_genes = genotype['gene'].tolist() if len(genotype) > 0 else []

        fakes = []
        for gene in all_wanted_genes:
            if gene not in contained_genes:
                fakes.append(fake_gene({'alleles': [], 'count': [], 'fc': [], 'fs': [], 'total_count': 0, 'kdiff': 0}, gene, subject_name))

        genotypes[subject_name] = pd.concat([genotype, pd.DataFrame(fakes)], ignore_index=True)

    # sort rows in each genotype by gene
    for genotype in genotypes.values():
        genotype.sort_values(by=['gene'], inplace=True)

    geno_path = make_output_file('tsv')
    genotypes = pd.concat(genotypes.values())
    genotypes.to_csv(geno_path, sep='\t')

    if format == 'pdf':
        attachment_filename = '%s_genotype.pdf' % species
    else:
        attachment_filename = None

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, rep_samples_by_dataset.keys(), gen_samples_by_dataset.keys(), locus_order=locus_order)

    output_path = make_output_file('html' if html else 'pdf')
    file_type = 'T' if html else 'F'
    cmd_line = ["-i", geno_path,
                "-o", output_path,
                "-t", file_type,
                "-k", str(params['f_kdiff']),
                "-c", chain,
                "-g", gene_order_file
                ]

    if run_rscript(HEATMAP_GENOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')



