# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest

from api.genomic.genomic import GENOMIC_SAMPLE_PATH
from api.reports.rep_genotype import fake_gene, process_igenotyper_genotype, process_vdjbase_genotype
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
    genotypes = []
    all_wanted_genes = set()

    if not chain:
        chain = g_chain

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        primer_trans, gene_subs = find_primer_translations(session)

        sample_list = []
        for sample_chunk in chunk_list(rep_samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list.extend(session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all())

        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            all_wanted_genes |= set(wanted_genes)
            for (name, genotype, patient_id) in sample_list:
                sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, genotype.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    continue

                genotype = pd.read_csv(sample_path, sep='\t', dtype=str)

                genotype = trans_df(genotype)

                # translate pipeline allele names to VDJbase allele names
                for col in ['alleles', 'GENOTYPED_ALLELES']:
                    genotype[col] = [translate_primer_alleles(x, y, primer_trans) for x, y in zip(genotype['gene'], genotype[col])]

                genotype['gene'] = [translate_primer_genes(x, gene_subs) for x in genotype['gene']]
                genotype = genotype[genotype.gene.isin(wanted_genes)]

                subject_name = name if len(rep_samples_by_dataset) == 1 else dataset + '_' + name

                if 'subject' not in genotype.columns.values:
                    genotype.insert(0, 'subject', subject_name)
                else:
                    genotype.subject = subject_name

                genotypes.append(genotype)

    for dataset in gen_samples_by_dataset.keys():
        session = genomic_dbs[species][dataset].session
        subjects = session.query(GenomicSubject).filter(GenomicSubject.identifier.in_(gen_samples_by_dataset[dataset])).all()

        if not subjects:
            continue

        annotation_method = subjects[0].annotation_method
        sample_list = [(subject.identifier, subject.annotation_path, subject.identifier) for subject in subjects]
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
        all_wanted_genes |= set(wanted_genes)

        # now that we have the filtered samples, build a richer list for genomic processing
        filtered_samples = [s[0] for s in sample_list]
        sample_list = [(subject.identifier, subject.name_in_study, subject.study.name, subject.annotation_path) for subject in subjects if subject.identifier in filtered_samples]

        if len(wanted_genes) > 0:
            for (subject_name, name_in_study, study_name, annotation_path) in sample_list:
                sample_path = os.path.join(GENOMIC_SAMPLE_PATH, annotation_path.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    continue        # this protects against missing genotype files, usually caused by the pipeline omitting samples listed in the yaml file

                if annotation_method == 'IGenotyper':
                    genotype = process_igenotyper_genotype(sample_path, subject_name, study_name, wanted_genes)
                elif annotation_method == 'VDJbase':
                    genotype = process_vdjbase_genotype(subject_name, wanted_genes, session, not params['f_pseudo_genes'])

                genotypes.append(genotype)

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    # add fakes to each genotype for missing genes

    for i in range(len(genotypes)):
        genotype = genotypes[i]
        subject_name = genotype.iloc[0]['subject']
        contained_genes = genotype['gene'].tolist()
        fakes = []
        for gene in all_wanted_genes:
            if gene not in contained_genes:
                fakes.append(fake_gene(['Deletion'], gene, subject_name))
        genotypes[i] = pd.concat([genotype, pd.DataFrame(fakes)], ignore_index=True)

    # sort rows in each genotype by gene
    for i in range(len(genotypes)):
        genotype = genotypes[i]
        genotype.sort_values(by=['gene'], inplace=True)

    geno_path = make_output_file('tsv')
    genotypes = pd.concat(genotypes)
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



