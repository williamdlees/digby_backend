# Haplotype heatmap for VDJbase samples

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, collate_samples, find_primer_translations, translate_primer_alleles, translate_primer_genes, \
    collate_gen_samples

from api.reports.report_utils import trans_df
from app import app, vdjbase_dbs, genomic_dbs
from db.genomic_db import Subject as GenomicSubject, Gene as GenomicGene, Sequence as GenomicSequence, SubjectSequence as GenomicSubjectSequence

from db.vdjbase_airr_model import Sample
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, get_multiple_order_file
from api.genomic.genomic import GENOMIC_SAMPLE_PATH
import pandas as pd
from receptor_utils.simple_bio_seq import read_csv


MULTIPLE_GENOTYPE_SCRIPT = "html_multiple_genotype_hoverText.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)

    if not chain:
        chain = g_chain

    genotypes = []
    all_wanted_genes = set()

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        primer_trans, gene_subs = find_primer_translations(session)
        sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(rep_samples_by_dataset[dataset])).all()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        if len(wanted_genes) > 0:
            all_wanted_genes |= set(wanted_genes)
            for (name, genotype, patient_id) in sample_list:
                sample_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, genotype.replace('samples/', ''))

                if not os.path.isfile(sample_path):
                    continue        # this protects against missing genotype files, usually caused by the pipeline omitting samples listed in the yaml file

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

    if len(genotypes) > 20:
        raise BadRequest('Please select at most 20 genotypes, or use the Genotype Heatmap report.')


    # add fakes to each genotype for missing genes

    for i in range(len(genotypes)):
        genotype = genotypes[i]
        subject_name = genotype.iloc[0]['subject']
        contained_genes = genotype['gene'].tolist()
        fakes = []
        for gene in all_wanted_genes:
            if gene not in contained_genes:
                fakes.append(fake_gene([], gene, subject_name))

        genotypes[i] = pd.concat([genotype, pd.DataFrame(fakes)], ignore_index=True)

    geno_path = make_output_file('tsv')
    genotypes = pd.concat(genotypes)
    genotypes.to_csv(geno_path, sep='\t')

    if format == 'pdf':
        attachment_filename = '%s_sampled_genotype.pdf' % (species)
    else:
        attachment_filename = None

    if not params['f_pseudo_genes']:
        pseudo = 'F'
    else:
        pseudo = 'T'

    locus_order = ('sort_order' in params and params['sort_order'] == 'Locus')
    gene_order_file = get_multiple_order_file(species, rep_samples_by_dataset.keys(), gen_samples_by_dataset.keys(), locus_order=locus_order)

    output_path = make_output_file('html' if html else 'pdf')

    file_type = 'T' if html else 'F'
    cmd_line = ["-i", geno_path,
                "-o", output_path,
                "-t", file_type,
                "-g", gene_order_file,
                "-c", chain]

    if run_rscript(MULTIPLE_GENOTYPE_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')


def fake_genotype(subject, genes):
    genotype = []
    for gene_name, alleles in genes.items():
        rec = fake_gene(alleles, gene_name, subject)
        genotype.append(rec)

    return genotype


def fake_gene(alleles, gene_name, subject):
    rec = {
        'subject': subject,
        'gene': gene_name,
        'alleles': ','.join(alleles),
        'counts': ','.join(['100' for _ in range(len(alleles))]),
        'total': f"{100 * len(alleles)}",
        'note': '',
        'kh': '1',
        'kd': '1',
        'kt': '1',
        'kq': '1',
        'k_diff': '1',
        'GENOTYPED_ALLELES': ','.join(alleles),
        'Freq_by_Clone': ';'.join(['1' for _ in range(len(alleles))]),
        'Freq_by_Seq': ';'.join(['1' for _ in range(len(alleles))]),
    }
    return rec

# why do we have this? Why not use process_vdjbase_genotype?
def process_igenotyper_genotype(sample_path, subject, study, wanted_genes):
    genotype = read_csv(sample_path)
    genes = {}

    # collect alleles and gene names
    for rec in genotype:
        gene = rec['vdjbase_allele'].split('*')[0]
        if gene not in genes:
            genes[gene] = []

        g = rec['vdjbase_allele'].split('*')[1]
        if g not in genes[gene]:
            genes[gene].append(g)

    geno_list = fake_genotype(subject, genes)
    return pd.DataFrame(geno_list)


def process_vdjbase_genotype(subject_name, wanted_genes, session, functional):
    alleles = session.query(GenomicSequence.name, GenomicGene.name)\
        .join(GenomicGene)\
        .join(GenomicSubjectSequence)\
        .join(GenomicSubject)\
        .filter(GenomicSubject.identifier == subject_name)\
        .filter(GenomicGene.name.in_(wanted_genes))\
        .filter(GenomicSubjectSequence.sequence_id == GenomicSequence.id)\
        .filter(GenomicSubjectSequence.subject_id == GenomicSubject.id)

    if functional:
        alleles = alleles.filter(GenomicSequence.functional == 'Functional')

    alleles = alleles.all()

    genes = {}
    for allele, gene in alleles:
        if gene not in genes:
            genes[gene] = []

        if '*' in allele:
            allele_name = allele.split('*')[1]
        elif '.' in allele:         # cirelli format
            if allele[-2:-1] == '.' and allele[-1].isalpha():
                allele_name = '01_' + allele[-1]
            elif '_' in allele.replace('LJI.Rh_', ''):
                allele_name = '01_' + allele.replace('LJI.Rh_', '').split('_')[-1]
            else:
                allele_name = '01'
        else:
            allele_name = allele    # don't expect this to happen

        print(f"{allele} -> {gene}*{allele_name}")

        if allele_name not in genes[gene]:
            genes[gene].append(allele_name)

    geno_list = fake_genotype(subject_name, genes)
    return pd.DataFrame(geno_list)
