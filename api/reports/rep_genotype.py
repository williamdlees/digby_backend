# Genotype List for AIRR-seq and Genomic samples

import os
import pandas as pd
from werkzeug.exceptions import BadRequest
from receptor_utils.simple_bio_seq import read_csv

from app import vdjbase_dbs, genomic_dbs
from api.reports.reports import run_rscript, send_report
from api.reports.report_utils import make_output_file, collate_samples, collate_gen_samples
from db.genomic_db import Subject as GenomicSubject, Gene as GenomicGene, Sequence as GenomicSequence, SubjectSequence as GenomicSubjectSequence
from db.vdjbase_airr_model import Sample
from db.vdjbase_model import Gene, Allele, AllelesSample
from api.vdjbase.vdjbase import apply_rep_filter_params, get_multiple_order_file


MULTIPLE_GENOTYPE_SCRIPT = "html_multiple_genotype_hoverText.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format not in ['pdf', 'html']:
        raise BadRequest('Invalid format requested')

    html = (format == 'html')
    chain, rep_samples_by_dataset = collate_samples(rep_samples)
    g_chain, gen_samples_by_dataset = collate_gen_samples(genomic_samples)

    if not chain:
        chain = g_chain

    genotypes = {}
    all_wanted_genes = set()

    for dataset in rep_samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(rep_samples_by_dataset[dataset])).all()
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
                genotype = process_genomic_genotype(subject_name, wanted_genes, session, not params['f_pseudo_genes'])
                genotypes[subject_name] = genotype

    if len(genotypes) == 0:
        raise BadRequest('No records matching the filter criteria were found.')

    if len(genotypes) > 20:
        raise BadRequest('Please select at most 20 genotypes, or use the Genotype Heatmap report.')


    # add fakes to each genotype for missing genes

    for subject_name, genotype in genotypes.items():
        contained_genes = genotype['gene'].tolist() if len(genotype) > 0 else []

        fakes = []
        for gene in all_wanted_genes:
            if gene not in contained_genes:
                fakes.append(fake_gene({'alleles': [], 'count': [], 'fc': [], 'fs': [], 'total_count': 0, 'kdiff': 0}, gene, subject_name))

        genotypes[subject_name] = pd.concat([genotype, pd.DataFrame(fakes)], ignore_index=True)

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
    for gene_name, gene_details in genes.items():
        rec = fake_gene(gene_details, gene_name, subject)
        genotype.append(rec)

    return genotype


def fake_gene(gene_details, gene_name, subject):
    rec = {
        'subject': subject,
        'gene': gene_name,
        'alleles': ','.join(gene_details['alleles']),
        'counts': ','.join(gene_details['count']),
        'total': gene_details['total_count'],
        'note': '',
        'kh': '1',
        'kd': '1',
        'kt': '1',
        'kq': '1',
        'k_diff': gene_details['kdiff'],
        'GENOTYPED_ALLELES': ','.join(gene_details['alleles']),
        'Freq_by_Clone': ';'.join(gene_details['fc']),
        'Freq_by_Seq': ';'.join(gene_details['fs']),
    }
    return rec

# why do we have this? Why not use process_genomic_genotype?
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


def process_genomic_genotype(subject_name, wanted_genes, session, functional):
    allele_query = session.query(GenomicSubject.identifier, GenomicSequence.name, GenomicGene.name)\
        .filter(GenomicSubject.identifier == subject_name) \
        .join(GenomicSubjectSequence, GenomicSubjectSequence.subject_id == GenomicSubject.id) \
        .join(GenomicSequence, GenomicSubjectSequence.sequence_id == GenomicSequence.id) \
        .filter(GenomicSequence.type.like('%REGION%')) \
        .join(GenomicGene, GenomicGene.id == GenomicSequence.gene_id)

    if len(wanted_genes) > 0:
        allele_query = allele_query.filter(GenomicGene.name.in_(wanted_genes))

    if functional:
        allele_query = allele_query.filter(GenomicSequence.functional == 'Functional')

    alleles = allele_query.all()

    genes = {}
    for _, allele, gene in alleles:
        if gene not in genes:
            genes[gene] = {'alleles': [], 'count': [], 'fc': [], 'fs': []}

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

        if allele_name not in genes[gene]['alleles']:
            genes[gene]['alleles'].append(allele_name)
            genes[gene]['count'].append('1')
            genes[gene]['fc'].append('1')
            genes[gene]['fs'].append('1')
            genes[gene]['kdiff'] = '1000'
            genes[gene]['total_count'] = 1
        else:
            print(f"Allele {allele_name} seen multiply in genotype")

    geno_list = fake_genotype(subject_name, genes)
    return pd.DataFrame(geno_list)


def process_repseq_genotype(sample_name, wanted_genes, session, functional):
    alleles = session.query(Sample.sample_name, Allele.name, Gene.name, AllelesSample) \
        .filter(Sample.sample_name == sample_name) \
        .join(AllelesSample, AllelesSample.sample_id == Sample.id) \
        .join(Allele, AllelesSample.allele_id == Allele.id) \
        .join(Gene, Allele.gene_id == Gene.id)

    if len(wanted_genes) > 0:
        alleles = alleles.filter(Gene.name.in_(wanted_genes))

    if functional:
        alleles = alleles.filter(Gene.pseudo_gene == False)

    alleles = alleles.all()

    genes = {}
    for sample_name, allele, gene, allelessample in alleles:
        if gene not in genes:
            genes[gene] = {'alleles': [], 'count': [], 'fc': [], 'fs': []}

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

        if allele_name not in genes[gene]['alleles']:
            genes[gene]['alleles'].append(allele_name)
            genes[gene]['count'].append(str(allelessample.count))
            genes[gene]['fc'].append(str(allelessample.freq_by_clone))
            genes[gene]['fs'].append(str(allelessample.freq_by_seq))
            genes[gene]['kdiff'] = str(allelessample.kdiff)
            genes[gene]['total_count'] = str(allelessample.total_count)
        else:
            print(f"Allele {allele_name} seen multiply in genotype")

    geno_list = fake_genotype(sample_name, genes)
    return pd.DataFrame(geno_list)
