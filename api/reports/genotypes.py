# Functions to extract genotypes by database query

from db.genomic_db import Subject as GenomicSubject, Sample as GenomicSample, Gene as GenomicGene, Sequence as GenomicSequence, SampleSequence as GenomicSampleSequence, SampleSequence as GenomicSampleSequence
from db.genomic_airr_model import Sample as GenomicSample, Sequence as GenomicSequence, G
from db.vdjbase_model import Gene, Allele, AllelesSample
from db.vdjbase_airr_model import Sample
from receptor_utils.simple_bio_seq import read_csv

import pandas as pd

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


def process_genomic_genotype(sample_name, wanted_genes, session, functional, fully_haplotyped):
    allele_query = session.query(GenomicSample.identifier, GenomicSequence.name, GenomicGene.name, GenomicSampleSequence.haplotype)\
        .filter(GenomicSample.identifier == sample_name) \
        .join(GenomicSampleSequence, GenomicSampleSequence.sample_id == GenomicSample.id) \
        .join(GenomicSequence, GenomicSampleSequence.sequence_id == GenomicSequence.id) \
        .filter(GenomicSequence.type.like('%REGION%')) \
        .join(GenomicGene, GenomicGene.id == GenomicSequence.gene_id)

    if len(wanted_genes) > 0:
        allele_query = allele_query.filter(GenomicGene.name.in_(wanted_genes))

    if functional:
        allele_query = allele_query.filter(GenomicSequence.functional == 'Functional')

    alleles = allele_query.all()

    genes = {}
    for _, allele, gene, haplotype in alleles:
        if gene not in genes:
            genes[gene] = {'alleles': [], 'count': [], 'fc': [], 'fs': [], 'haplotypes': []}

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
            genes[gene]['haplotypes'].extend(haplotype.split(','))
            genes[gene]['count'].append('1')
            genes[gene]['fc'].append('1')
            genes[gene]['fs'].append('1')
            genes[gene]['kdiff'] = '1000'
            genes[gene]['total_count'] = 1
        else:
            print(f"Allele {allele_name} seen multiply in genotype")

    if fully_haplotyped:
        dropped = []
        for gene_name in list(genes.keys()):
            if 'h0' in genes[gene_name]['haplotypes'] or len(genes[gene_name]['haplotypes']) < 2:
                dropped.append(genes[gene_name])
                del genes[gene_name]

    geno_list = fake_genotype(sample_name, genes)
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
