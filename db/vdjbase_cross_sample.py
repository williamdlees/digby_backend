#
# Cross-sample processing - once sample updates are complete
#
import sqlite3

import pandas as pd
import os, math

from sqlalchemy import func
from db.vdjbase_model import Allele, AllelesSample, Gene, GenesDistribution, AllelesPattern
from db.vdjbase_airr_model import Sample
import re
from db.vdjbase_formats import *
from sqlalchemy.sql import func
from db.vdjbase_projects import compound_genes


def update_alleles_appearance(session):
    """
    This function update the appearances count according to the number of patients.
    """
    result = ['Updating allele appearance counts']
    alleles = session.query(Allele)

    for allele in alleles:
        max_kdiff = session.query(func.max(AllelesSample.kdiff)).filter(AllelesSample.allele_id == allele.id).one_or_none()[0]
        allele.max_kdiff = max_kdiff if max_kdiff is not None else 0


        allele.appears = session.query(AllelesSample.patient_id).filter(AllelesSample.hap == 'geno').filter(AllelesSample.allele_id == allele.id).distinct().count()
        if allele.similar is not None and allele.similar != '':
            sims = allele.similar.split(', ')
            for sim in sims:
                sim = sim.replace('|', '')
                allele.appears += session.query(AllelesSample.patient_id)\
                    .join(Allele)\
                    .filter(AllelesSample.hap == 'geno')\
                    .filter(Allele.name.ilike(sim))\
                    .distinct().count()

    session.commit()
    return result

def calculate_gene_frequencies(ds_dir, session):
    # upload the gene frequencies of the samples by calculating them from the genotype file.

    result = ['Updating overall gene frequencies']
    genes = session.query(Gene)
    gene_ids_by_name = dict(zip([gene.study_title for gene in genes], [gene.id for gene in genes]))

    samples = session.query(Sample)

    for sample in samples:
        if sample.genotype is None or not os.path.isfile(os.path.join(ds_dir, sample.genotype)):
            print('skipping genotype processing for sample %s - no file' % sample.study_title)
            continue

        genotype = pd.read_csv(os.path.join(ds_dir, sample.genotype), sep='\t')

        frequencies_by_clone = {}
        frequencies_by_seq = {}

        # counting the appearance of the genes and the total of each family
        family_total_seq = {}
        family_total_clone = {}
        for gene, count_seq, count_clone in zip(genotype[GENE_COLUMN], genotype[FREQ_BY_SEQ], genotype[FREQ_BY_CLONE]):
            family = gene[:4]
            if family not in family_total_seq.keys():
                family_total_seq[family] = 0
                family_total_clone[family] = 0

            if isinstance(count_seq, float):
                if math.isnan(count_seq):
                    count_seq = 0
                count_seq = str(count_seq)

            if isinstance(count_clone, float):
                if math.isnan(count_clone):
                    count_clone = 0
                count_clone = str(count_clone)

            frequencies_by_seq[gene] = 0
            for x in str(count_seq).split(INT_SEP):
                frequencies_by_seq[gene] += int(float(x))
            family_total_seq[family] += frequencies_by_seq[gene]

            frequencies_by_clone[gene] = 0
            for x in str(count_clone).split(INT_SEP):
                frequencies_by_clone[gene] += int(float(x))
            family_total_clone[family] += frequencies_by_clone[gene]

        # calculate the frequency of each gene acording to the family
        for gene in frequencies_by_seq.keys():
            family = gene[:4]

            if family_total_seq[family] == 0:
                print('family_total_seq is zero for sample %s' % sample.study_title)
            else:
                frequencies_by_seq[gene] /= float(family_total_seq[family])

            if family_total_clone[family] == 0:
                print('family_total_clone is zero for sample %s' % sample.study_title)
            else:
                frequencies_by_clone[gene] /= float(family_total_clone[family])

        for gene in frequencies_by_seq.keys():
            # if we can't find an id for the gene in the file, look for a compound gene

            gene_id = None

            if gene in gene_ids_by_name:
                gene_id = gene_ids_by_name[gene]
            elif gene in compound_genes:
                gene_id = gene_ids_by_name[compound_genes[gene]]

            gd = GenesDistribution(
                frequency=frequencies_by_clone[gene],
                gene_id=gene_id,
                patient_id=sample.patient_id,
                sample_id=sample.id,
                count_by_clones=True,
            )
            session.add(gd)

            gd = GenesDistribution(
                frequency=frequencies_by_seq[gene],
                gene_id=gene_id,
                patient_id=sample.patient_id,
                sample_id=sample.id,
                count_by_clones=False,
            )
            session.add(gd)

    session.commit()
    return result


def calculate_patterns(session):
    # This function calculate the alleles that meets condition of pattern.
    result = ['Calculating patterns']

    # take all V genes
    genes = session.query(Gene).filter(Gene.name.like('%V%')).all()

    for gene in genes:
        alleles = session.query(Allele).filter(Allele.gene_id == gene.id).filter(Allele.name.notlike('%Del')).all()
        patterns = session.query(Allele).filter(Allele.gene_id == gene.id).filter(Allele.name.notlike('%Del')).filter(Allele.is_single_allele == 0).all()

        for pattern in patterns:
            pattern_seq = pattern.seq.replace("n","[a,g,c,t,n]")

            for allele in alleles:

                if (pattern.id != allele.id):
                    if (re.match(pattern_seq, allele.seq)):
                        # check if pattern already exists
                        if not session.query(AllelesPattern)\
                                .filter(AllelesPattern.allele_in_p_id == allele.id)\
                                .filter(AllelesPattern.pattern_id == pattern.id)\
                                .count():
                            ap = AllelesPattern(
                                allele_in_p_id=allele.id,
                                pattern_id=pattern.id
                            )
                            session.add(ap)
    session.commit()
    return result
