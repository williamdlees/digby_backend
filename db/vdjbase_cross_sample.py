#
# Cross-sample processing - once sample updates are complete
#
import sqlite3

import pandas as pd
import os, math

from sqlalchemy import func
from db.vdjbase_model import Allele, AllelesSample, Gene, GenesDistribution, AllelesPattern
from db.vdjbase_airr_model import Sample, Patient
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
    gene_ids_by_name = dict(zip([gene.name for gene in genes], [gene.id for gene in genes]))

    freqs = session.query(Sample.id, Patient.id, func.sum(AllelesSample.freq_by_seq), func.sum(AllelesSample.freq_by_clone), Gene.name, Sample.sample_name) \
        .filter(Patient.id == Sample.patient_id) \
        .filter(AllelesSample.sample_id == Sample.id) \
        .filter(AllelesSample.allele_id == Allele.id) \
        .filter(Allele.gene_id == Gene.id) \
        .group_by(*(Sample.id, Gene.id)) \
        .all()

    for sample_id, sample_name in list(set([(x[0], x[5]) for x in freqs])):
        frequencies_by_seq = {}
        frequencies_by_clone = {}
        family_total_seq = {}
        family_total_clone = {}

        for _, patient_id, count_seq, count_clone, gene_name, _ in [x for x in freqs if x[0] == sample_id]:
            family = gene_name[:4]
            if family not in family_total_seq.keys():
                family_total_seq[family] = 0
                family_total_clone[family] = 0

            frequencies_by_seq[gene_name] = count_seq
            frequencies_by_clone[gene_name] = count_clone
            family_total_seq[family] += count_seq
            family_total_clone[family] += count_clone


        # calculate the frequency of each gene acording to the family
        for gene in frequencies_by_seq.keys():
            family = gene[:4]

            if family_total_seq[family] == 0:
                print('family_total_seq for family %s is zero for sample %s' % (family, sample_name))
            else:
                frequencies_by_seq[gene] /= float(family_total_seq[family])

            if family_total_clone[family] == 0:
                print('family_total_clone for family %s is zero for sample %s' % (family, sample_name))
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
                patient_id=patient_id,
                sample_id=sample_id,
                count_by_clones=True,
            )
            session.add(gd)

            gd = GenesDistribution(
                frequency=frequencies_by_seq[gene],
                gene_id=gene_id,
                patient_id=patient_id,
                sample_id=sample_id,
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
