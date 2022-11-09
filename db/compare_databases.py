# Simple comparisons of two databases

import argparse
from vdjbase_db import ContentProvider
from db.vdjbase_model import Allele, AllelesSample, Gene, SNP, HaplotypesFile, SamplesHaplotype
from db.vdjbase_airr_model import Sample


def main():
    parser = argparse.ArgumentParser(description='Compare two VDJbase databases')
    parser.add_argument('db1', help='First database')
    parser.add_argument('db2', help='Second database')
    args = parser.parse_args()

    db1 = ContentProvider(args.db1)
    db2 = ContentProvider(args.db2)

    alleles1 = db1.session.query(Allele.name, Allele.pipeline_name, Allele.similar, Allele.seq).all()
    alleles1 = {v[0]: {
        'pipeline_name': v[1].split(','),
        'similar': [x for x in v[2].split('|') if len(x) > 0],
        'seq': v[3]
    } for v in alleles1}

    alleles2 = db2.session.query(Allele.name, Allele.pipeline_name, Allele.similar, Allele.seq).all()
    alleles2 = {v[0]: {
        'pipeline_name': v[1].split(','),
        'similar': [x for x in v[2].split('|') if len(x) > 0],
        'seq': v[3]
    } for v in alleles2}

    allele_names1 = set(list(alleles1.keys()))
    allele_names2 = set(list(alleles2.keys()))

    if len(allele_names1 - allele_names2):
        print(f"{len(allele_names1 - allele_names2)} alleles in db1 are not in db2: {', '.join(sorted(list(allele_names1 - allele_names2)))}")

    if len(allele_names2 - allele_names1):
        print(f"{len(allele_names2 - allele_names1)} alleles in db2 are not in db1: {', '.join(sorted(list(allele_names2 - allele_names1)))}")

    if len(allele_names1 ^ allele_names2) == 0:
        print('The list of allele names in db1 is the same as that in db2.\n')

    joint_allele_names = allele_names1 & allele_names2

    diffs = False
    for allele_name in joint_allele_names:
        if alleles1[allele_name]['seq'] != alleles2[allele_name]['seq']:
            print(f"the sequence of {allele_name} in db1 does not match that in db2")
            diffs = True

    if not diffs:
        print('Sequences in all alleles listed in both databases agree')

    diffs = False
    for allele_name in joint_allele_names:
        sims1 = set(alleles1[allele_name]['similar'])
        sims2 = set(alleles2[allele_name]['similar'])
        if len(sims1 ^ sims2) > 0:
            print(f"the similar_alleles of {allele_name} in db1 do not match those in db2 ({', '.join(list(sims1))}, {', '.join(list(sims2))})")
            diffs = True

    if not diffs:
        print('Similars in all alleles listed in both databases agree')

    diffs = False
    for allele_name in joint_allele_names:
        pips1 = set(alleles1[allele_name]['pipeline_name'])
        pips2 = set(alleles2[allele_name]['pipeline_name'])
        if len(pips1 ^ pips2) > 0:
            print(f"the pipeline names of {allele_name} in db1 does not match those in db2 ({', '.join(list(pips1))}, {', '.join(list(pips2))})")
            diffs = True

    if not diffs:
        print('Pipeline names in all alleles listed in both databases agree')

    samples1 = db1.session.query(Sample).all()
    samples2 = db1.session.query(Sample).all()

    sample1_names = set([sample.sample_name for sample in samples1])
    sample2_names = set([sample.sample_name for sample in samples2])

    if len(sample1_names - sample2_names):
        print(f"{len(sample1_names - sample2_names)} samples in db1 are not in db2: {', '.join(sorted(list(sample1_names - sample2_names)))}")

    if len(sample2_names - sample1_names):
        print(f"{len(sample2_names - sample1_names)} samples in db2 are not in db1: {', '.join(sorted(list(sample2_names - sample1_names)))}")

    if len(sample1_names ^ sample2_names) == 0:
        print('The list of sample names in db1 is the same as that in db2.\n')

    sample_alleles1 = {
        sample.sample_name: set([x.allele.name for x in sample.alleles]) for sample in samples1
    }

    sample_alleles2 = {
        sample.sample_name: set([x.allele.name for x in sample.alleles]) for sample in samples2
    }

    joint_sample_names = list(sample1_names & sample2_names)

    diffs = False
    for sample_name in joint_sample_names:
        if len(sample_alleles1[sample_name] - sample_alleles2[sample_name]):
            print(f"{len(sample_alleles1[sample_name] - sample_alleles2[sample_name])} alleles in db1 sample {sample_name} are not listed in db2: {', '.join(sorted(list(sample_alleles1[sample_name] - sample_alleles2[sample_name])))}")
            diffs = True

        if len(sample_alleles2[sample_name] - sample_alleles1[sample_name]):
            print(f"{len(sample_alleles2[sample_name] - sample_alleles1[sample_name])} alleles in db1 sample {sample_name} are not listed in db2: {', '.join(sorted(list(sample_alleles2[sample_name] - sample_alleles1[sample_name])))}")
            diffs = True

    if not diffs:
        print('All samples present in both db1 and db2 have the same allele assignments')

if __name__ == '__main__':
    main()
