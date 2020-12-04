# Manage novel alleles assigned by vdjbase

import csv
import os.path
from datetime import datetime


def check_for_dataset(novels, species, dataset):
    if species not in novels:
        novels[species] = {}
    if dataset not in novels[species]:
        novels[species][dataset] = []


def load_novel_alleles(file):
    novels = {}

    if not os.path.isfile(file):
        print('Warning: novel file %s not found: creating new novel allele database' % file)
        return novels

    with open(file, 'r') as fi:
        reader = csv.DictReader(fi)
        for row in reader:
            check_for_dataset(novels, row['species'], row['dataset'])
            novels[row['species']][row['dataset']].append(row)

    return novels


def save_novel_alleles(novels, file):
    headers = ['species', 'dataset', 'name', 'index', 'sequence', 'origin', 'date_assigned']

    with open(file, 'w', newline='') as fo:
        writer = csv.DictWriter(fo, fieldnames=headers)
        writer.writeheader()
        for species in novels.keys():
            for dataset in novels[species].keys():
                for row in novels[species][dataset]:
                    writer.writerow(row)


def find_novel_allele(novels, species, dataset, seq):
    check_for_dataset(novels, species, dataset)
    for row in novels[species][dataset]:
        if row['sequence'] == seq:
            return row

    return None


def find_next_index(novels, species, dataset):
    check_for_dataset(novels, species, dataset)
    max_index = 0
    for row in novels[species][dataset]:
        max_index = max(max_index, row['index'])
    return max_index + 1


def assign_novel_gene(novels, species, dataset, seq, prefix, family, origin):
    allele = find_novel_allele(novels, species, dataset, seq)

    if allele is not None:
        return allele

    index = find_next_index(novels, species, dataset)
    allele = {
        'species': species,
        'dataset': dataset,
        'name': '%s%s.%d' % (prefix, family, index),
        'index': index,
        'sequence': seq,
        'origin': origin,
        'date_assigned': datetime.now().isoformat(timespec='seconds')
    }

    novels[species][dataset].append(allele)
    return allele


