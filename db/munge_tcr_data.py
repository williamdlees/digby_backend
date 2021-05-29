# hack to move TCR files onto correct layout, naming etc

import yaml
import os
from os import listdir
from os.path import isfile, join
from shutil import copyfile
import csv

def munge(project, dir, yml_file):
    with open(os.path.join(dir, yml_file)) as fi:
        pdata = yaml.load(fi)

    os.mkdir(os.path.join('samples', project))

    for subject in pdata[project]['Samples'].keys():
        os.mkdir(os.path.join('samples', project, subject))

    geno_path = os.path.join(dir, 'genotypes')
    for file in listdir(geno_path):
        name_els = file.split('_')
        file = join(geno_path, file)
        if isfile(file) and project in file:
            sample = name_els[0] + '_' + name_els[1] + '_' + name_els[2]
            target = 'samples/' + project + '/' + sample + '/' + sample + '_genotype.tsv'
            convert_geno(file, target)

    haplo_path = os.path.join(dir, 'hapotypes')
    for file in listdir(haplo_path):
        name_els = file.split('_')
        file = join(haplo_path, file)
        if isfile(file) and project in file:
            sample = name_els[0] + '_' + name_els[1] + '_' + name_els[2]
            haplo_gene = name_els[4] + '_' + (name_els[5] if len(name_els) == 6 else '')
            haplo_gene = 'TRB' + haplo_gene.split('.')[0]
            if haplo_gene == 'TRBD2':
                haplo_gene = 'TRBD2_1_2'
            # fudge sample name to S1 in some projects
            if project != 'P19':
                sample = sample.replace('S2', 'S1')
                target = 'samples/' + project + '/' + sample + '/' + sample + '_gene-' + haplo_gene + '_haplotype.tsv'
                copyfile(file, target)


def convert_geno(source, target):
    with open(source, 'r') as fi, open(target, 'w', newline='') as fo:
        reader = csv.DictReader(fi, delimiter='\t')
        writer = None

        for row in reader:
            if writer is None:
                writer = csv.DictWriter(fo, delimiter='\t', fieldnames=reader.fieldnames)
                writer.writeheader()
            if len(row) > 1:
                if row['counts'] == 'NA' or row['GENOTYPED_ALLELES'] == 'Deletion':
                    row['Freq_by_Clone'] = 'NA'
                    row['Freq_by_Seq'] = 'NA'
                else:
                    la = len(row['GENOTYPED_ALLELES'].split(','))
                    row['Freq_by_Clone'] = ';'.join(row['counts'].split(',')[:la])
                    row['Freq_by_Seq'] = ';'.join(row['counts'].split(',')[:la])
                writer.writerow(row)


munge('P4', '../TCR_DS1', 'DS1.yml')
munge('P16', '../TCR_DS2', 'P16.yml')
munge('P17', '../TCR_DS2', 'P17.yml')
munge('P18', '../TCR_DS2', 'P18.yml')
munge('P19', '../TCR_DS3', 'P19.yml')


