# hack to move TCR files onto correct layout, naming etc

import yaml
import os
from os import listdir
from os.path import isfile, join
from shutil import copyfile

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
            copyfile(file, target)

    haplo_path = os.path.join(dir, 'haplotypes')
    for file in listdir(haplo_path):
        name_els = file.replace('_MC1', '').split('_')      # remove MC1 from some names in P20 June run
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


munge('P4', '../TCR_June_2021_2/DS1', 'DS1.yml')
munge('P16', '../TCR_June_2021_2/DS2', 'P16.yml')
munge('P17', '../TCR_June_2021_2/DS2', 'P17.yml')
munge('P18', '../TCR_June_2021_2/DS2', 'P18.yml')
munge('P19', '../TCR_June_2021_2/DS3', 'P19.yml')
munge('P20', '../TCR_June_2021_2/DS4', 'P20.yml')


