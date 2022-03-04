# TODO - I think we could do with some serious refactoring of the following functions
# the startnng point might be to have one place where the allele name we are working with
# - read from the genotype - is parsed into its components and other attributes (eg the
# 'ambiguous' alleles/genes it cites - are assembled.
import csv
import math
import os
import re
import pandas as pd
from sqlalchemy import func

from db.vdjbase_airr_model import Sample
from db.vdjbase_exceptions import DbCreationError
from db.vdjbase_formats import GENE_COLUMN, GENOTYPE_KDIFF_COLUMN, GENOTYPED_ALLELES_COLUMN, FREQ_BY_CLONE, INT_SEP
from db.vdjbase_model import Allele, AllelesSample, Gene, SNP, HaplotypesFile, SamplesHaplotype
from db.vdjbase_projects import compound_genes


def process_genotypes(ds_dir, species, dataset, session):
    result = ['Processing genotype files']
    samples = session.query(Sample).all()

    # Read names assigned to ambiguous alleles by the pipeline
    # pipeline_names carries the translation, which is used when reading the genotypes

    pipeline_names = {}   # allele name given pipeline name
    allele_names = {}     # pipeline name given allele name
    if os.path.isfile(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv')):
        with open(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv'), 'r') as fi:
            reader = csv.DictReader(fi)
            for row in reader:
                if len(row['PATTERN']) == 1:
                    row['PATTERN'] = '0' + row['PATTERN']   # thanks Excel
                pipeline_name = row['GENE'] + '*' + row['PATTERN']
                if len(row['ALLELES']) == 1:
                    allele = '0' + row['ALLELES']
                elif '.' in row['ALLELES']:
                    allele_list = row['ALLELES'].replace(' ', '').split(',')
                    for a in allele_list:
                        if '.' not in a:
                            result.append('Error in ambiguous allele definition for %s: each item must have a dot.' % row['GENE'])
                            break
                    gene_num = allele_list[0].split('.')[0]
                    row['GENE'] = row['GENE'].split('-')[0] + '-' + gene_num
                    alleles = []
                    for item in allele_list:
                        a, n = item.split('.')
                        if a != gene_num:
                            alleles.append('%s.%s' % (a, n))
                        else:
                            alleles.append(n)
                    allele = '_'.join(alleles)
                elif ',' in row['ALLELES'] or ' ' in row['ALLELES']:
                    allele = row['ALLELES'].replace(' ', '').replace(',', '_')
                else:
                    allele = row['ALLELES']

                allele_name = row['GENE'] + '*' + allele
                allele_names[allele_name] = pipeline_name
                pipeline_names[pipeline_name] = allele_name

                gene_name = allele_name.split('*')[0]
                pipeline_gene_name = pipeline_name.split('*')[0]

                if gene_name != pipeline_gene_name:
                    add_compound_gene(session, allele_name, pipeline_gene_name)

    for sample in samples:
        genotype_file = os.path.join('samples', sample.study.study_title, sample.patient.study_title, sample.study_title + '_geno_H_binom.tab').replace('\\', '/')

        if not os.path.isfile(os.path.join(ds_dir, genotype_file)):
            genotype_file = os.path.join('samples', sample.study.study_title, sample.study_title, sample.study_title + '_genotype.tsv').replace('\\', '/')  # new directory layout

        if not os.path.isfile(os.path.join(ds_dir, genotype_file)):
            genotype_file = os.path.join('samples', sample.study.study_title, sample.study_title, sample.study_title + '.tsv').replace('\\', '/')  # another new directory layout

        if os.path.isfile(os.path.join(ds_dir, genotype_file)):
            sample.genotype = genotype_file
            sample_genotype(os.path.join(ds_dir, genotype_file), sample.id, sample.patient.id, pipeline_names, allele_names, session)
            session.commit()
        else:
            result.append('Error: no genotype file for sample %s' % sample.study_title)

    return result


allele_pattern = re.compile("^[0-9]{2}$")
mut_pattern = re.compile("^[A,G,T,C,a,g,t,c][0-9]+[A,G,T,C,a,g,t,c]$")
ext_mut_pattern = re.compile("^[0-9]+[A,G,T,C,a,g,t,c]+[0-9]+$")
gene_pattern = re.compile(".+\.[0-9]{2}$")
nt_pattern = re.compile("[A,G,T,C,a,g,t,c]+")


def sample_genotype(inputfile, sample_id, patient_id, pipeline_names, allele_names, session):
    """
    Upload genotype for sample.
    :param inputfile: the path to the genotype
    :param sample_id: sample id
    """
    print(inputfile)
    haplo = "geno"
    genotype = pd.read_csv(inputfile, sep="\t")

    for index, row in genotype.iterrows():
        gene = row[GENE_COLUMN]
        kdiff = row[GENOTYPE_KDIFF_COLUMN]

        for index, allele in enumerate(row[GENOTYPED_ALLELES_COLUMN].split(",")):
            if (allele == "Unk") | ("NR" in allele):
                continue
            elif len(str(allele)) == 1:
                allele = "0" + str(allele)
            elif ("del" in allele.lower()):
                allele = "Del"


            # check if the allele exists in the genotype according to the clone size
            if allele != "Del":
                # check for bug found in TRB genotypes
                if len(str(row[FREQ_BY_CLONE]).split(INT_SEP)) <= index:
                    print('Error: FREQ_BY_CLONE for gene %s does not have enough values' % gene)
                    continue
                clone_size = int(float(str(row[FREQ_BY_CLONE]).split(INT_SEP)[index]))
                if not clone_size:
                    continue

            # parse allele, handling old ambiguous style (01_02) and new bp style
            ambig_alleles = []
            allele_snps = []
            ext_snps = []
            base_allele = allele

            if '_' in allele:

                base_allele = allele.split('_')[0]
                rep = allele.split("_")[1:]

                for r in rep:
                    if re.search(allele_pattern, r):
                        ambig_alleles.append(r)
                    elif re.search(mut_pattern, r):
                        allele_snps.append(r)
                    elif re.search(ext_mut_pattern, r):
                        ext_snps.append(r)
                    else:
                        base_allele += '_' + r     # this is a piece that isn't part of the tigger naming scheme and must therefore
                                                        # be part of the base allele name

            base_allele_name = gene + "*" + base_allele
            pipeline_name = gene + "*" + allele

            if base_allele_name in pipeline_names:
                pipeline_name = base_allele_name
                base_allele_name = pipeline_names[base_allele_name]
            elif 'bp' in base_allele_name:        # will only happen if one of the 'bp' alleles has not been put in ambiguous_allele_names.csv
                # fudge for time being
                pipeline_name = base_allele_name
                base_allele_name = base_allele_name.replace('bp', '')
                pipeline_names[pipeline_name] = base_allele_name
                print('%s is not listed in ambiguous_allele_names.csv: assuming it corresponds to %s' % (pipeline_name, base_allele_name))
            elif 'ap' in base_allele_name:        # will only happen if one of the 'bp' alleles has not been put in ambiguous_allele_names.csv
                # fudge for time being
                pipeline_name = base_allele_name

                if '_' in base_allele:
                    base_allele_name = gene + "*" + base_allele.split('_')[0]

                if base_allele_name in pipeline_names:
                    base_allele_name = pipeline_names[base_allele_name]
                else:
                    base_allele_name = base_allele_name.replace('ap', '')

                    if base_allele_name in pipeline_names:
                        base_allele_name = pipeline_names[base_allele_name]
                    else:
                        pipeline_names[pipeline_name] = base_allele_name
                print('%s is not listed in ambiguous_allele_names.csv: assuming it corresponds to %s' % (pipeline_name, base_allele_name))

            this_allele_name = base_allele_name

            if len(ambig_alleles) > 0:
                this_allele_name += '_' + '_'.join(ambig_alleles)
                if len(pipeline_name) > 0:
                    pipeline_name += '_' + '_'.join(ambig_alleles)

            if len(allele_snps) > 0:
                this_allele_name += '_' + '_'.join(allele_snps)
                if len(pipeline_name) > 0:
                    pipeline_name += '_' + '_'.join(allele_snps)

            if len(ext_snps) > 0:
                this_allele_name += '_' + '_'.join(ext_snps)
                if len(pipeline_name) > 0:
                    pipeline_name += '_' + '_'.join(ext_snps)

            add2sample(this_allele_name, base_allele_name, sample_id, haplo, patient_id, kdiff, pipeline_name, allele_snps, session)


# Each unique sequence is only present on one single row in Allele. If multiple allele names
# correspond to the same sequence, they are listed in that row in the 'similar' field.
# This function returns that row, given the allele name

def find_allele_or_similar(allele_name, session):
    try:
        allele = session.query(Allele).filter(Allele.name == allele_name).one_or_none()
    except Exception as e:
        print('Multiple rows found for allele %s: check fasta files in reference directory.' % allele_name)

    if allele is not None:
        return allele

    allele = session.query(Allele).filter(Allele.similar.ilike(("%|" + allele_name + "|%"))).one_or_none()
    return allele


# Add a row to AllelesSample reflecting the presence of this allele in the sample.
# If the allele is not present in Allele already, add it there
# if the allele has a 'pipeline name', translate it to the 01_02 form

def add2sample (allele_name, base_allele_name, sample_id, haplo, pid, kdiff, pipeline_name, allele_snps, session):
    kdiff = float(kdiff)
    if math.isnan(kdiff):
        kdiff = 0.0

    alleles_sample = AllelesSample(
        hap=haplo,
        kdiff=kdiff,
        sample_id=sample_id,
        patient_id=pid,
    )

    allele = find_allele_or_similar(allele_name, session)
    if allele is None:
        allele = new_allele(allele_name, base_allele_name, pipeline_name, allele_snps, session)

    if len(pipeline_name) > 0:
        if len(allele.pipeline_name) == 0:
            allele.pipeline_name = pipeline_name
        elif pipeline_name not in allele.pipeline_name:
            allele.pipeline_name = allele.pipeline_name + ', ' + pipeline_name

    alleles_sample.allele_id = allele.id

    asc = session.query(AllelesSample)\
        .filter(AllelesSample.allele_id == alleles_sample.allele_id)\
        .filter(AllelesSample.sample_id == alleles_sample.sample_id)\
        .filter(AllelesSample.hap == haplo).count()

    if asc == 0:
        session.add(alleles_sample)

# allele_name is the name of the allele as it will appear in VDJbase.
# base_allele_name is the name without snps. It can include ambiguous alleles.
# pipeline_name is the name of the allele as found in the genotype file
def new_allele(allele_name, base_allele_name, pipeline_name, allele_snps, session):
    ambiguous_alleles = [allele_name.split("*")[1].split("_")[0]]

    # check that we can find the first allele mentioned in allele_name. This may have a format
    # that's a little more complex than a simple number, for example IGHV4-61*01_S9382, but will
    # precede any snps or ambiguous allele references

    reps = allele_name.split('*')[1].split('_')
    allele_parts = [reps[0]]
    reps = reps[1:]

    for r in reps:
        if re.search(allele_pattern, r) or re.search(gene_pattern, r) or re.search(mut_pattern, r) or re.search(ext_mut_pattern, r):
            break

        allele_parts.append(r)

    first_allele_name = allele_name.split('*')[0] + '*' + '_'.join(allele_parts)
    first_allele = find_allele_or_similar(first_allele_name, session)

    if first_allele is None:
        # check for duplicate (D) allele
        allele_D_name = first_allele_name.split("*")[0] + "D*" + first_allele_name.split("*")[1].split("_")[0]
        first_allele = find_allele_or_similar(allele_D_name, session)

    if first_allele is None:
        raise DbCreationError('Error processing allele %s: base allele not in reference set' % first_allele_name)

    # Check for a compound gene in the pipeline name, set gene id accordingly
    base_gene = pipeline_name.split('*')[0]
    if base_gene in compound_genes:
        gene_id = session.query(Gene.id).filter(Gene.name == compound_genes[base_gene]).one_or_none()[0]
    else:
        gene_id = first_allele.gene_id

    seq = first_allele.seq
    seq = "".join(seq.split())

    # For each gene/allele referenced in the name as being indistinguishable, check that
    # we have the sequence in the reference set. Construct a merged 'ambiguous' sequence
    # which masks out any sites differing between all the mentioned alleles with ns. Also
    # take account in the sequence of any SNPs in the name.
    # Example: TRBV3-1*01_02_2.01_2.02_2.03_a234g. Here there are 4 alleles to check and
    # merge, plus one mutation. Check for syntax errors in the allele designation. War
    # (and ignore) if any incorrectly-formatted terms are encountered.

    final_allele_name = base_allele_name
    if allele_name != base_allele_name:
        rep = allele_name.replace(base_allele_name + '_', '').lower().split('_')
        is_novel_allele = False
        for r in rep:
            if re.search(allele_pattern, r) or re.search(gene_pattern, r):
                if r not in ambiguous_alleles:
                    final_allele_name += "_" + r
                    allele_seq_1 = seq

                    if '.' not in r:
                        base_name = allele_name.split("*")[0]
                        a_name = r
                    else:
                        base_name = allele_name.split('-')[0] + '-' + r.split('.')[0]
                        a_name = r.split('.')[1]

                    allele2_name = base_name + "*" + a_name
                    allele2 = find_allele_or_similar(allele2_name, session)

                    if allele2 is None:
                        allele2 = find_allele_or_similar(base_name + "D*" + a_name, session)

                    if allele2 is None:
                        raise DbCreationError('Error processing allele %s: base allele %s not in reference set' % (allele_name, allele2_name))

                    allele_seq_2 = allele2.seq
                    seq = ""
                    length_diff = len(allele_seq_1) - len(allele_seq_2)

                    if length_diff > 0:
                        allele_seq_1 += "n"*length_diff
                    elif length_diff < 0:
                        allele_seq_2 += "n" * (length_diff * -1)

                    for nuc1, nuc2 in zip(allele_seq_1, allele_seq_2):
                        if nuc1 == nuc2:
                            seq += nuc1
                        else:
                            seq += "n"

                ambiguous_alleles.append(r)

            elif re.search(mut_pattern, r):
                final_allele_name += "_" + r
                place = int(r[1:][:len(r) - 2]) - 1
                seq = seq[:place] + r[len(r) - 1] + seq[place + 1:]
                is_novel_allele = True

            elif re.search(ext_mut_pattern, r):
                final_allele_name += '_' + r
                nts = re.search(nt_pattern, r)
                first = int(r[:nts.span()[0]])
                last = int(r[nts.span()[1]:])
                new_seq = seq[:first - 1] + r[nts.span()[0]:nts.span()[1]]
                if len(seq) > len(new_seq):
                    new_seq += seq[len(new_seq) - len(seq):]
                seq = new_seq
                is_novel_allele = True
            else:
                print(f'Error processing allele {allele_name}: unrecognised term in allele designation: {r} was ignored')

    same_seq_allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()    # we only expect one - one row in Allele for each sequence

    if same_seq_allele is not None:
        allele = None
        temp = same_seq_allele.study_title
        sim = same_seq_allele.similar
        if (temp == final_allele_name):
            allele = same_seq_allele
        else:
            # if we're adding a gene-ambiguous allele, it should become the allele name

            if '.' in final_allele_name.split('*')[1]:
                final_allele_name, temp = temp, final_allele_name

            if sim is not None and len(sim) > 0:
                if (final_allele_name in sim):
                    allele = same_seq_allele
                else:
                    sim = sim + ", |" + final_allele_name + "|"
            else:
                sim = "|" + final_allele_name + "|"

        if not allele:
            allele_table = Allele.__table__
            stmt = allele_table.update().where(allele_table.c.seq == seq).values(similar=sim, name = temp)
            session.execute(stmt)
            session.commit()
            allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()

    else:
        if pipeline_name == '':
            breakpoint()

        allele = Allele(
            name=final_allele_name,
            seq=seq,
            seq_len=str(len(seq)),
            gene_id=gene_id,
            is_single_allele=len(ambiguous_alleles) == 1,
            appears=0,
            low_confidence=False,
            novel=is_novel_allele,
            max_kdiff=0.0,
            similar='',
            pipeline_name=pipeline_name,
            closest_ref = first_allele,
        )

        session.add(allele)
        session.flush()
        session.refresh(allele)

        # add the snps

        for snp in allele_snps:
            s = SNP(
                from_base=snp[0],
                pos=int(snp[1:-1]),
                to_base=snp[-1],
                allele=allele
            )
            session.add(s)

    return allele


# Construct and add a 'compound gene', e.g. TRBV5/6, based on an allele name, eg TRBV6-5*01_6.01_6.02_6.03
def add_compound_gene(session, vdjbase_allele_name, pipeline_gene_name):
    if pipeline_gene_name in compound_genes:
        return

    # construct the compound gene from the allele extensions

    root = vdjbase_allele_name.split('-')[0]
    exts = vdjbase_allele_name.split('_')[1:]
    nums = []
    num = vdjbase_allele_name.split('-')[1]
    num = num.split('*')[0]
    nums.append(num)
    for ext in exts:
        if '.' in ext:
            num = ext.split('.')[0]
            if num not in nums:
                nums.append(num)

    vdjbase_gene_name = root + '-' + '/'.join(nums)
    compound_genes[pipeline_gene_name] = vdjbase_gene_name

    max_locus_order = session.query(func.max(Gene.locus_order)).one_or_none()[0]
    max_alpha_order = session.query(func.max(Gene.locus_order)).one_or_none()[0]
    species = session.query(Gene.species).filter(Gene.locus_order == max_locus_order).one_or_none()[0]

    print('Adding %s' % vdjbase_gene_name)

    g = Gene(
        name=vdjbase_gene_name,
        type=vdjbase_gene_name[:4],
        family=vdjbase_gene_name.split('-')[0],
        species=species,
        locus_order=max_locus_order+1,
        alpha_order=max_alpha_order+1,
        pseudo_gene=False
    )
    session.add(g)
    session.commit()


# Add Allele records for *Del
def add_deleted_alleles(session):
    result = ['Adding deleted alleles']
    genes = session.query(Gene).all()

    for gene in genes:
        a = Allele(
            name=gene.study_title + '*Del',
            seq='',
            seq_len='0',
            gene_id=gene.id,
            is_single_allele=False,
            appears=0,
            low_confidence=False,
            novel=0,
            max_kdiff=0.0,
            similar='',
            pipeline_name='',
        )
        session.add(a)
    session.commit()
    return result


def process_haplotypes_and_stats(ds_dir, species, dataset, session):
    result = ['Processing haplotype files']
    samples = session.query(Sample).all()

    for sample in samples:
        sample_dir = os.path.join('samples', sample.study.study_title, sample.patient.study_title) #old format

        if not os.path.isdir(os.path.join(ds_dir, sample_dir)):
            sample_dir = os.path.join('samples', sample.study.study_title, sample.study_title) #new format

        if os.path.isdir(os.path.join(ds_dir, sample_dir)):
            for filename in os.listdir(os.path.join(ds_dir, sample_dir)):
                if sample.study_title in filename:
                    if 'haplotype.' in filename:
                        haplo_gene = filename.replace('_haplotype.tab', '')
                        haplo_gene = haplo_gene.replace('_haplotype.tsv', '')
                        haplo_gene = haplo_gene.split('_gene-')[1]
                        process_haplotype(os.path.join(sample_dir, filename).replace('\\', '/'), sample, haplo_gene, session)
                    elif 'ogrdb_plots' in filename:
                        sample.genotype_report = os.path.join(sample_dir, filename).replace('\\', '/')
                    elif 'ogrdb_report' in filename:
                        sample.genotype_stats = os.path.join(sample_dir, filename).replace('\\', '/')

            if sample.genotype_report is None:
                print("No genotype report for sample %s" % sample.study_title)
            if sample.genotype_stats is None:
                print("No genotype stats for sample %s" % sample.study_title)
        else:
            print("No sample directory for sample %s" % sample.study_title)

    session.flush()
    return result


def process_haplotype(filename, sample, haplo_gene, session):
    gene = "-".join(haplo_gene.split("-")[:-1])
    alleles = haplo_gene.split("-")[-1]
    allele1 = gene + "_" + alleles.split("_")[0]
    allele2 = gene + "_" + alleles.split("_")[1]

    hf = HaplotypesFile(
        by_gene=haplo_gene,
        allele_col1=allele1,
        allele_col2=allele2,
        file=filename,
    )
    session.add(hf)
    session.flush()

    sh = SamplesHaplotype(
        samples_id=sample.id,
        haplotypes_file_id=hf.id,
    )
    session.add(sh)

