# TODO - I think we could do with some serious refactoring of the following functions
# the starting point might be to have one place where the allele name we are working with
# - read from the genotype - is parsed into its components and other attributes (eg the
# 'ambiguous' alleles/genes it cites - are assembled.
import csv
import math
import os
import re
import pandas as pd
from sqlalchemy import func
import glob
from collections import namedtuple
from receptor_utils import simple_bio_seq as simple

from db.vdjbase_airr_model import Sample
from db.vdjbase_exceptions import DbCreationError
from db.vdjbase_model import Allele, AllelesSample, Gene, SNP, HaplotypesFile, SamplesHaplotype
from db.vdjbase_projects import compound_genes

# transaction log of alleles created
new_alleles = {}


def process_genotypes(ds_dir, species, dataset, session):
    """
    Process genotype files for samples in the given directory.

    This function reads genotype files for each sample in the given directory, finds the assigned ambiguous alleles,
    and parses the genotypes.

    It logs allele audit information and returns a list of log messages.

    :param ds_dir: The directory containing the genotype files.
    :type ds_dir: str
    :param species: The species of the samples.
    :type species: str
    :param dataset: The dataset (locus) of the samples.
    :type dataset: str
    :param session: The database session to use.
    :type session: sqlalchemy.orm.session.Session
    :return: A list of log messages.
    :rtype: list of str
    """

    result = ['Processing genotype files']
    samples = session.query(Sample).all()

    # Read names assigned to ambiguous alleles by the pipeline
    # pipeline_names carries the translation, which is used when reading the genotypes

    allele_names, pipeline_names = read_ambiguous_alleles_file(ds_dir, result, session)

    for sample in list(samples):
        sample.genotype, sample.asc_genotype = find_genotype_files(ds_dir, sample, result)

        if sample.asc_genotype or sample.genotype:
            sample_genotype(sample, pipeline_names, session)
        else:
            result.append('Error: no genotype file for sample %s - removing from sample list' % sample.sample_name)
            session.delete(sample)

        # fix up filenames for the live database
        if sample.genotype:
            sample.genotype = sample.genotype.replace('\\', '/').replace(ds_dir + '/', '')
        if sample.asc_genotype:
            sample.asc_genotype = sample.asc_genotype.replace('\\', '/').replace(ds_dir + '/', '')
            
    session.commit()

    # dump audit log

    with open('allele_audit_log.csv', 'w') as fo:
        fo.write('allele_name,pipeline_names,similar_names\n')
        for new_allele_name in sorted(new_alleles.keys()):
            rec = new_alleles[new_allele_name]
            fo.write(f"{rec['allele_name']}, {';'.join(sorted(rec['pipeline_names']))}, {';'.join(sorted(rec['similar']))}\n")

    return result


def find_genotype_files(ds_dir, sample, result):
    """
    Finds genotype files for the given sample in the specified directory, and returns the file paths for the tigger and asc files.

    :param ds_dir: The directory where the sample genotype files are located.
    :type ds_dir: str
    :param sample: The sample for which to find the genotype files.
    :type sample: Sample object
    :param result: A list to which any warning messages will be appended.
    :type result: list
    :return: A tuple containing the file paths for the tigger and asc genotype files, respectively. If no file is found, the corresponding variable in the tuple will be None.
    :rtype: tuple

    This function searches for genotype files for the given sample in the 'samples' directory of the specified ds_dir.
    It first looks for files with 'geno' in their name. If a file is found with 'asc' in its name, it is considered as the asc genotype file,
    and all other geno files are considered as tigger genotype files. If no file is found with 'asc' in its name, then the first geno file
    found is considered as the tigger genotype file.

    If multiple files of either type are found, a warning message will be appended to the result list, and the function will choose the
    first file.
    """
    geno_files = glob.glob(os.path.join(ds_dir, 'samples', sample.study.study_name, sample.sample_name, '*_geno*'))
    asc_files = [x for x in geno_files if 'asc' in os.path.basename(x)]
    tigger_files = [x for x in geno_files if 'asc' not in os.path.basename(x)]

    asc_file = None
    tigger_file = None

    if asc_files:
        if len(asc_files) > 1:
            result.append('Multiple asc genotype files found. Picking first.')
        asc_file = asc_files[0]

    if tigger_files:
        if len(tigger_files) > 1:
            result.append('Multiple tigger genotype files found. Picking first.')
        tigger_file = tigger_files[0]

    return tigger_file, asc_file


def read_ambiguous_alleles_file(ds_dir, result, session):
    """
    Reads an ambiguous allele definition file from the given directory, creates two dictionaries to map allele names and pipeline names,
    and updates the database session with compound genes if necessary.

    :param ds_dir: The directory where the ambiguous allele definition file is located.
    :type ds_dir: str
    :param result: A list to which any error messages will be appended.
    :type result: list
    :param session: A session object to update the database with compound genes.
    :type session: database session
    :return: A tuple containing two dictionaries: allele_names and pipeline_names.
    :rtype: tuple

    The allele_names dictionary maps allele names to pipeline names, while the pipeline_names dictionary maps pipeline names to allele names.

    Each row in the ambiguous allele definition file should have the following columns:
    - GENE: The gene name.
    - PATTERN: The pattern of the gene.
    - ALLELES: A comma-separated list of allele names for the gene.

    If the length of the PATTERN or ALLELES field is 1, a '0' is added before the value. If the ALLELES field contains dots, it is split into a list of allele names,
    and the gene number is extracted. If the alleles belong to a different gene, the gene name is updated by replacing the original gene number with the new gene number.

    If an error is encountered during the parsing of the ambiguous allele definition file, it will be appended to the result list.

    If the gene name in the allele_name variable does not match the gene name in the pipeline_name variable, the add_compound_gene function will be called to
    add a compound gene to the database session.
    """
    pipeline_names = {}  # allele name given pipeline name
    allele_names = {}  # pipeline name given allele name
    if os.path.isfile(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv')):
        with open(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv'), 'r') as fi:
            reader = csv.DictReader(fi)
            for row in reader:
                if len(row['PATTERN']) == 1:
                    row['PATTERN'] = '0' + row['PATTERN']  # thanks Excel
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
    return allele_names, pipeline_names


allele_pattern = re.compile("^[0-9]{2}$")
mut_pattern = re.compile("^[A,G,T,C,a,g,t,c][0-9]+[A,G,T,C,a,g,t,c]$")
ext_mut_pattern = re.compile("^[0-9]+[A,G,T,C,a,g,t,c]+[0-9]+$")
gene_pattern = re.compile(".+\.[0-9]{2}$")
nt_pattern = re.compile("[A,G,T,C,a,g,t,c]+")


def parse_allele(allele):
    """
    Parse allele names, handling old ambiguous style (01_02) and new bp style.

    Note that ambiguity is with respect to the sequence *before* SNPs are taken into consideration,
    i.e. the SNPs apply to all the alleles listed as ambiguous. If different SNPs were allowed
    to be associated with individual unmutated alleles, every allele in the reference set would
    need to be included in the name.

    :param allele: The allele name to be parsed.
    :type allele: str
    :return: A tuple containing:
        - base_allele: the base allele name, without any ambiguity or SNP information.
        - ambig_alleles: a list of ambiguous alleles, if any are present in the allele name.
        - allele_snps: a list of SNPs associated with the allele, if any are present in the allele name.
        - ext_snps: a list of external SNPs associated with the allele, if any are present in the allele name.
    :rtype: tuple[str, list[str], list[str], list[str]]
    """
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
            elif re.search(gene_pattern, r):
                ambig_alleles.append(r)
            elif re.search(mut_pattern, r):
                allele_snps.append(r)
            elif re.search(ext_mut_pattern, r):
                ext_snps.append(r)
            else:
                base_allele += '_' + r  # this is a piece that isn't part of the tigger naming scheme and must therefore
                # be part of the base allele name
    return base_allele, ambig_alleles, allele_snps, ext_snps


def parse_ambiguous_allele(allele, gene, pipeline_names):
    """
    Parse an allele, taking into account ambiguous alleles listed in pipeline_names.

    :param allele: The name of the allele.
    :param gene: The name of the gene.
    :param pipeline_names: A dictionary of ambiguous allele names.

    :returns: A tuple containing the parsed allele SNPs, base allele name, pipeline name, and this allele name.
    """
    base_allele, ambig_alleles, allele_snps, ext_snps = parse_allele(allele)
    base_allele_name = gene + "*" + base_allele
    pipeline_name = gene + "*" + allele

    if base_allele_name in pipeline_names:
        pipeline_name = base_allele_name
        base_allele_name = pipeline_names[base_allele_name]

        # If the base allele might be ambiguous, parse it out. But we don't expect SNPs in a pipeline name lookup,
        # they stay as they were in the original pipeline name
        if '_' in base_allele_name:
            base_allele, ambig_alleles, _, _ = parse_allele(base_allele_name.split('*')[1])
            base_allele_name = base_allele_name.split('*')[0] + "*" + base_allele

    this_allele_name = base_allele_name

    if len(ambig_alleles) > 0:
        this_allele_name += '_' + '_'.join(ambig_alleles)

    if len(allele_snps) > 0:
        this_allele_name += '_' + '_'.join(allele_snps)
        if len(pipeline_name) > 0 and '_'.join(allele_snps) not in pipeline_name:
            pipeline_name += '_' + '_'.join(allele_snps)

    if len(ext_snps) > 0:
        this_allele_name += '_' + '_'.join(ext_snps)
        if len(pipeline_name) > 0 and '_'.join(ext_snps) not in pipeline_name:
            pipeline_name += '_' + '_'.join(ext_snps)

    return allele_snps, base_allele_name, pipeline_name, this_allele_name


# Read a genotype. Add the contents of each row to the database
def sample_genotype(sample, pipeline_names, session):
    """
    Read a genotype file for the given sample and add the contents of each row to the database.

    :param sample: The sample object for which to process the genotype.
    :type sample: Sample

    :param pipeline_names: A dictionary of ambiguous allele names and their pipeline-assigned names.
    :type pipeline_names: dict

    :param session: The database session object to use for database operations.
    :type session: sqlalchemy.orm.session.Session

    :return: A list of strings containing information about the processing of the genotype file.
    :rtype: list
    """
    processed_gene_types = []

    if sample.asc_genotype:
        processed_gene_types = process_asc_genotype(sample, processed_gene_types, pipeline_names, session)

    process_tigger_genotype(sample, processed_gene_types, pipeline_names, session)


def process_asc_genotype(sample, processed_gene_types, pipeline_names, session):
    """
    Processes an ASC genotype file for a sample, adds the resulting data to the database, and returns a list of
    processed gene types.

    :param sample: A `Sample` object representing the sample to process.
    :type sample: Sample

    :param processed_gene_types: A list of gene types that have already been processed.
    :type processed_gene_types: list

    :param pipeline_names: A list of pipeline names to use when parsing ambiguous allele names.
    :type pipeline_names: list

    :param session: The database session to use for adding the processed data.
    :type session: Session

    :return: The updated list of processed gene types, including any new gene types processed during this function call.
    :rtype: list
    """
    my_processed_types = []
    print(sample.asc_genotype)
    genotype = pd.read_csv(sample.asc_genotype, sep='\t')
    for index, row in genotype.iterrows():
        gene = row["gene"]
        gene_type = gene[3]

        if gene_type in processed_gene_types:
            continue

        if gene_type not in my_processed_types:
            my_processed_types.append(gene_type)

        kdiff = 0

        allele_counts = {}
        ac = str(row["counts"]).split(',')
        total_count = sum([int(x) if x.isnumeric() else 0 for x in ac])
        absolute_fraction = {}
        af = str(row["absolute_fraction"]).split(',')

        for index, allele in enumerate(str(row["imgt_alleles"]).split(",")):
            if not ((allele == "Unk") or ("NR" in allele) or ("del" in allele.lower())) and index < len(ac):
                allele_counts[allele] = int(ac[index])
                absolute_fraction[allele] = float(af[index])

        for index, allele in enumerate(str(row["genotyped_imgt_alleles"]).split(",")):
            if (allele == "Unk") or ("NR" in allele):
                continue
            elif len(str(allele)) == 1:
                allele = "0" + str(allele)
            elif ("del" in allele.lower()):
                allele = "Del"

            count = allele_counts[allele] if allele in allele_counts else 0
            freq_by_clone = count
            freq_by_seq = count

            # check for multiple (ambiguous) allele names in the call and convert format to VDJbase if found
            # a bit of a weird split, as / can be in the name, as well as the separators between names

            cell_type = allele[:2]
            allele_calls = allele.split('/' + cell_type)

            if len(allele_calls) > 1:
                for i in range(1, len(allele_calls)):
                    allele_calls[i] = cell_type + allele_calls[i]

                allele = asc_name_to_vdjbase(allele_calls, session)

            allele_part = '*'.join(allele.split('*')[1:])
            gene_part = allele.split('*')[0]
            allele_snps, base_allele_name, pipeline_name, this_allele_name = parse_ambiguous_allele(allele_part, gene_part, pipeline_names)

            try:
                add2sample(this_allele_name, base_allele_name, sample.id, sample.patient.id, kdiff, pipeline_name, allele_snps, freq_by_clone, freq_by_seq, count,
                    total_count, session)
            except DbCreationError as e:
                print(e)
                
    processed_gene_types.extend(my_processed_types)
    return processed_gene_types


# Convert an ASC genotype style name (e.g. IGHV3-23D*01/IGHV3-23*01) to a VDJbase-style name (IGHV3-23*01_3.23D_01)

def asc_name_to_vdjbase(allele_calls, session):
    """
    Convert an ASC genotype style name to a VDJbase-style name, e.g. IGHV3-23D*01/IGHV3-23*01 to IGHV3-23*01_3.23D_01.

    :param allele_calls: list of ASC-style genotype names to convert
    :type allele_calls: list[str]
    :param session: SQLAlchemy session object for database queries
    :type session: sqlalchemy.orm.Session
    :return: the converted VDJbase-style name, or None if there is an error
    :rtype: str or None
    :raises DbCreationError: if the input allele_calls contain unexpected ambiguous allele format, extension-style snp, or different SNPs
    """
    parsed_calls = []

    for allele_call in allele_calls:
        base_allele, ambig_alleles, allele_snps, ext_snps = parse_allele(allele_call)
        parsed_calls.append((base_allele, ambig_alleles, allele_snps, ext_snps))

    # we don't expect any ext_snps in these names

    ext_snps = [p[3] for p in parsed_calls if len(p[3]) > 0]
    if len(ext_snps) > 0:
        raise DbCreationError(f"Error processing allele call {allele_calls}: unexpected extension-style snp found")
        return None

    # we don't expect any ambiguous alleles in these names

    ambig_alleles = [p[1] for p in parsed_calls if len(p[1]) > 0]
    if len(ambig_alleles) > 0:
        raise DbCreationError(f"Error processing allele call {allele_calls}: unexpected ambiguous allele format found")
        return None

    # we expect SNPs of all alleles to be the same

    allele_snps = list(set(['_'.join(p[2]) for p in parsed_calls if len(p[2]) > 0]))
    if len(allele_snps) > 1:
        raise DbCreationError(f"Error processing allele call {allele_calls}: alleles have different SNPs")
        return None

    allele_snps = allele_snps[0] if allele_snps else ''

    # remove any alleles that are marked as 'similar' and hence don't have their own row in the database

    alleles = [find_allele_or_similar(x[0], session) for x in parsed_calls]
    alleles = list(set(alleles))
    allele_names = [allele.name for allele in alleles]

    def gene_number(name):
        if '-' in name:
            num = ''.join(name.split('-')[1:])
        else:
            num = name[4:]         # e.g. for J genes

        if '*' in num:
            num = num.split('*')[0]

        return num

    allele_names.sort(key=lambda n: int(gene_number(n)) if gene_number(n).isnumeric() else 1000)
    vdjbase_name = allele_names[0]

    for allele_name in allele_names[1:]:
        vdjbase_name += f"_{gene_number(allele_name)}.{allele_name.split('*')[1]}"

    if allele_snps:
        vdjbase_name += '_' + allele_snps

    return vdjbase_name


def process_tigger_genotype(sample, processed_gene_types, pipeline_names, session):
    """
    Process the TIGGER genotype data for a given sample.
    ignore any rows whose gene type is listed in processed_gene_types.
    this lets us, for example, take Vs from asc, and Ds and Js from tigger

    :param sample: the sample object
    :type sample: obj
    :param processed_gene_types: a list of processed gene types
    :type processed_gene_types: list
    :param pipeline_names: a list of pipeline names
    :type pipeline_names: list
    :param session: the database session
    :type session: obj
    """
    print(sample.genotype)
    genotype = simple.read_csv(sample.genotype, delimiter='\t')
    for row in genotype:
        gene = row["gene"]
        gene_type = gene[3]

        if gene_type in processed_gene_types:
            continue

        if row["GENOTYPED_ALLELES"] == "NA" or row["GENOTYPED_ALLELES"] == '':
            continue

        allele_scores = {}
        if 'k_diff' in row:
            for allele in row["GENOTYPED_ALLELES"].split(","):
                allele_scores[allele] = row["k_diff"]
        elif 'z_score' in row:
            scores = str(row["z_score"]).split(',')
            for index, allele in enumerate(row["GENOTYPED_ALLELES"].split(",")):
                if index < len(scores):
                    allele_scores[allele] = scores[index]
                else:
                    allele_scores[allele] = '0'
                    print(f"Error: z_score for gene {gene} does not have enough values")
            # print(f"{','.join(allele_scores.keys())} {','.join(allele_scores.values())} ")

        # allele counts come from un-genotyped column
        allele_counts = {}
        ac = row["counts"].split(',')
        total_count = sum([int(x) if x.isnumeric() else 0 for x in ac])
        for index, allele in enumerate(str(row["alleles"]).split(",")):
            if not ((allele == "Unk") or ("NR" in allele) or ("del" in allele.lower())) and index < len(ac):
                allele_counts[allele] = ac[index]
        for index, allele in enumerate(str(row["GENOTYPED_ALLELES"]).split(",")):
            if (allele == "Unk") or ("NR" in allele):
                continue
            elif len(str(allele)) == 1:
                allele = "0" + str(allele)
            elif ("del" in allele.lower()):
                allele = "Del"

            # check if the allele exists in the genotype according to the clone size
            freq_by_clone = 0
            freq_by_seq = 0
            count = 0
            score = 0
            if allele != "Del":
                # check for bug found in TRB genotypes
                if len(str(row["Freq_by_Clone"]).split(";")) <= index:
                    print('Error: FREQ_BY_CLONE for gene %s does not have enough values' % gene)
                    continue

                # if freq_by_clone and freq_by_seq are zero, the allele has made it into the genotype but
                # there is no unambiguous support for it, so ogrdbstats won't report it. Let's live with
                # the inconsistency for the time being: hopefully will be cleared up by the ref book

                freq_by_clone = int(float(str(row["Freq_by_Clone"]).split(";")[index]))

                if len(str(row["Freq_by_Seq"]).split(";")) <= index:
                    print('Error: FREQ_BY_SEQ for gene %s does not have enough values' % gene)
                    continue

                freq_by_seq = int(float(str(row["Freq_by_Seq"]).split(";")[index]))
                count = int(allele_counts[allele]) if allele in allele_counts else 0
                score = float(allele_scores[allele]) if allele in allele_scores else 0

            allele_snps, base_allele_name, pipeline_name, this_allele_name = parse_ambiguous_allele(allele, gene, pipeline_names)

            try:
                add2sample(this_allele_name, base_allele_name, sample.id, sample.patient.id, score, pipeline_name, allele_snps, freq_by_clone, freq_by_seq, count,
                           total_count, session)
            except DbCreationError as e:
                print(e)
            

def find_allele_or_similar(allele_name, session):
    """
    Each unique sequence is only present on one single row in Allele. If multiple allele names
    correspond to the same sequence, they are listed in that row in the 'similar' field.
    This function returns that row, given the allele name

    :param allele_name:
    :type allele_name:
    :param session:
    :type session:
    :return:
    :rtype:
    """
    try:
        allele = session.query(Allele).filter(Allele.name == allele_name).one_or_none()
    except Exception as e:
        print('Multiple rows found for allele %s: check fasta files in reference directory.' % allele_name)

    if allele is not None:
        return allele

    allele = session.query(Allele).filter(Allele.similar.ilike(("%|" + allele_name + "|%"))).one_or_none()
    return allele


def add2sample(allele_name, base_allele_name, sample_id, pid, kdiff, pipeline_name, allele_snps, freq_by_clone, freq_by_seq, count, total_count, session):
    """
    Add a row to AllelesSample reflecting the presence of this allele in the sample.
    If the allele is not present in Allele already, add it there
    if the allele has a 'pipeline name', translate it to the 01_02 form

    :param allele_name: The name of the allele to add.
    :type allele_name: str
    :param base_allele_name: The base name of the allele.
    :type base_allele_name: str
    :param sample_id: The ID of the sample to which the allele belongs.
    :type sample_id: int
    :param pid: The patient ID.
    :type pid: int
    :param kdiff: The kdiff value.
    :type kdiff: float
    :param pipeline_name: The name of the pipeline.
    :type pipeline_name: str
    :param allele_snps: The SNPs of the allele.
    :type allele_snps: str
    :param freq_by_clone: The allele frequency by clone.
    :type freq_by_clone: float
    :param freq_by_seq: The allele frequency by sequence.
    :type freq_by_seq: float
    :param count: The count of the allele.
    :type count: int
    :param total_count: The total count of the allele.
    :type total_count: int
    :param session: The database session.
    :type session: SQLAlchemy session object
    :return: None
    """
    kdiff = float(kdiff)
    if math.isnan(kdiff):
        kdiff = 0.0

    allele = find_allele_or_similar(allele_name, session)
    if allele is None:
        allele = new_allele(allele_name, base_allele_name, pipeline_name, allele_snps, session)

    if allele_name not in new_alleles:      # i.e. it was in the reference set
        new_alleles[allele_name] = {
            'allele_name': allele_name,
            'pipeline_names': [pipeline_name],
            'similar': []
        }

    if len(pipeline_name) > 0:
        if allele.pipeline_name == '':
            allele.pipeline_name = pipeline_name
        else:
            pns = allele.pipeline_name.split(', ')
            if pipeline_name not in pns:
                pns.append(pipeline_name)
                allele.pipeline_name = ', '.join(pns)

        if pipeline_name not in new_alleles[allele_name]['pipeline_names']:
            new_alleles[allele_name]['pipeline_names'].append(pipeline_name)

    alleles_sample = AllelesSample(
        hap='geno',
        kdiff=kdiff,
        freq_by_clone=freq_by_clone,
        freq_by_seq=freq_by_seq,
        count=count,
        total_count=total_count,
        sample_id=sample_id,
        patient_id=pid,
        allele_id=allele.id
    )

    asc = session.query(AllelesSample)\
        .filter(AllelesSample.allele_id == alleles_sample.allele_id)\
        .filter(AllelesSample.sample_id == alleles_sample.sample_id)\
        .filter(AllelesSample.hap == 'geno').count()

    if asc == 0:
        session.add(alleles_sample)



def new_allele(allele_name, base_allele_name, pipeline_name, allele_snps, session):
    """
    Create a new Allele object in the database and return it. If the sequence of the allele is the same as an existing
    allele, return that instead.

    :param allele_name: the name of the allele as it will appear in VDJbase
    :type allele_name: str
    :param base_allele_name: the name without snps. It can include ambiguous alleles.
    :type base_allele_name: str
    :param pipeline_name: the name of the allele as found in the genotype file
    :type pipeline_name: str
    :param allele_snps: list of SNPs for the new allele
    :type allele_snps: List[str]
    :param session: the SQLAlchemy session for database operations
    :type session: Any
    :return: the newly created Allele object
    :rtype: Allele
    :raises DbCreationError: if the base allele is not found in the reference set
    """
    xx = f"allele_name: {allele_name}, base_allele_name: {base_allele_name}"

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
        raise DbCreationError('Error processing allele %s: base allele %s not in reference set. Allele not added.' % (allele_name, first_allele_name))

    # Check for a compound gene in the pipeline name, set gene id accordingly
    base_gene = pipeline_name.split('*')[0]
    if base_gene in compound_genes:
        gene_id = session.query(Gene.id).filter(Gene.name == compound_genes[base_gene]).one_or_none()[0]
    else:
        gene_id = first_allele.gene_id

    seq = first_allele.seq
    seq = "".join(seq.split())

    final_allele_name, is_novel_allele, seq = create_merged_sequence(allele_name, ambiguous_alleles, base_allele_name, seq, session)
    same_seq_allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()    # we only expect one row in Allele for each sequence

    if same_seq_allele is not None:
        allele = None
        temp = same_seq_allele.name
        sim = same_seq_allele.similar

        if (temp == final_allele_name):
            allele = same_seq_allele
        else:
            # if we're adding a gene-ambiguous allele, it should become the allele name

            if '.' in final_allele_name.split('*')[1]:
                final_allele_name, temp = temp, final_allele_name

            # TODO: do we really need these bars?
            if sim is not None and len(sim) > 0:
                if (final_allele_name in sim):
                    allele = same_seq_allele
                else:
                    sim = sim + ", |" + final_allele_name + "|"
            else:
                sim = "|" + final_allele_name + "|"

        if not allele:
            allele_table = Allele.__table__
            print(f"adding similar allele {temp}")

            if temp not in new_alleles:         # has changed its name through being gene_ambiguous
                new_alleles[temp] = {
                    'allele_name': temp,
                    'pipeline_names': [pipeline_name],
                    'similar': [x for x in sim.split('|') if len(x) > 0]
                }
            else:
                new_alleles[temp]['similar'].append(final_allele_name)

            stmt = allele_table.update().where(allele_table.c.seq == seq).values(similar=sim, name=temp)
            session.execute(stmt)
            session.commit()
            allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()

    else:
        ambig = 'ambiguous' if len(ambiguous_alleles) > 1 else ''
        print(f"adding {ambig} allele {final_allele_name}")

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

        new_alleles[final_allele_name] = {
            'allele_name': final_allele_name,
            'pipeline_names': [pipeline_name],
            'similar': []
        }

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


def create_merged_sequence(allele_name, ambiguous_alleles, base_allele_name, seq, session):
    """
    Creates a merged 'ambiguous' sequence for a given allele name by checking that the sequence for each allele referenced
    in the name as being indistinguishable is in the reference set. This function also takes into account any SNPs in the
     name. If any incorrectly-formatted terms are encountered, a warning is issued and the term is ignored.

    Example: TRBV3-1*01_02_2.01_2.02_2.03_a234g. Here there are 4 alleles to check and
    merge, plus one mutation. Check for syntax errors in the allele designation.

    Warn (and ignore) if any incorrectly-formatted terms are encountered.

    :param allele_name: the name of the allele to be merged
    :type allele_name: str
    :param ambiguous_alleles: a list of ambiguous alleles
    :type ambiguous_alleles: list
    :param base_allele_name: the base name of the allele
    :type base_allele_name: str
    :param seq: the sequence of the allele
    :type seq: str
    :param session: a session object for connecting to the database
    :type session: sqlalchemy.orm.session.Session
    :return: a tuple of the final allele name, a boolean indicating whether the allele is novel, and the merged sequence
    :rtype: tuple
    :raises DbCreationError: if the base allele is not in the reference set
    """
    final_allele_name = base_allele_name

    if allele_name != base_allele_name:
        rep = allele_name.replace(final_allele_name + '_', '').lower().split('_')
        is_novel_allele = False

        for r in rep:
            if re.search(allele_pattern, r) or re.search(gene_pattern, r):
                allele_seq_1 = seq

                if r not in ambiguous_alleles:
                    ambiguous_alleles.append(r)

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
                    allele_seq_1 += "n" * length_diff
                elif length_diff < 0:
                    allele_seq_2 += "n" * (length_diff * -1)

                for nuc1, nuc2 in zip(allele_seq_1, allele_seq_2):
                    if nuc1 == nuc2:
                        seq += nuc1
                    else:
                        seq += "n"

                final_allele_name += '_' + r

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

    return final_allele_name, is_novel_allele, seq


# Construct and add a 'compound gene', e.g. TRBV5/6, based on an allele name, eg TRBV6-5*01_6.01_6.02_6.03
def add_compound_gene(session, vdjbase_allele_name, pipeline_gene_name):
    """
    Add a compound gene to the database based on an allele name and pipeline gene name.

    :param session: The database session to use.
    :type session: sqlalchemy.orm.session.Session
    :param vdjbase_allele_name: The name of the VDJbase allele to use as the basis for the compound gene.
    :type vdjbase_allele_name: str
    :param pipeline_gene_name: The name of the pipeline gene to add the compound gene to.
    :type pipeline_gene_name: str
    :return: None
    :rtype: None

    Construct and add a 'compound gene', e.g. TRBV5/6, based on an allele name, eg TRBV6-5*01_6.01_6.02_6.03

    If the pipeline gene name already has a compound gene associated with it in the database, this function does nothing.

    After constructing the gene name, this function adds a new Gene object to the database with the appropriate attributes,
    including the newly constructed gene name. The locus_order and alpha_order attributes are set to the maximum values already
    present in the database plus one. The pseudo_gene attribute is set to False.
    """
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
    """
    Add *Del Allele records for each gene

    :param session: A SQLAlchemy session object.
    :return: A list of strings, with one item indicating the start of the function operation.
    :rtype: list of str
    """
    result = ['Adding deleted alleles']
    genes = session.query(Gene).all()

    for gene in genes:
        a = Allele(
            name=gene.name + '*Del',
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
    """
    Process haplotype files and update genotype information for each sample in the given dataset.

    :param ds_dir: A string representing the path to the dataset directory.
    :param species: A string representing the species of the samples.
    :param dataset: A string representing the dataset name.
    :param session: A SQLAlchemy session object.

    :return: A list of strings containing the processing result message.

    :raises OSError: If the dataset directory does not exist.

    :raises Exception: If any error occurs while processing the haplotype files.
    """
    result = ['Processing haplotype files']
    samples = session.query(Sample).all()

    for sample in samples:
        sample_dir = os.path.join('samples', sample.study.study_name, sample.patient.patient_name) #old format

        if not os.path.isdir(os.path.join(ds_dir, sample_dir)):
            sample_dir = os.path.join('samples', sample.study.study_name, sample.sample_name) #new format

        if os.path.isdir(os.path.join(ds_dir, sample_dir)):
            for filename in os.listdir(os.path.join(ds_dir, sample_dir)):
                if sample.sample_name in filename:
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
                print("No genotype report for sample %s" % sample.sample_name)
            if sample.genotype_stats is None:
                print("No genotype stats for sample %s" % sample.sample_name)
        else:
            print("No sample directory for sample %s" % sample.sample_name)

    session.flush()
    return result


def process_haplotype(filename, sample, haplo_gene, session):
    gene = "-".join(haplo_gene.split("-")[:-1])
    alleles = haplo_gene.split("-")[-1]
    allele1 = gene + "_" + alleles.split("_")[0]
    allele2 = gene + "_" + alleles.split("_")[1]

    hf =   HaplotypesFile(
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

