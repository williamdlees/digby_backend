# Code to import projects, driven by the contents of projects.yml
import copy
import os
import os.path
from collections import defaultdict
from glob import glob
import re
from db.vdjbase_exceptions import *
import yaml
from db.vdjbase_model import Study, GenoDetection, SeqProtocol, TissuePro, Patient, Sample, Allele, AllelesSample, Gene, HaplotypesFile, SamplesHaplotype
import pandas as pd
import math
from db.vdjbase_formats import *
import csv


from sqlalchemy import update

def import_studies(ds_dir, species, dataset, session):
    result = []
    study_data = {}

    study_file = os.path.join(ds_dir, 'projects.yml')
    if os.path.isfile(study_file):
        with open(study_file, 'r') as fi:
            study_data = yaml.safe_load(fi)
    else:
        # try to consolidate data from sample-level files
        study_data = consolidate_metadata(ds_dir)

    for sd in study_data.values():
        try:
            s = Study(
                name=sd['Project'],
                institute=sd['Institute'],
                researcher=sd['Researcher'],
                num_subjects=sd['Number of Subjects'],
                num_samples=sd['Number of Samples'],
                reference=sd['Reference'],
                contact=sd['Contact'],
                accession_id=sd['Accession id'],
                accession_reference=sd['Accession reference'],
            )
            session.add(s)

            study_rec = session.query(Study).filter(Study.name == sd['Project']).one_or_none()

            for gd in sd['Genotype Detections'].values():
                g = session.query(GenoDetection).filter(GenoDetection.name == gd['Name']).one_or_none()

                if g is None:
                    g = GenoDetection(
                        name=gd['Name'],
                        prepro_tool=gd['Pre-processing'] if gd['Pre-processing'] else '',
                        aligner_tool=gd['Aligner Tool'],
                        aligner_ver=gd['Aligner Version'],
                        aligner_reference=gd['Germline Reference'],
                        geno_tool=gd['Genotyper Tool'],
                        geno_ver=gd['Genotyper Version'],
                        haplotype_tool=gd['Haplotyper Tool'],
                        haplotype_ver=gd['Haplotyper Version'],
                        single_assignment=gd['Single Assignment'],
                        detection=gd['Repertoire or Germline'],
                    )
                    session.add(g)

            for sp in sd['Sequence Protocol'].values():
                p = session.query(SeqProtocol).filter(SeqProtocol.name == sp['Name']).one_or_none()

                if p is None:
                    p = SeqProtocol(
                        name = sp['Name'],
                        umi=sp['UMI'],
                        sequencing_length=sp['Sequencing_length'],
                        primers_3_location=sp['Primer 3 location'],
                        primers_5_location=sp['Primer 5 location'],
                        sequencing_platform=sp['Sequencing_platform'],
                        helix=sp['Helix'],
                    )
                    session.add(p)

            for tp in sd['Tissue Processing'].values():
                t = session.query(TissuePro).filter(TissuePro.name == tp['Name']).one_or_none()

                if t is None:
                    t = TissuePro(
                        name=tp['Name'],
                        species=tp['Species'],
                        tissue=tp['Tissue'],
                        cell_type=tp['Cell Type'],
                        sub_cell_type=tp['Sub Cell Type'],
                        isotype=tp['Isotype'],
                    )
                    session.add(t)

            for pa in sd['Subjects'].values():
                p = Patient(
                    name=pa['Name'],
                    sex=pa['Sex'],
                    ethnic=pa['Ethnic'],
                    status=pa['Health Status'],
                    cohort=pa['Cohort'],
                    age=pa['Age'],
                    study_id=study_rec.id,
                    name_in_paper=pa['Original name'],
                    country=pa['Country'],
                    igsnper_sample_id=0,
                )
                session.add(p)

            for sa in sd['Samples'].values():
                gd_id = session.query(GenoDetection.id).filter(GenoDetection.name == sa['Genotype Detection Name']).one_or_none()
                if gd_id is None:
                    result.append(['sample %s: Genotype Detection record not found' % sa['Name']])
                sp_id = session.query(SeqProtocol.id).filter(SeqProtocol.name == sa['Sequence Protocol Name']).one_or_none()
                if sp_id is None:
                    result.append(['sample %s: Sequence Protocol record not found' % sa['Name']])
                tp_id = session.query(TissuePro.id).filter(TissuePro.name == sa['Tissue Processing Name']).one_or_none()
                if tp_id is None:
                    result.append(['sample %s: Tissue Processing record not found' % sa['Name']])
                pa_id = session.query(Patient.id).filter(Patient.name == sa['Subject Name']).one_or_none()
                if pa_id is None:
                    result.append(['sample %s: Subject record not found' % sa['Name']])
                s = Sample(
                    name=sa['Name'],
                    chain=sa['Chain'],
                    row_reads=sa['Reads'],
                    genotype='',
                    genotype_graph='',
                    date=sa['Date'],
                    samples_group=sa['Sample Group'],
                    geno_detection_id=gd_id[0],
                    patient_id=pa_id[0],
                    seq_protocol_id=sp_id[0],
                    study_id=study_rec.id,
                    tissue_pro_id=tp_id[0],
                    genotype_stats='',
                )
                session.add(s)
        except Exception as e:
            print('Exception processing metadata file')

    session.commit()
    result.append('Import completed!')
    return result


# raise DbCreationError('project file %s does not exist' % study_file)

# enumerate dirs and paths under the specified directory
def listdp(dir):
    dirs = [os.path.split(name)[0] for name in glob(os.path.join(dir, '*/'))]
    return zip([os.path.split(name)[1] for name in dirs], dirs)

# Consolidate metadata from yml files in each sample directory
def consolidate_metadata(export_dir):
    metadata = {}
    results = []
    for project_name, project_path in listdp(os.path.join(export_dir, 'samples')):
        for sample_name, sample_path in listdp(project_path):
            metadata_file = os.path.join(sample_path, sample_name + '.yml')
            if not os.path.isfile(metadata_file):
                print('metadata file %s not found: sample will not be included!' % metadata_file)
                results.append('metadata file %s not found: sample will not be included!' % metadata_file)
                continue

            try:
                with open(metadata_file, 'r') as fi:
                    study_data = yaml.safe_load(fi)

                if project_name not in study_data:
                    raise Exception('Project name %s not in sample %s metadata. Sample will not be included.' % (sample_name, project_name))

                if not isinstance(study_data[project_name]['Accession id'], str):
                    raise Exception('Accession id in metadata for sample %s, project %s is not a string. Sample will not be included.' % (sample_name, project_name))

                # fix up some potentially missing fields

                if study_data[project_name]['Samples'][sample_name]['Sequence Protocol Name'] is None:
                    study_data[project_name]['Samples'][sample_name]['Sequence Protocol Name'] = list(study_data[project_name]['Sequence Protocol'].keys())[0]

                if study_data[project_name]['Samples'][sample_name]['Tissue Processing Name'] is None:
                    study_data[project_name]['Samples'][sample_name]['Tissue Processing Name'] = list(study_data[project_name]['Tissue Processing'].keys())[0]

                if study_data[project_name]['Samples'][sample_name]['Reads'] is None:
                    study_data[project_name]['Samples'][sample_name]['Reads'] = 0

                if study_data[project_name]['Samples'][sample_name]['Chain'] is None:
                    study_data[project_name]['Samples'][sample_name]['Chain'] = 'IGH'

                # fix for fields changed with 16 oct release

                for gk, gd in study_data[project_name]['Genotype Detections'].items():
                    if 'Aligner Tool' not in gd and 'Alignment Tool' in gd:
                            study_data[project_name]['Genotype Detections'][gk]['Aligner Tool'] = gd['Alignment Tool']
                            if gd['Alignment reference v'] == gd['Alignment reference d'] and gd['Alignment reference v'] == gd['Alignment reference j']:
                                study_data[project_name]['Genotype Detections'][gk]['Germline Reference'] = gd['Alignment reference v']
                            else:
                                study_data[project_name]['Genotype Detections'][gk]['Germline Reference'] = '(%s (v), %s (d), %s (j))' %\
                                                     (gd['Alignment reference v'], gd['Alignment reference d'], gd['Alignment reference j'])
                    if isinstance(gd['Single Assignment'], str):
                        study_data[project_name]['Genotype Detections'][gk]['Single Assignment'] = ('t' in gd['Single Assignment'] or 'T' in gd['Single Assignment'])

                for sk, sp in study_data[project_name]['Sequence Protocol'].items():
                    if 'Primer.3.location' in sp and 'Primer 3 location' not in sp:
                        study_data[project_name]['Sequence Protocol'][sk]['Primer 3 location'] = sp['Primer.3.location']
                    if 'Primer.5.location' in sp and 'Primer 5 location' not in sp:
                        study_data[project_name]['Sequence Protocol'][sk]['Primer 5 location'] = sp['Primer.5.location']

                for sk, sp in study_data[project_name]['Subjects'].items():
                    if 'Health.Status' in sp and 'Health Status' not in sp:
                        study_data[project_name]['Subjects'][sk]['Health Status'] = sp['Health.Status']
                    if 'Original.name' in sp and 'Original name' not in sp:
                        study_data[project_name]['Subjects'][sk]['Original name'] = sp['Original.name']

                for tk, tp in study_data[project_name]['Tissue Processing'].items():
                    if 'Cell.Type' in tp and 'Cell Type' not in tp:
                        study_data[project_name]['Tissue Processing'][tk]['Cell Type'] = tp['Cell.Type']
                    if 'Sub.Cell.Type' in tp and 'Sub Cell Type' not in tp:
                        study_data[project_name]['Tissue Processing'][tk]['Sub Cell Type'] = tp['Sub.Cell.Type']

                if project_name not in metadata:
                    metadata[project_name] = copy.deepcopy(study_data[project_name])
                else:
                    if sample_name not in study_data[project_name]['Samples']:
                        raise Exception('Wrong sample name not in in sample %s metadata. Sample will not be included.' % sample_name)

                    metadata[project_name]['Samples'][sample_name] = copy.deepcopy(study_data[project_name]['Samples'][sample_name])

                    for section in ('Sequence Protocol', 'Subjects', 'Tissue Processing', 'Genotype Detections'):
                        for item in study_data[project_name][section].keys():
                            if item not in metadata[project_name][section]:
                                metadata[project_name][section][item] = copy.deepcopy(study_data[project_name][section][item])


            except Exception as e:
                print('Exception processing metadata file %s: %s' % (metadata_file, e))

    with open(os.path.join(export_dir, 'consolidated.yml'), 'w') as fo:
        fo.write(yaml.dump(metadata, default_flow_style=False))

    return metadata

def process_genotypes(ds_dir, species, dataset, session):
    result = ['Processing genotype files']
    samples = session.query(Sample).all()

    # Read names assigned to ambiguous alleles by the pipeline

    pipeline_names = {}   # allele name given pipeline name
    allele_names = {}     # pipeline name given allele name
    if os.path.isfile(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv')):
        with open(os.path.join(ds_dir, 'reference/ambiguous_allele_names.csv'), 'r') as fi:
            reader = csv.DictReader(fi)
            for row in reader:
                if len(row['ALLELES']) == 1:
                    allele = '0' + row['ALLELES']
                elif ',' in row['ALLELES'] or ' ' in row['ALLELES']:
                    allele = row['ALLELES'].replace(' ', '').replace(',', '_')
                else:
                    allele = row['ALLELES']

                pipeline_name = row['GENE'] + '*' + row['PATTERN']
                allele_name = row['GENE'] + '*' + allele
                allele_names[allele_name] = pipeline_name
                pipeline_names[pipeline_name] = allele_name

    for sample in samples:
        genotype_file = os.path.join('samples', sample.study.name, sample.patient.name, sample.name + '_geno_H_binom.tab').replace('\\', '/')

        if not os.path.isfile(os.path.join(ds_dir, genotype_file)):
            genotype_file = os.path.join('samples', sample.study.name, sample.name, sample.name + '_genotype.tsv').replace('\\', '/')  # new directory layout

        if not os.path.isfile(os.path.join(ds_dir, genotype_file)):
            genotype_file = os.path.join('samples', sample.study.name, sample.name, sample.name + '.tsv').replace('\\', '/')  # another new directory layout

        if os.path.isfile(os.path.join(ds_dir, genotype_file)):
            sample.genotype = genotype_file
            sample_genotype(os.path.join(ds_dir, genotype_file), sample.id, sample.patient.id, pipeline_names, allele_names, session)
            session.commit()
        else:
            result.append('Error: no genotype file for sample %s' % sample.name)

    return result



def sample_genotype(inputfile, sample_id, patient_id, pipeline_names, allele_names, session):
    """
    Upload genotype to sample.
    :param inputfile: the path to the genotype
    :param sample_id: sample id
    """
    haplo = "geno"
    genotype = pd.read_csv(inputfile, sep="\t")
    allele_pattern = re.compile("^[0-9]{2}$")
    mut_pattern = re.compile("^[A,G,T,C,a,g,t,c][0-9]+[A,G,T,C,a,g,t,c]$")

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

            # check if the allele exist in the genotype according to the clone size
            if allele != "Del":
                clone_size = int(row[FREQ_BY_CLONE].split(INT_SEP)[index])
                if not clone_size:
                    continue

            # parse allele, handling old ambiguous style (01_02) and new bp style
            ambig_alleles = ''
            allele_snps = ''
            base_allele = allele

            if '_' in allele:
                ambig_alleles = []
                allele_snps = []

                base_allele = allele.split('_')[0]
                rep = allele.split("_")[1:]

                for r in rep:
                    if re.search(allele_pattern, r):
                        ambig_alleles.append(r)
                    elif re.search(mut_pattern, r):
                        allele_snps.append(r)

                ambig_alleles = '_'.join(ambig_alleles)
                allele_snps = '_'.join(allele_snps)

            base_allele_name = gene + "*" + base_allele
            pipeline_name = ''

            if base_allele_name in pipeline_names:
                pipeline_name = base_allele_name
                base_allele_name = pipeline_names[base_allele_name]
            elif 'bp' in base_allele_name:        # will only happen if one of the 'bp' alleles has not been put in ambiguous_allele_names.csv
                # fudge for time being
                pipeline_name = base_allele_name
                base_allele_name = base_allele_name.replace('bp', '')
                pipeline_names[pipeline_name] = base_allele_name
                print('%s is not listed in ambiguous_allele_names.csv: assuming it corresponds to %s' % (pipeline_name, base_allele_name))

            if len(ambig_alleles) > 0:
                base_allele_name += '_' + ambig_alleles
                if len(pipeline_name) > 0:
                    pipeline_name += '_' + ambig_alleles

            if len(allele_snps) > 0:
                base_allele_name += '_' + allele_snps
                if len(pipeline_name) > 0:
                    pipeline_name += '_' + allele_snps

            add2sample(base_allele_name, sample_id, haplo, patient_id, kdiff, pipeline_name, session)

# Each unique sequence is only present on one single row in Allele. If multiple allele names
# correspond to the same sequence, they are listed in that row in the 'similar' field.
# This function returns that row, given the allele name

def find_allele_or_similar(allele_name, session):
    allele = session.query(Allele).filter(Allele.name == allele_name).one_or_none()

    if allele is not None:
        return allele

    allele = session.query(Allele).filter(Allele.similar.ilike(("%|" + allele_name + "|%"))).one_or_none()
    return allele if allele is not None else None


# Add a row to AllelesSample reflecting the presence of this allele in the sample.
# If the allele is not present in Allele already, add it there
# if the allele has a 'pipeline name', translate it to the 01_02 form

def add2sample (allele_name, sample_id, haplo, pid, kdiff, pipeline_name, session):
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
        allele = new_allele(allele_name, session)

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

def new_allele(allele_name, session):
    base_allele = find_allele_or_similar(allele_name.split("_")[0], session)

    if base_allele is None:
        # check for duplicate (D) allele
        allele_D_name = allele_name.split("*")[0] + "D*" + allele_name.split("*")[1].split("_")[0]
        base_allele = find_allele_or_similar(allele_D_name, session)

    if base_allele is None:
        raise DbCreationError('Error processing allele %s: base allele not in reference set' % allele_name)

    allele_pattern = re.compile("^[0-9]{2}$")
    mut_pattern = re.compile("^[A,G,T,C,a,g,t,c][0-9]+[A,G,T,C,a,g,t,c]$")
    gid = base_allele.gene_id
    seq = base_allele.seq
    seq = "".join(seq.split())

    final_allele_name = allele_name.split("_")[0]
    rep = allele_name.lower().split("_")[1:]
    is_novel_allele = False
    ambiguous_alleles = [allele_name.split("*")[1].split("_")[0]]
    for r in rep:
        if re.search(allele_pattern, r):
            if r not in ambiguous_alleles:
                final_allele_name += "_" + r
                allele_seq_1 = seq
                allele2_name = allele_name.split("*")[0] + "*" + r # r = allele for example 01

                allele2 = find_allele_or_similar(allele2_name, session)

                if allele2 is None:
                    allele2_D_name = allele_name.split("*")[0] + "D*" + r # r = allele for example 01
                    allele2 = find_allele_or_similar(allele2_D_name, session)

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

    same_seq_allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()    # we only expect one - one row in Allele for each sequence

    if same_seq_allele is not None:
        allele = None
        temp = same_seq_allele.name
        sim = same_seq_allele.similar
        if (temp == final_allele_name):
            allele = same_seq_allele
        else:
            if sim is not None and len(sim) > 0:
                if (final_allele_name in sim):
                    allele = same_seq_allele
                else:
                    sim = sim + ", |" + final_allele_name + "|"
            else:
                sim = "|" + final_allele_name + "|"

        if not allele:
            allele_table = Allele.__table__
            stmt = allele_table.update().where(allele_table.c.seq == seq).values(similar=sim)
            session.execute(stmt)
            allele = session.query(Allele).filter(Allele.seq == seq).one_or_none()

    else:
        allele = Allele(
            name=final_allele_name,
            seq=seq,
            seq_len=str(len(seq)),
            gene_id=gid,
            is_single_allele=len(ambiguous_alleles)==1,
            appears=0,
            low_confidence=False,
            novel=is_novel_allele,
            max_kdiff=0.0,
            similar='',
            pipeline_name='',
        )
        session.add(allele)
        session.flush()
        session.refresh(allele)

    return allele


# Add Allele records for *Del
def add_deleted_alleles(session):
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
    result = ['Processing haplotype files']
    samples = session.query(Sample).all()

    for sample in samples:
        sample_dir = os.path.join('samples', sample.study.name, sample.patient.name) #old format

        if not os.path.isdir(os.path.join(ds_dir, sample_dir)):
            sample_dir = os.path.join('samples', sample.study.name, sample.name) #new format

        if os.path.isdir(os.path.join(ds_dir, sample_dir)):
            for filename in os.listdir(os.path.join(ds_dir, sample_dir)):
                if sample.name in filename:
                    if 'haplotype.' in filename:
                        haplo_gene = filename.replace('_haplotype.tab', '').replace('_haplotype.tsv', '').split('_gene-')[1]
                        process_haplotype(os.path.join(sample_dir, filename).replace('\\', '/'), sample, haplo_gene, session)
                    elif 'ogrdb_plots' in filename:
                        sample.genotype_report = os.path.join(sample_dir, filename).replace('\\', '/')
                    elif 'ogrdb_report' in filename:
                        sample.genotype_stats = os.path.join(sample_dir, filename).replace('\\', '/')

            if sample.genotype_report is None:
                print("No genotype report for sample %s" % sample.name)
            if sample.genotype_stats is None:
                print("No genotype stats for sample %s" % sample.name)
        else:
            print("No sample directory for sample %s" % sample.name)

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
        haplotypes_files_id=hf.id,
    )
    session.add(sh)


