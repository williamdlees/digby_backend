from app import app, vdjbase_dbs
import os.path
import os

import yaml

from db.vdjbase_airr_model import SeqProtocol, Study, TissuePro, Patient, Sample
from db.vdjbase_model import *


def export_metadata():
    export_dir = os.path.join(app.config['EXPORT_DIR'], 'vdjbase_metadata')

    if not os.path.isdir(export_dir):
        os.mkdir(export_dir)

    for species in vdjbase_dbs.keys():
        if not os.path.isdir(os.path.join(export_dir, species)):
            os.mkdir(os.path.join(export_dir, species))
        for dataset in vdjbase_dbs[species]:
            if not os.path.isdir(os.path.join(export_dir, species, dataset)):
                os.mkdir(os.path.join(export_dir, species, dataset))
            if os.path.isfile(os.path.join(export_dir, species, dataset, 'projects.yml')):
                os.remove(os.path.join(export_dir, species, dataset, 'projects.yml'))
            session = vdjbase_dbs[species][dataset].get_session()
            studies = session.query(Study).all()

            meta = {}
            for study in studies:
                samples = session.query(Sample)\
                    .join(SeqProtocol) \
                    .join(GenoDetection) \
                    .join(TissuePro) \
                    .filter(Sample.study == study).all()
                patients = session.query(Patient).filter(Patient.study == study).all()
                rec = {
                    'Project': study.study_title,
                    'Researcher': study.submitted_by,
                    'Institute': study.lab_address,
                    'Number of Subjects': len(patients),
                    'Number of Samples': len(samples),
                    'Reference': study.pub_ids,
                    'Contact': study.study_contact,
                    'Accession id': study.study_id,
                    'Accession reference': study.accession_reference,
                    'Subjects': {},
                    'Sequence Protocol': {},
                    'Tissue Processing': {},
                    'Genotype Detections': {},
                    'Samples': {},
                }

                for patient in patients:
                    rec['Subjects'][patient.patient_name] = {
                        'Name': patient.patient_name,
                        'Original name': patient.name_in_paper,
                        'Sex': patient.sex,
                        'Ethnic': patient.ethnicity,
                        'Country': patient.ancestry_population,
                        'Health Status': patient.disease_diagnosis_label,
                        'Age': patient.age,
                        'Cohort': patient.study_group_description
                    }

                for sample in samples:
                    rec['Samples'][sample.sample_name] = {
                        'Name': sample.sample_name,
                        'Chain': sample.chain,
                        'Date': sample.date,
                        'Reads': sample.row_reads,
                        'Sample Group': sample.samples_group,
                        'Subject Name': sample.patient.patient_name,
                        'Sequence Protocol Name': sample.seq_protocol.name,
                        'Tissue Processing Name': sample.tissue_pro.name,
                        'Genotype Detection Name': sample.geno_detection.name,
                    }

                    if sample.seq_protocol.name not in rec['Sequence Protocol']:
                        rec['Sequence Protocol'][sample.seq_protocol.name] = {
                            'Name': sample.seq_protocol.name,
                            'Sequencing_platform': sample.seq_protocol.sequencing_platform,
                            'Sequencing_length': sample.seq_protocol.read_length,
                            'UMI': sample.seq_protocol.umi != 0,
                            'Helix': sample.seq_protocol.helix,
                            'Primer 5 location': sample.seq_protocol.reverse_pcr_primer_target_location,
                            'Primer 3 location': sample.seq_protocol.forward_pcr_primer_target_location,
                        }
                    if sample.tissue_pro.name not in rec['Tissue Processing']:
                        rec['Tissue Processing'][sample.tissue_pro.name] = {
                            'Name': sample.tissue_pro.name,
                            'Species': sample.tissue_pro.cell_species_label,
                            'Tissue': sample.tissue_pro.tissue_label,
                            'Cell Type': sample.tissue_pro.cell_subset_label,
                            'Sub Cell Type': sample.tissue_pro.sub_cell_type,
                            'Isotype': sample.tissue_pro.cell_phenotype,
                        }
                    if sample.geno_detection.name not in rec['Genotype Detections']:
                        rec['Genotype Detections'][sample.geno_detection.name] = {
                            'Name': sample.geno_detection.name,
                            'Repertoire or Germline': sample.geno_detection.detection,
                            'Pre-processing': sample.geno_detection.prepro_tool,
                            'Aligner Tool': sample.geno_detection.aligner_tool,
                            'Aligner Version': sample.geno_detection.aligner_ver,
                            'Germline Reference': sample.geno_detection.aligner_reference,
                            'Genotyper Tool': sample.geno_detection.geno_tool,
                            'Genotyper Version': sample.geno_detection.geno_ver,
                            'Haplotyper Tool': sample.geno_detection.haplotype_tool,
                            'Haplotyper Version': sample.geno_detection.haplotype_ver,
                            'Single Assignment': sample.geno_detection.single_assignment,
                        }

                meta[study.study_title] = rec

            with open(os.path.join(export_dir, species, dataset, 'projects.yml'), 'w') as fo:
                yaml.dump(meta, stream=fo, allow_unicode=True)

    return 'Export complete!'
