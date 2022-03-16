# Download data for rep-seq samples
from Bio import SeqIO
from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report
from api.reports.report_utils import make_output_file, chunk_list
from app import app, vdjbase_dbs
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype, Allele, AllelesSample, Gene, AlleleConfidenceReport
from db.vdjbase_airr_model import GenoDetection, SeqProtocol, Study, TissuePro, Patient, Sample
import csv
import zipfile
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params, find_vdjbase_sequences, \
    sequence_filters, find_vdjbase_samples
from sqlalchemy import func
import pandas as pd
from api.vdjbase.vdjbase import sample_info_filters
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SAMPLE_CHUNKS = 500

def zipdir(path, ziph, arc_root):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            path = os.path.join(root, file)
            ziph.write(path, arcname=path.replace(arc_root, ''))


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if 'Sample info' in params['type']:
        samples_by_dataset = {}
        for rep_sample in rep_samples:
            if rep_sample['dataset'] not in samples_by_dataset:
                samples_by_dataset[rep_sample['dataset']] = []
            samples_by_dataset[rep_sample['dataset']].append(rep_sample['sample_name'])

        attribute_query = []
        headers = []

        for name, filter in sample_info_filters.items():

            if filter['model'] is not None:
                attribute_query.append(filter['field'])
                headers.append(name)

        rows = []
        for dataset in samples_by_dataset.keys():
            session = vdjbase_dbs[species][dataset].session

            for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
                sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
                sample_list = [s[0] for s in sample_list]

                results = session.query(*attribute_query)\
                    .join(GenoDetection, GenoDetection.id == Sample.geno_detection_id)\
                    .join(Patient, Patient.id == Sample.patient_id)\
                    .join(SeqProtocol)\
                    .join(TissuePro)\
                    .join(Study, Sample.study_id == Study.id)\
                    .filter(Sample.sample_name.in_(sample_list)).all()

                rows.extend(results)

        outfile = make_output_file('csv')
        with open(outfile, 'w', newline='') as fo:
            writer = csv.writer(fo, dialect='excel')
            writer.writerow(headers)
            for row in rows:
                writer.writerow(row)

        return send_report(outfile, 'csv', attachment_filename='sample_info.csv')

    elif 'Sample files' in params['type']:
        outfile = make_output_file('zip')
        with zipfile.ZipFile(outfile, 'w', zipfile.ZIP_DEFLATED) as fo:
            samples_by_dataset = {}
            for rep_sample in rep_samples:
                if rep_sample['dataset'] not in samples_by_dataset:
                    samples_by_dataset[rep_sample['dataset']] = []
                samples_by_dataset[rep_sample['dataset']].append(rep_sample['sample_name'])

            added_files = []            # handle multiple samples in same dir etc
            added_dirs = []
            for dataset in samples_by_dataset.keys():
                print('adding dataset')
                session = vdjbase_dbs[species][dataset].session
                for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
                    sample_list = session.query(Sample.genotype, Sample.igsnper_plot_path).filter(Sample.sample_name.in_(sample_chunk)).all()
                    for p1, p2 in sample_list:
                        if p1 is not None and len(p1) > 0:
                            sample_dir = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, os.path.dirname(p1.replace('samples/', '')))
                            if sample_dir not in added_dirs:
                                zipdir(sample_dir, fo, os.path.join(VDJBASE_SAMPLE_PATH, species))        # sample files
                                added_dirs.append(sample_dir)
                        if p2 is not None and len(p2) > 0:
                            igsnper_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, p2)
                            if igsnper_path not in added_files:
                                fo.write(igsnper_path, arcname=igsnper_path.replace(os.path.join(VDJBASE_SAMPLE_PATH, species), ''))
                                added_files.append(igsnper_path)

        return send_report(outfile, 'zip', attachment_filename='sample_data.zip')

    elif 'Ungapped' in params['type'] or 'Gapped' in params['type']:
        required_cols = ['name', 'seq', 'dataset']
        seqs = find_sequences(params, rep_samples, species, required_cols)

        recs = []
        for seq in seqs:
            id = '%s|%s|%s' % (seq['name'], species, seq['dataset'])
            recs.append(SeqRecord(Seq(seq['seq'] if 'Gapped' in params['type'] else seq['seq'].replace('.', '')), id=id, description=''))

        outfile = make_output_file('fasta')
        SeqIO.write(recs, outfile, "fasta")
        return send_report(outfile, 'fasta', attachment_filename='%s_sequences.fasta' % species)

    elif 'Gene info' in params['type']:
        headers = []
        for name, att_filter in sequence_filters.items():
            if att_filter['model'] is not None:
                headers.append(name)

        headers.append('dataset')
        rows = find_sequences(params, rep_samples, species, headers)

        outfile = make_output_file('csv')
        with open(outfile, 'w', newline='') as fo:
            writer = csv.DictWriter(fo, dialect='excel', fieldnames=headers)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        return send_report(outfile, 'csv', attachment_filename='sequence_info.csv')

    raise BadRequest('No output from report')


def find_sequences(params, rep_samples, species, required_cols):
    # because the run api is samples-oriented, we have to do a little work to recover the datasets. We don't need
    # the calculated samples. Not ideal, could consider changing the api if we keep bumping in to this
    datasets = []
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in datasets:
            datasets.append(rep_sample['dataset'])
    seqs = find_vdjbase_sequences(species, datasets, required_cols, params['filters'])
    return seqs

