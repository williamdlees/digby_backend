# Download data for genomic samples
from Bio import SeqIO
from werkzeug.exceptions import BadRequest

from api.genomic.genomic import genomic_sequence_filters, find_genomic_samples, genomic_subject_filters, \
    find_genomic_sequences
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file
from app import app
from db.genomic_db import Subject, Sample
import csv
import zipfile
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def zipdir(path, ziph, arc_root):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            path = os.path.join(root, file)
            ziph.write(path, arcname=path.replace(arc_root, ''))


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(genomic_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if 'Sample info' in params['type']:
        headers = genomic_subject_filters.keys()

        attribute_query = []
        for name, filter in genomic_subject_filters.items():
            if filter['model'] is not None:
                attribute_query.append(filter['field'])

        rows = find_genomic_samples(attribute_query, species, genomic_datasets, params['filters'])

        outfile = make_output_file('csv')
        with open(outfile, 'w', newline='') as fo:
            writer = csv.DictWriter(fo, dialect='excel', fieldnames=headers)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        return send_report(outfile, 'csv', attachment_filename='sample_info.csv')

    elif 'Sample files' in params['type']:
        outfile = make_output_file('zip')
        with zipfile.ZipFile(outfile, 'w', zipfile.ZIP_DEFLATED) as fo:
            added_dirs = []
            sample_paths = find_genomic_samples([Sample.annotation_path], species, genomic_datasets, params['filters'])
            sample_paths = ['/'.join(['study_data/Genomic', s['annotation_path'].split('Genomic')[1]]) for s in sample_paths]
            sample_paths = [os.path.join(app.config['STATIC_PATH'], s) for s in sample_paths]
            for sample_path in sample_paths:
                sample_dir = os.path.dirname(sample_path)
                if sample_dir not in added_dirs:
                    zipdir(sample_dir, fo, app.config['STATIC_PATH'])        # sample files
                    added_dirs.append(sample_dir)

        return send_report(outfile, 'zip', attachment_filename='sample_data.zip')

    elif 'Ungapped' in params['type'] or 'Gapped' in params['type']:
        seq_name = 'sequence' if 'Ungapped' in params['type'] else 'gapped_sequence'
        required_cols = ['name', seq_name, 'dataset']
        seqs = find_genomic_sequences(required_cols, genomic_datasets, species, params['filters'])

        recs = []
        for seq in seqs:
            if len(seq[seq_name]) > 0:
                id = '%s|%s|%s' % (seq['name'], species, seq['dataset'])
                recs.append(SeqRecord(Seq(seq[seq_name]), id=id, description=''))

        outfile = make_output_file('fasta')
        SeqIO.write(recs, outfile, "fasta")
        return send_report(outfile, 'fasta', attachment_filename='%s_sequences.fasta' % species)

    elif 'Gene info' in params['type']:
        headers = genomic_sequence_filters.keys()
        rows = find_genomic_sequences(headers, genomic_datasets, species, params['filters'])

        outfile = make_output_file('csv')
        with open(outfile, 'w', newline='') as fo:
            writer = csv.DictWriter(fo, dialect='excel', fieldnames=headers)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

        return send_report(outfile, 'csv', attachment_filename='sequence_info.csv')

    raise BadRequest('No output from report')


