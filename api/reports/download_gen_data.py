# Download data for genomic samples
import csv
import zipfile
import os
import glob
from werkzeug.exceptions import BadRequest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from api.genomic.genomic import find_genomic_samples, find_genomic_sequences
from db.genomic_api_query_filters import genomic_sequence_filters, genomic_sample_filters
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file
from app import app
from db.genomic_airr_model import Sample


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
        headers = genomic_sample_filters.keys()

        attribute_query = []
        for name, filter in genomic_sample_filters.items():
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
            # add reference sequence
            for genomic_dataset in genomic_datasets:
                for fn in [species + '*.gff3', species + '*.fasta']:
                    files = glob.glob(os.path.join(app.config['STATIC_PATH'], 'study_data', 'Genomic', 'samples', species, genomic_dataset, fn))
                    for file in files:
                        fo.write(file, arcname=file.replace(app.config['STATIC_PATH'], ''))

            # add samples
            added_dirs = []
            sample_paths = find_genomic_samples([Sample.annotation_path], species, genomic_datasets, params['filters'])
            sample_paths = [s for s in sample_paths if '.csv' in s['annotation_path']] # remove any null paths
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
        seqs, _, _ = find_genomic_sequences(required_cols, genomic_datasets, species, params['filters'])

        recs = []
        for seq in seqs:
            if len(seq[seq_name]) > 0:
                id = '%s|%s|%s' % (seq['name'], species, seq['dataset'])
                recs.append(SeqRecord(Seq(seq[seq_name]), id=id, description=''))

        outfile = make_output_file('fasta')
        SeqIO.write(recs, outfile, "fasta")
        return send_report(outfile, 'fasta', attachment_filename='%s_sequences.fasta' % species)

    elif 'Gene info' in params['type']:
        headers = list(genomic_sequence_filters.keys())
        rows, headers, alias_dict = find_genomic_sequences(headers, genomic_datasets, species, params['filters'])

        for h in headers:
            if h in alias_dict:
                headers[headers.index(h)] = alias_dict[h]

        outfile = make_output_file('csv')
        with open(outfile, 'w', newline='') as fo:
            writer = csv.DictWriter(fo, dialect='excel', fieldnames=headers)
            writer.writeheader()
            for row in rows:
                for r in list(row.keys()):
                    if r in alias_dict:
                        row[alias_dict[r]] = row[r]
                        del row[r]
                writer.writerow(row)

        return send_report(outfile, 'csv', attachment_filename='sequence_info.csv')

    raise BadRequest('No output from report')


