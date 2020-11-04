# Create or update IgSNPer records for a VDJbase dataset
import subprocess

from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH
from app import vdjbase_dbs, app
from db.vdjbase_model import Sample, Patient, Study, Gene
from api.reports.reports import make_output_file
import os.path
import shutil
from db.vdjbase_db import ContentProvider
from glob import glob


def do_igsnper(species, dataset):
    export_dir = os.path.join(app.config['EXPORT_DIR'], 'vdjbase_metadata')
    ds_dir = os.path.abspath(os.path.join(export_dir, species, dataset))
    if not os.path.isfile(os.path.join(ds_dir, 'db.sqlite3')):
        return 'No database found'

    igsnper_dir = os.path.join(ds_dir, 'igsnper')
    db = ContentProvider(os.path.join(ds_dir, 'db.sqlite3'))

    # remove any existing igsnper related database fields

    db.session.query(Gene).update({Gene.igsnper_plot_path: ''}, synchronize_session=False)
    db.session.query(Sample).update({Sample.igsnper_plot_path: ''}, synchronize_session=False)
    db.session.query(Patient).update({Patient.igsnper_sample_id: 0}, synchronize_session=False)
    db.session.commit()


    # Create table of tigger files

    tigger_file_name = make_output_file('txt')
    with open(tigger_file_name, 'w') as fo:
        header = "TiggerFilePath     ProjectID     SubjectID\n"
        fo.write(header)

        samples = db.session.query(Sample.genotype, Sample.name, Study.name, Patient.name).join(Study, Sample.study_id == Study.id).join(Patient, Sample.patient_id == Patient.id).all()

        for sample in samples:
            if 'S1' in sample[1]:
                tigger_file_path = os.path.join(VDJBASE_SAMPLE_PATH, species, dataset, sample[0].replace('samples/', ''))

                if os.path.isfile(tigger_file_path):
                    fo.write('%s     %s     %s\n' % (tigger_file_path, sample[2], sample[3]))

    cmd_line = ['python', os.path.join(app.config['IGSNPER_PATH'], 'ig_snper.py')]

    if os.path.isdir(igsnper_dir):
        shutil.rmtree(igsnper_dir, ignore_errors=True, onerror=None)

    args = ['-o', igsnper_dir, '-c', tigger_file_name]

    cmd_line.extend(args)
    print("Running IgSNPer: '%s'\n" % ' '.join(cmd_line))
    proc = subprocess.Popen(cmd_line, cwd=app.config['IGSNPER_PATH'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for stdout_line in iter(proc.stdout.readline, b''):
        print(stdout_line.decode("utf-8"))
    proc.stdout.close()
    return_code = proc.wait()


    for fn in glob(os.path.join(igsnper_dir, 'html_reports', '*.html')):
        gene = os.path.splitext(os.path.basename(fn))[0]
        if db.session.query(Gene).filter(Gene.name == gene).count() == 1:
            db.session.query(Gene).filter(Gene.name == gene).update({Gene.igsnper_plot_path: 'igsnper/html_reports/%s.html' % gene})
        else:
            print('Igsnpper identified gene %s, which is not listed in the Genes table.' % gene)

    for fn in glob(os.path.join(igsnper_dir, '*_processed/*.txt')):
        subject = os.path.splitext(os.path.basename(fn))[0]
        study, individual = subject.split('_')
        sample_name = subject + '_S1'

        sample_query = db.session.query(Sample).filter(Sample.name == sample_name)
        if sample_query.count() == 1:
            sample_query.update({Sample.igsnper_plot_path: 'igsnper/%s_processed/%s.txt' % (study, subject)})
            sample_id = sample_query.one_or_none().id
            db.session.query(Patient).filter(Patient.name == subject).update({Patient.igsnper_sample_id: sample_id})
        else:
            print('Cant find the record for sample %s' % sample_name)

    db.session.commit()

    db.close()


    return("IgSNP completed!")

