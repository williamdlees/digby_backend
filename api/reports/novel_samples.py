# Novel allele support in AIRR-seq repertoires

from werkzeug.exceptions import BadRequest
from sqlalchemy import or_
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file, chunk_list
from app import vdjbase_dbs
from db.vdjbase_model import AllelesSample, Allele, HaplotypeEvidence, Gene
from db.vdjbase_airr_model import Sample, SeqProtocol
from receptor_utils.simple_bio_seq import write_csv
from api.vdjbase.vdjbase import apply_rep_filter_params


SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if format != 'xls':
        raise BadRequest('Invalid format requested')

    rep_samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in rep_samples_by_dataset:
            rep_samples_by_dataset[rep_sample['dataset']] = []
        rep_samples_by_dataset[rep_sample['dataset']].append(rep_sample['sample_name'])

    results = []

    for dataset, sample_list in rep_samples_by_dataset.items():
        session = vdjbase_dbs[species][dataset].get_session()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        novels = session.query(Allele, AllelesSample, Sample)\
            .join(AllelesSample, AllelesSample.allele_id == Allele.id)\
            .join(Sample, AllelesSample.sample_id == Sample.id)\
            .join(SeqProtocol, Sample.seq_protocol_id == SeqProtocol.id)\
            .join(Gene, Allele.gene_id == Gene.id)\
            .filter(Allele.novel == True)\
            .filter(Gene.name.in_(wanted_genes))\
            .filter(Sample.sample_name.in_(sample_list))

        
        if 'Only' in params['read_length']:
            novels = novels.filter(or_(SeqProtocol.complete_sequences.ilike('full'), SeqProtocol.complete_sequences.ilike('complete')))

        novels = novels.all()

        for novel, allelesample, sample in novels:
            result = {}
            result['name'] = novel.name
            result['sample'] = sample.sample_name
            result['novel_count'] = allelesample.count
            result['total_count'] = allelesample.total_count
            results.append(result)

    output_path = make_output_file('csv')
    write_csv(output_path, results)
    return send_report(output_path, 'csv', f'{species}_novels_in_samples.csv')


