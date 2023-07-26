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
        session = vdjbase_dbs[species][dataset].session
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)

        novels = session.query(Allele)\
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

        for novel in novels:
            haplos = session.query(HaplotypeEvidence)\
                .join(Sample)\
                .filter(HaplotypeEvidence.allele_id == novel.id)\
                .filter(HaplotypeEvidence.sample_id == Sample.id)\
                .filter(Sample.sample_name.in_(sample_list))\
                .all()

            for haplo in haplos:
                counts = haplo.counts.split('),')
                novel_allele = novel.name.split('*')[1].upper()
                novel_count = 0
                total_count = 0

                try:
                    for count in counts:
                        acount = count.split('(')[1]
                        acount = acount.replace(')', '').replace(' ', '')
                        max_count = max([int(c) for c in acount.split(',')])
                        if novel_allele + ' ' in count:
                            novel_count = max_count
                        total_count += max_count

                except:
                    print('Error in parsing haplotype counts: %s' % haplo.counts)
                    novel_count = 0

                result = {}
                result['name'] = novel.name
                result['sample'] = haplo.sample.sample_name
                result['hap_gene'] = haplo.hap_gene
                result['counts'] = haplo.counts
                result['novel_count'] = novel_count
                result['total_count'] = total_count
                results.append(result)

    output_path = make_output_file('csv')
    write_csv(output_path, results)
    return send_report(output_path, 'csv', f'{species}_allele_support.csv')


