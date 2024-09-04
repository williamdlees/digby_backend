# Services related to vdjbase repseq-based data sets

import datetime
import decimal
import os.path
from os.path import isfile
import json
from math import ceil
import pickle

from flask import request
from flask_restx import Resource, reqparse

from api.reports.report_utils import make_output_file
from api.restx import api
from sqlalchemy import inspect, func, or_
from sqlalchemy import null as sa_null
from sqlalchemy_filters import apply_filters
from werkzeug.exceptions import BadRequest
from db.filter_list import apply_filter_to_list
from api.system.system import digby_protected

from app import vdjbase_dbs, app, genomic_dbs
from db.vdjbase_api_query_filters import sample_info_filters, sequence_filters
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype, Allele, AllelesSample, Gene, AlleleConfidenceReport, HaplotypeEvidence
from db.vdjbase_airr_model import GenoDetection, SeqProtocol, Study, TissuePro, Patient, Sample, DataPro
from db.genomic_db import Gene as GenomicGene

from api.reports.genotypes import process_repseq_genotype

VDJBASE_SAMPLE_PATH = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/samples')


# Return SqlAlchemy row as a dict, using correct column names
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


ns = api.namespace('repseq', description='Genes and annotations inferred from RepSeq data')

# for use by genomic api

def get_vdjbase_species():
    return list(vdjbase_dbs.keys())


def get_genomic_species():
    return list(genomic_dbs.keys())


@ns.route('/species')
class SpeciesApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the list of species for which information is held """
        return list(set(vdjbase_dbs.keys()) | set(get_genomic_species()))


def find_datasets(species):
    datasets = []
    for name in vdjbase_dbs[species].keys():
        datasets.append({'dataset': name, 'description': vdjbase_dbs[species][name].description})
    return datasets


@ns.route('/all_novels')
class NovelsApi(Resource):
    @digby_protected()
    def get(self):
        """ Returns the list all novel alleles across all datasets """
        ret = {}
        species = get_vdjbase_species()
        for sp in species:
            if sp in vdjbase_dbs:
                for ds_name, ds in vdjbase_dbs[sp].items():
                    session = ds.session
                    novels = session.query(Allele).filter(Allele.novel == 1).all()

                    if novels:
                        if sp not in ret:
                            ret[sp] = {}
                        if ds_name not in ret[sp]:
                            ret[sp][ds_name] = {}
                        for novel in novels:
                            ret[sp][ds_name][novel.name] = (novel.seq.replace('.', ''), novel.appears)

        return ret


@ns.route('/novels/<string:species>/<string:dataset>')
class NovelsSpApi(Resource):
    @digby_protected()
    def get(self, species, dataset):
        """ Returns details of full-length novel alleles in a single dataset"""
        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        novels = session.query(Allele)\
            .join(AllelesSample)\
            .join(Sample)\
            .join(SeqProtocol, Sample.seq_protocol_id == SeqProtocol.id)\
            .filter(Allele.novel == True)\
            .filter(Allele.is_single_allele == True)\
            .filter(or_(SeqProtocol.complete_sequences.ilike('full'), SeqProtocol.complete_sequences.ilike('complete')))\
            .all()

        results = {}
        for novel in novels:
            result = {
                'name': novel.name,
                'subject_count': len(set([sample.patient_id for sample in novel.samples])),
                'j_haplotypes': 0,
                'd_haplotypes': 0,
                'hetero_alleleic_j_haplotypes': 0,
                'example': '',
                'sequence': '',
            }

            haplos = session.query(HaplotypeEvidence)\
                .join(Sample)\
                .filter(HaplotypeEvidence.allele_id == novel.id)\
                .all()

            d_haps = 0
            j_haps = 0
            hetero_haps = 0
            best_hap = {'gene_type': '', 'hetero': False, 'count': 0, 'example': novel.samples[0].sample_name}

            for haplo in haplos:
                gene_type = 'D' if 'D' in haplo.hap_gene else 'J'
                if gene_type == 'D':
                    d_haps += 1
                else:
                    j_haps += 1

                counts = haplo.counts.split('),')
                hetero = len(counts) > 1

                if hetero and gene_type == 'J':
                    hetero_haps += 1

                novel_allele = novel.name.split('*')[1].upper()
                novel_count = 0

                try:
                    for count in counts:
                        if novel_allele + ' ' in count:
                            count = count.split('(')[1]
                            novel_counts = count.replace(')', '').replace(' ', '')
                            novel_count = max([int(c) for c in novel_counts.split(',')])
                except:
                    print('Error in parsing haplotype counts: %s' % haplo.counts)
                    novel_count = 0

                if best_hap['gene_type'] == '' \
                    or (best_hap['gene_type'] == 'D' and gene_type == 'J') \
                    or (best_hap['hetero'] and not hetero) \
                    or best_hap['count'] < novel_count:
                        best_hap = {'gene_type': gene_type, 'hetero': hetero, 'count': novel_count, 'example': haplo.sample.sample_name}

            result['j_haplotypes'] = j_haps
            result['d_haplotypes'] = d_haps
            result['hetero_alleleic_j_haplotypes'] = hetero_haps
            result['example'] = best_hap['example']
            result['sequence'] = novel.seq

            results[result['name']] = result

        return results


@ns.route('/ref_seqs/<string:species>')
class DataSetAPI(Resource):
    @digby_protected()
    def get(self, species):
        """ Returns the list of datasets available for the selected species """
        if species in vdjbase_dbs:
            return find_datasets(species)
        else:
            return list()


@ns.route('/dataset_info/<string:species>/<string:dataset>')
class DataSetInfoAPI(Resource):
    @digby_protected()
    def get(self, species, dataset):
        """Returns information and statistics on the dataset"""
        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        stats = {}

        stats['description'] = vdjbase_dbs[species][dataset].description
        stats['total_subjects'] = session.query(Patient.id).count()
        stats['total_samples'] = session.query(Sample.id).count()
        stats['sex_count'] = session.query(Patient.sex, func.count(Patient.sex)).group_by(Patient.sex).order_by(func.count(Patient.sex).desc()).all()
        stats['study_count'] = session.query(Study.study_name, func.count(Sample.sample_name)).join(Sample, Sample.study_id == Study.id).group_by(Study.study_name).order_by(func.count(Sample.sample_name).desc()).all()
        stats['condition_count'] = session.query(Patient.disease_diagnosis_label, func.count(Patient.disease_diagnosis_label)).group_by(Patient.disease_diagnosis_label).order_by(func.count(Patient.disease_diagnosis_label).desc()).all()
        stats['celltype_count'] = session.query(TissuePro.sub_cell_type, func.count(TissuePro.sub_cell_type)).group_by(TissuePro.sub_cell_type).order_by(func.count(TissuePro.sub_cell_type).desc()).all()
        stats['tissue_count'] = session.query(TissuePro.tissue_label, func.count(TissuePro.tissue_label)).group_by(TissuePro.tissue_label).order_by(func.count(TissuePro.tissue_label).desc()).all()
        studies = session.query(
            Study.study_name,
            Study.study_id,
            Study.lab_address,
            Study.num_subjects,
            Study.num_samples,
            Study.pub_ids,
            Study.accession_reference
        ).all()
        stats['studies'] = []

        for row in studies:
            row = row._asdict()
            row['subjects_in_vdjbase'] = session.query(Study.id).join(Patient, Patient.study_id == Study.id).filter(Study.study_name == row['study_name']).count()
            row['samples_in_vdjbase'] = session.query(Study.id).join(Sample, Sample.study_id == Study.id).filter(Study.study_name == row['study_name']).count()
            row['study_id'] = link_convert(row['study_id'])
            row['pub_ids'] = link_convert(row['pub_ids'])
            stats['studies'].append(row)

        stats['studies'].sort(key = lambda p: int(p['study_name'].replace('P', '')))

        return stats


# Convert some frequently occurring accession ids to links
def link_convert(item):
    res = []

    if not item:
        return None

    for el in item.split(','):
        for prefix in ['BioProject:', 'SRA:']:
            if prefix in el:
                el = el.replace(prefix, '')
        el = el.strip()

        if 'PMID:' in el:
            el = el.replace('PMID:', '')
            el = el.strip()
            el = f'https://pubmed.ncbi.nlm.nih.gov/{el}'
        elif ' ' not in el:   # general safety measure for things that won't translate
            if 'PRJNA' in el:
                el = f'https://www.ncbi.nlm.nih.gov/bioproject/?term={el}'
            elif 'PRJEB' in el:
                el = f'https://www.ebi.ac.uk/ena/browser/view/{el}'
            elif 'PRJCA' in el:
                el = f'https://ngdc.cncb.ac.cn/bioproject/browse/{el}'
            elif 'SRP' in el:
                el = f'https://www.ncbi.nlm.nih.gov/sra/?term={el}'

        res.append(el)

    return ','.join(res)




rep_sample_bool_values = {
    'synthetic': ('Synthetic', '(blank)'),
    'single_assignment': ('Single', '(blank)'),
}

@ns.route('/sample_info/<string:species>/<string:dataset>/<string:sample>')
class SampleInfoApi(Resource):
    @digby_protected()
    def get(self, species, dataset, sample):
        """ Returns information on the selected sample """

        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        return get_sample_info(species, dataset, sample), 200


@ns.route('/all_samples_info/<string:species>/<string:dataset>')
class AllSamplesInfoApi(Resource):
    @digby_protected()
    def get(self, species, dataset):
        """ Returns information on all samples """
        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        cache_filename = f"{app.config['OUTPUT_PATH']}/airrseq_all_samples_info_{species}_{dataset}.pickle"
        if os.path.isfile(cache_filename):
            # check that the file is newer than the last revision dates of the databases
            last_revision_time = max([genomic_dbs[species][dataset].created for dataset in genomic_dbs[species].keys()])
            if datetime.datetime.fromtimestamp(os.path.getmtime(cache_filename)) > last_revision_time:
                try:
                    with open(cache_filename, 'rb') as f:
                        return pickle.load(f), 200
                except Exception as e:
                    app.logger.error(f"Error loading cache file {cache_filename}: {e}")
                    os.remove(cache_filename)
            else:
                os.remove(cache_filename)

        metadata_list = []

        for dataset in vdjbase_dbs[species].keys():
            session = vdjbase_dbs[species][dataset].session
            samples = session.query(Sample.sample_name).all()
            for sample in samples:
                metadata_list.append(get_sample_info(species, dataset, sample[0]))

        if not metadata_list:
            return None, 404

        with open(cache_filename, 'wb') as f:
            pickle.dump(metadata_list, f, pickle.HIGHEST_PROTOCOL)

        return metadata_list, 200


def get_sample_info(species, dataset, sample):
    session = vdjbase_dbs[species][dataset].session
    attribute_query = []

    for col in sample_info_filters.keys():
        if sample_info_filters[col]['field'] is not None:
            attribute_query.append(sample_info_filters[col]['field'])

    info = session.query(*attribute_query)\
        .join(GenoDetection, GenoDetection.id == Sample.geno_detection_id)\
        .join(Patient, Patient.id == Sample.patient_id)\
        .join(SeqProtocol, SeqProtocol.id == Sample.seq_protocol_id)\
        .join(TissuePro, TissuePro.id == Sample.tissue_pro_id)\
        .join(DataPro, DataPro.id == Sample.data_pro_id) \
        .join(Study, Sample.study_id == Study.id)\
        .filter(Sample.sample_name == sample).one_or_none()

    if info:
        info = info._asdict()

        for k,v in info.items():
            if v:
                if isinstance(v, (datetime.datetime, datetime.date)):
                    info[k] = v.isoformat()

        haplotypes = session.query(HaplotypesFile.by_gene_s).join(SamplesHaplotype).join(Sample).filter(Sample.sample_name==sample).order_by(HaplotypesFile.by_gene_s).all()
        info['haplotypes'] = [(h[0]) for h in haplotypes]

    return info


filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int, location='args')
filter_arguments.add_argument('page_size', type=int, location='args')
filter_arguments.add_argument('filter', type=str, location='args')
filter_arguments.add_argument('sort_by', type=str, location='args')
filter_arguments.add_argument('cols', type=str, location='args')


@ns.route('/samples/<string:species>/<string:dataset>')
class SamplesApi(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of samples in the selected dataset """

        if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
            return list(), 404

        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols']) if args['cols'] else list(sample_info_filters.keys())
        if 'sample_name' not in required_cols:
            required_cols.append('sample_name')
        if 'dataset' not in required_cols:
            required_cols.append('dataset')

        for col in required_cols:
            if col not in sample_info_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

        if 'genotype' in required_cols:                # needed to compose paths to files
            for field in ('patient_name', 'study_name'):
                if field not in required_cols:
                    required_cols.append(field)
            required_cols.append('genotypes')
            required_cols.append('genotype_stats')
            required_cols.append('genotype_report')

        attribute_query = [sample_info_filters['sample_id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'sample_id' and sample_info_filters[col]['field'] is not None:
                attribute_query.append(sample_info_filters[col]['field'])

        attribute_query.append(Sample.id)

        filter = json.loads(args['filter']) if args['filter'] else []
        ret = find_vdjbase_samples(attribute_query, species, dataset.split(','), filter)

        total_size = len(ret)

        uniques = {}

        for f in required_cols:
            uniques[f] = []

        uniques['dataset'] = dataset.split(',')

        # special column for names by dataset

        uniques['names_by_dataset'] = {}
        filter_applied = len(filter) > 0
        if filter_applied:
            for dataset in uniques['dataset']:
                uniques['names_by_dataset'][dataset] = []

        for s in ret:
            for f in required_cols:
                if sample_info_filters[f]['field'] is not None and 'no_uniques' not in sample_info_filters[f]:
                    el = s[f]
                    if isinstance(el, (datetime.datetime, datetime.date)):
                        el = el.date().isoformat()
                    elif isinstance(el, bool):
                        if f in rep_sample_bool_values:
                            el = rep_sample_bool_values[f][0 if el else 1]
                    if (not isinstance(el, str) or len(el) > 0) and el not in uniques[f]:
                        uniques[f].append(el)
                    if (sample_info_filters[f]['field'].type.python_type is str) and len(el) == 0 and '(blank)' not in uniques[f]:
                        uniques[f].append('(blank)')

            if filter_applied:
                uniques['names_by_dataset'][s['dataset']].append(s['sample_name'])

        def name_sort_key(name):
            name = name.split('_')
            for i in range(len(name)):
                name[i] = name[i][1:].zfill(4)
            return name

        def num_sort_key(x):
            if x is None or x == '':
                return -1
            else:
                try:
                    return float(x)
                except:
                    return 0

        for f in required_cols:
            try:
                if 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'underscore':
                    uniques[f].sort(key=name_sort_key)
                elif 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'numeric':
                    uniques[f].sort(key=num_sort_key)
                elif f != 'names_by_dataset':
                    uniques[f].sort(key=lambda x: (x is None or x == '', x))
            except:
                pass

        if 'haplotypes' in required_cols:
            uniques['haplotypes'] = []
            for dset in dataset.split(','):
                session = vdjbase_dbs[species][dset].session
                haplotypes = session.query(HaplotypesFile.by_gene_s).distinct().order_by(HaplotypesFile.by_gene_s).all()
                x = [(h[0]) for h in haplotypes]
                uniques['haplotypes'].extend(x)
            uniques['haplotypes'] = list(set(uniques['haplotypes']))

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None) else []
        if len(sort_specs) == 0:
            sort_specs = [{'field': 'name', 'order': 'asc'}]

        for spec in sort_specs:
            f = spec['field']
            if f in sample_info_filters.keys():
                if 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'underscore':
                    ret = sorted(ret, key=lambda x: name_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
                elif 'sort' in sample_info_filters[f] and sample_info_filters[f]['sort'] == 'numeric':
                    ret = sorted(ret, key=lambda x: num_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
                else:
                    ret = sorted(ret, key=lambda x: ((x[f] is None or x[f] == ''), x[f]), reverse=(spec['order'] == 'desc'))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first:first + args['page_size']]

        for rec in ret:
            for k, v in rec.items():
                if isinstance(v, (datetime.datetime, datetime.date)):
                    rec[k] = v.date().isoformat()
                elif isinstance(v, decimal.Decimal):
                    rec[k] = '%0.2f' % v

        if 'genotypes' in required_cols:
            for r in ret:
                r['genotypes'] = {}
                r['genotypes']['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['sample_name'], 'sort_order': 'Locus'})

                r['genotypes']['path'] = app.config['BACKEND_LINK']
                sp = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset']]) + '/'
                r['genotypes']['tigger'] = sp + r['genotype'].replace('samples', '') if r['genotype_stats'] else ''
                r['genotypes']['ogrdbstats'] = sp + r['genotype_stats'].replace('samples', '') if r['genotype_stats'] else ''
                r['genotypes']['ogrdbplot'] = sp + r['genotype_report'].replace('samples', '') if r['genotype_report'] else ''
                del r['genotype_stats']
                del r['genotype_report']

                session = vdjbase_dbs[species][r['dataset']].session
                igsnper_path = session.query(Sample.igsnper_plot_path, Sample.genotype_report).filter(Sample.sample_name == r['sample_name']).one_or_none()

                if igsnper_path is not None and igsnper_path[0] is not None:
                    r['genotypes']['igsnper'] = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], igsnper_path[0]])
                else:
                    r['genotypes']['igsnper'] = ''


        if 'haplotypes' in required_cols:
            for r in ret:
                session = vdjbase_dbs[species][r['dataset']].session
                haplotypes = session.query(Sample.sample_name, func.group_concat(HaplotypesFile.by_gene_s), func.group_concat(HaplotypesFile.file))
                h = haplotypes.filter(Sample.sample_name == r['sample_name'])\
                    .join(SamplesHaplotype, SamplesHaplotype.samples_id == Sample.id)\
                    .join(HaplotypesFile, HaplotypesFile.id == SamplesHaplotype.haplotypes_file_id)\
                    .one_or_none()
                if h is not None and h[1] is not None:
                    r['haplotypes'] = {}
                    r['haplotypes']['path'] = app.config['BACKEND_LINK']
                    for (hap, filename) in zip(h[1].split(','), h[2].split(',')):
                        filename = filename.replace('samples/', '')
                        fp = os.path.join(VDJBASE_SAMPLE_PATH, species, r['dataset'], filename)
                        sp = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], filename])
                        if isfile(fp):
                            r['haplotypes'][hap] = {}
                            r['haplotypes'][hap]['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['sample_name'], 'hap_gene': hap, 'sort_order' : 'Locus'})
                            r['haplotypes'][hap]['rabhit'] = sp
                else:
                    r['haplotypes'] = ''

        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
        }, 200

def find_vdjbase_samples(attribute_query, species, datasets, filter):
    hap_filters = None
    allele_filters = None
    dataset_filters = []
    filter_spec = []

    for f in list(filter):
        try:
            if f['field'] == 'haplotypes':
                hap_filters = f
            elif f['field'] == 'allele':
                allele_filters = f
            elif f['field'] == 'dataset':
                dataset_filters.append(f)
            else:
                f['model'] = sample_info_filters[f['field']]['model']
                if 'fieldname' in sample_info_filters[f['field']]:
                    f['field'] = sample_info_filters[f['field']]['fieldname']
                if f['field'] in rep_sample_bool_values:
                    value = []
                    for v in f['value']:
                        value.append('1' if v == rep_sample_bool_values[f['field']][0] else '0')
                    f['value'] = value
                elif '(blank)' in f['value']:
                    value_specs = [
                        {'model': sample_info_filters[f['field']]['model'], 'field': f['field'], 'op': 'is_null', 'value': ''},
                        {'model': sample_info_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': ''},
                    ]

                    for v in f['value']:
                        if v != '(blank)':
                            value_specs.append({'model': sample_info_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': v})
                    
                    f = {'or': value_specs}

                filter_spec.append(f)

        except Exception as e:
            raise BadRequest(f'Bad filter string: {f}: {e}')
    ret = []

    if len(dataset_filters) > 0:
        apply_filter_to_list(datasets, dataset_filters)

    for dset in datasets:
        session = vdjbase_dbs[species][dset].session

        query = session.query(*attribute_query)\
            .join(GenoDetection, Sample.geno_detection_id == GenoDetection.id)\
            .join(Patient, Sample.patient_id == Patient.id)\
            .join(SeqProtocol, Sample.seq_protocol_id == SeqProtocol.id)\
            .join(TissuePro, Sample.tissue_pro_id == TissuePro.id)\
            .join(Study, Sample.study_id == Study.id)
        query = apply_filters(query, filter_spec)

        if hap_filters:
            hap_samples = session.query(Sample.sample_name.distinct()).join(SamplesHaplotype).join(HaplotypesFile).filter(
                HaplotypesFile.by_gene_s.in_(hap_filters['value']))
            query = query.filter(Sample.sample_name.in_(hap_samples))

        if allele_filters:
            allele_samples = session.query(Sample.sample_name.distinct()).join(AllelesSample,
                                                                        Sample.id == AllelesSample.sample_id).join(
                Allele, Allele.id == AllelesSample.allele_id).filter(Allele.name.in_(allele_filters['value'])).all()
            if allele_samples is None:
                allele_samples = []
            query = query.filter(Sample.sample_name.in_([s[0] for s in allele_samples]))

        res = query.all()

        for r in res:
            s = r._asdict()
            for k, v in s.items():
                if isinstance(v, (datetime.datetime, datetime.date)):
                    s[k] = v.date().isoformat()
                if v is None:
                    s[k] = ''
            s['dataset'] = dset
            s['id'] = '%s.%d' % (dset, s['id'])
            ret.append(s)
    return ret


rep_sequence_bool_values = {
    'is_single_allele': ('Unambiguous', '(blank)'),
    'low_confidence': ('Low Confidence', '(blank)'),
    'novel': ('Novel', '(blank)'),
}

@ns.route('/sequences/<string:species>/<string:dataset>')
class SequencesApi(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of sequences in the selected datasets """

        if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
            return list()

        args = filter_arguments.parse_args(request)
        if (args['page_size'] and args['page_size'] <= 0) or (args['page_number'] and args['page_number'] < 0):
            return list()

        if args['cols'] is None or not args['cols']:
            required_cols = list(sequence_filters.keys())
        else:
            required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in sequence_filters.keys():
                print('bad column in request: %s' % col)
                return list(), 404

        if 'gene_name' in required_cols:
            required_cols.append('igsnper_plot_path')

        datasets = dataset.split(',')
        ret = find_vdjbase_sequences(species, datasets, required_cols, json.loads(args['filter']) if args['filter'] else [])
        total_size = len(ret)

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in ret:
            for f in required_cols:
                if sequence_filters[f]['field'] is not None and 'no_uniques' not in sequence_filters[f]:
                    el = s[f]
                    if isinstance(el, (datetime.datetime, datetime.date)):        # convert and add to uniques. Can't convert the returned records until we sort
                        el = el.date().isoformat()
                    elif isinstance(el, decimal.Decimal):
                        el = '%0.2f' % el
                    elif isinstance(el, bool):
                        if f in rep_sequence_bool_values:
                            el = rep_sequence_bool_values[f][0 if el else 1]
                    if (not isinstance(el, str) or len(el) > 0) and el not in uniques[f]:
                        uniques[f].append(el)
                    if isinstance(el, str) and len(el) == 0 and '(blank)' not in uniques[f]:
                        uniques[f].append('(blank)')

        uniques['dataset'] = dataset.split(',')

        # For gene order, it would be a good idea if datasets from the same species agreed!
        # really need some way of merging gene order from several sets

        session = vdjbase_dbs[species][datasets[0]].session
        gene_order = session.query(Gene.name, Gene.alpha_order).all()
        gene_order = {x[0]: x[1] for x in gene_order}

        def allele_sort_key(name):
            if '*' in name:
                gene = name.split('*')
            else:
                gene = (name, '')

            return((gene_order[gene[0]] if gene[0] in gene_order else 999, gene[1]))

        def num_sort_key(x):
            if x is None or x == '':
                return -1
            else:
                try:
                    return float(x)
                except:
                    return 0

        for f in required_cols:
            try:
                if 'sort' in sequence_filters[f] and sequence_filters[f]['sort'] == 'allele':
                    uniques[f].sort(key=allele_sort_key)
                elif 'sort' in sequence_filters[f] and sequence_filters[f]['sort'] == 'numeric':
                    uniques[f].sort(key=num_sort_key)
                else:
                    uniques[f] = list(set([x for x in sorted(uniques[f], key=lambda x: (x is None or x == '', x))]))
            except:
                pass

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        if len(sort_specs) > 0:
            for spec in sort_specs:
                if spec['field'] in sequence_filters.keys():
                    f = spec['field']
                    if 'sort' in sequence_filters[f] and sequence_filters[f]['sort'] == 'numeric':
                        ret = sorted(ret, key=lambda x: num_sort_key(x[f]), reverse=(spec['order'] == 'desc'))
                    else:
                        ret = sorted(ret, key=lambda x: ((x[f] is None or x[f] == ''), x[f]), reverse=(spec['order'] == 'desc'))
        else:
            ret = sorted(ret, key=lambda x : allele_sort_key(x['name']))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first : first + args['page_size']]

        for rec in ret:
            for k, v in rec.items():
                if isinstance(v, (datetime.datetime, datetime.date)):
                    rec[k] = v.date().isoformat()
                elif isinstance(v, decimal.Decimal):
                    rec[k] = '%0.2f' % v



        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size']) if args['page_size'] else 1
        }


def find_vdjbase_sequences(species, datasets, required_cols, seq_filter):
    sample_id_filter = None
    filter_spec = []
    dataset_filters = []
    if seq_filter:
        for f in seq_filter:
            try:
                if f['field'] == 'sample_id':
                    sample_id_filter = f
                elif f['field'] == 'dataset':
                    #                        f['model'] = 'dataset_table'
                    dataset_filters.append(f)
                else:
                    if f['field'] == 'similar':
                        value = []
                        for v in f['value']:
                            if v != '(blank)':
                                value.append('|' + v.replace(',', '|') + '|')
                            else:
                                value.append(v)
                        f['value'] = value
                    f['model'] = sequence_filters[f['field']]['model']
                    if 'fieldname' in sequence_filters[f['field']]:
                        f['field'] = sequence_filters[f['field']]['fieldname']
                    if f['field'] in rep_sequence_bool_values:
                        value = []
                        for v in f['value']:
                            value.append('1' if v == rep_sequence_bool_values[f['field']][0] else '0')
                        f['value'] = value
                    elif '(blank)' in f['value']:
                        value_specs = [
                            {'model': sequence_filters[f['field']]['model'], 'field': f['field'], 'op': 'is_null', 'value': ''},
                            {'model': sequence_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': ''},
                        ]

                        for v in f['value']:
                            if v != '(blank)':
                                value_specs.append({'model': sequence_filters[f['field']]['model'], 'field': f['field'], 'op': '==', 'value': v})
                        
                        f = {'or': value_specs}

                    filter_spec.append(f)
            except Exception as e:
                raise BadRequest(f'Bad filter string: {f}: {e}')

    if 'notes_count' in required_cols and 'notes' not in required_cols:
        required_cols.append('notes')
    ret = []
    if len(dataset_filters) > 0:
        apply_filter_to_list(datasets, dataset_filters)
    for dset in datasets:
        session = vdjbase_dbs[species][dset].session

        attribute_query = []

        for col in required_cols:
            if col not in sequence_filters or 'field' not in sequence_filters[col]:
                breakpoint()
            if sequence_filters[col]['field'] is not None:
                attribute_query.append(sequence_filters[col]['field'])

        query = session.query(*attribute_query).join(Gene, Allele.gene_id == Gene.id)
        query = apply_filters(query, filter_spec)

        if 'notes' in required_cols:
            query = query.outerjoin(AlleleConfidenceReport, Allele.id == AlleleConfidenceReport.allele_id).group_by(
                Allele.id)

        res = query.all()

        required_names = []
        appears = {}

        if sample_id_filter is not None:
            required_ids = []
            if dset in sample_id_filter['value']:
                for id in sample_id_filter['value'][dset]:
                    required_ids.append(id)

            alleles_with_samples = session.query(Allele) \
                .join(AllelesSample) \
                .join(Sample) \
                .filter(Sample.sample_name.in_(required_ids)).all()

            for a in alleles_with_samples:
                appearances = session.query(AllelesSample.patient_id) \
                    .filter(AllelesSample.hap == 'geno') \
                    .filter(AllelesSample.allele_id == a.id) \
                    .join(Sample) \
                    .filter(Sample.sample_name.in_(required_ids)) \
                    .filter(AllelesSample.hap == 'geno') \
                    .distinct().count()

                if a.similar is not None and a.similar != '':
                    sims = a.similar.split(', ')
                    for sim in sims:
                        sim = sim.replace('|', '')
                        appearances += session.query(AllelesSample.patient_id) \
                            .join(Allele) \
                            .filter(AllelesSample.hap == 'geno') \
                            .filter(Allele.name.ilike(sim)) \
                            .distinct().count()

                required_names.append(a.name)
                appears[a.name] = appearances

        for r in res:
            if len(required_names) == 0 or r.name in required_names:
                s = r._asdict()

                if len(required_names) > 0:
                    s['appears'] = appears[r.name]

                for k, v in s.items():
                    if k == 'similar' and v is not None:
                        s[k] = v.replace('|', '')
                s['dataset'] = dset

                if 'igsnper_plot_path' in s and s['igsnper_plot_path'] is not None and len(s['igsnper_plot_path']) > 0:
                    s['igsnper_plot_path'] = '/'.join(
                        [app.config['BACKEND_LINK'], 'static/study_data/VDJbase/samples', species, dset,
                         s['igsnper_plot_path']])
                else:
                    s['igsnper_plot_path'] = ''

                ret.append(s)

    return ret


@ns.route('/all_subjects_genotype/<string:species>')
class AllSubjectsGenotypeApi(Resource):
    @digby_protected()
    def get(self, species):
        """ Return genotypes for all subjects of the specified species in the specified data type """

        if species not in vdjbase_dbs:
            return None, 404
          
        cache_filename = f"{app.config['OUTPUT_PATH']}/airrseq_all_subjects_genotype_{species}.pickle"
        if os.path.isfile(cache_filename):
            # check that the file is newer than the last revision dates of the databases
            last_revision_time = max([genomic_dbs[species][dataset].created for dataset in genomic_dbs[species].keys()])
            if datetime.datetime.fromtimestamp(os.path.getmtime(cache_filename)) > last_revision_time:
                try:
                    with open(cache_filename, 'rb') as f:
                        return pickle.load(f), 200
                except Exception as e:
                    app.logger.error(f"Error loading cache file {cache_filename}: {e}")
                    os.remove(cache_filename)
            else:
                os.remove(cache_filename)

        all_subjects = []
        for dataset in vdjbase_dbs[species].keys():
            session = vdjbase_dbs[species][dataset].session
            subjects = session.query(Patient.patient_name).all()
            all_subjects.extend(subjects)

        all_subjects = sorted(list(set([s[0] for s in all_subjects])))

        genotype_sets = []

        for subject_name in all_subjects:
            genotypes = []
            for dataset in vdjbase_dbs[species].keys():
                genotype = single_genotype(species, dataset, subject_name)
                if genotype:
                    genotypes.append(genotype)
            if genotypes:
                genotype_sets.append({
                    'subject_name': subject_name,
                    'GenotypeSet': {
                        'receptor_genotype_set_id': 'Genomic_genotype_set_' + subject_name,
                        'genotype_class_list': genotypes
                    }
                })

        if not genotype_sets:
            return None, 404

        with open(cache_filename, 'wb') as f:
            pickle.dump(genotype_sets, f, pickle.HIGHEST_PROTOCOL)

        return genotype_sets, 200


@ns.route('/genotype/<string:species>/<string:subject_name>')
class GenotypeApi(Resource):
    @digby_protected()
    def get(self, species, subject_name):
        """ Returns the inferred genotype (in MiAIRR format) of the specified sample """

        if species not in vdjbase_dbs:
            return None, 404

        genotypes = []

        for dataset in vdjbase_dbs[species].keys():
            genotype = single_genotype(species, dataset, subject_name)
            if genotype:
                genotypes.append(genotype)

        if not genotypes:
            return None, 404

        ret = {
            'GenotypeSet': {
                'receptor_genotype_set_id': 'Tigger_genotype_set_' + subject_name,
                'genotype_class_list': genotypes
            }

        }

        return ret, 400


def single_genotype(species, dataset, subject_name):
    session = vdjbase_dbs[species][dataset].session
    samples = session.query(Sample).join(Patient, Sample.patient_id == Patient.id).filter(Patient.patient_name == subject_name).all()

    if len(samples) == 0:
        return None

    sample = samples[0]     # TODO more intelligent way to select sample??

    genotype = process_repseq_genotype(sample.sample_name, [], session, False)
    germline_set = {
        'V': sample.geno_detection.aligner_reference_v,
        'D': sample.geno_detection.aligner_reference_d,
        'J': sample.geno_detection.aligner_reference_j,
    }
    documented = []
    undocumented = []
    deleted = []
    for row in genotype.itertuples():
        gene_type = row.gene[3]

        if gene_type not in germline_set.keys():
            continue

        if row.alleles == 'Del':
            deleted.append({'label': row.gene, 'germline_set_ref': germline_set[gene_type], 'phasing': 0})
        for allele in row.GENOTYPED_ALLELES.split(','):
            allele_name = row.gene + '*' + allele
            res = session.query(Allele.seq, Allele.novel).filter(Allele.name == allele_name).one_or_none()
            if res:
                seq, novel = res

                if novel:
                    undocumented.append({'allele_name': allele_name, 'germline_set_ref': germline_set[gene_type], 'sequence': seq, 'phasing': 0})
                else:
                    documented.append({'label': allele_name, 'germline_set_ref': germline_set[gene_type], 'phasing': 0})
    ret = {
        'receptor_genotype_id': 'Tigger_genotype_' + sample.sample_name + '_' + dataset,
        'locus': dataset,
        'documented_alleles': documented,
        'undocumented_alleles': undocumented,
        'deleted_genes': deleted,
        'inference_process': 'repertoire_sequencing',
        'genotyping_tool': sample.geno_detection.geno_tool,
        'genotyping_tool_version': sample.geno_detection.geno_ver,
    }
    return ret


def find_rep_filter_params(species, datasets):
    if species not in vdjbase_dbs or set(datasets).difference(set(vdjbase_dbs[species])):
        return ([], [])

    genes = []
    gene_types = []
    haplotypes = []

    for dataset in datasets:
        session = vdjbase_dbs[species][dataset].session
        g_q = session.query(Gene.name, Gene.type).all()
        genes.extend((res[0] for res in g_q))
        gene_types.extend((res[1] for res in g_q))
        h_q = session.query(HaplotypesFile.by_gene).distinct()
        haplotypes.extend(res[0] for res in h_q)

    genes = sorted(set(genes))
    gene_types = sorted(set(gene_types))
    haplotypes = sorted(set(haplotypes))

    params = [
        {
          "id": "f_pseudo_genes",
          "type": "boolean",
          "label": "Include pseudogenes"
        },
    ]

    params.extend([
        {
            "id": "f_gene_types",
            "type": "multi_select",
            "label": "Only process selected gene types (leave blank for all)",
            "options": gene_types
        },
        {
            "id": "f_genes",
            "type": "multi_select",
            "label": "Only process selected genes (leave blank for all)",
            "options": genes

        }
    ])

    return (params, haplotypes)

# Test if session refers to an AIRR-seq database

def is_session_airrseq(session):
    for species in vdjbase_dbs.keys():
        for dataset in vdjbase_dbs[species].keys():
            if vdjbase_dbs[species][dataset] == session:
                return True

# Apply filter params to a list of samples in the context of a specific dataset
# wanted_genes is returned in the required search order

def apply_rep_filter_params(params, sample_list, session):
    if is_session_airrseq(session):
        gene = Gene
    else:
        gene = GenomicGene

    if 'per_sample' in params:
        sample_list = filter_per_sample(params['per_sample'], sample_list)
    if 'sort_order' in params:
        if params['sort_order'] == 'Alphabetic':
            gq = session.query(gene).order_by(gene.alpha_order)
        else:
            gq = session.query(gene).order_by(gene.locus_order)
    else:
        gq = session.query(gene).order_by(gene.alpha_order)

    if len(params['f_gene_types']) > 0:
        gq = gq.filter(gene.type.in_(params['f_gene_types']))
    if len(params['f_genes']) > 0:
        gq = gq.filter(gene.name.in_(params['f_genes']))
    if not params['f_pseudo_genes']:
        gq = gq.filter(gene.pseudo_gene == 0)
    wanted_genes = gq.all()
    wanted_genes = [g.name for g in wanted_genes]
    return sample_list, wanted_genes


def filter_per_sample(per_sample, sample_list):
    if (per_sample != 'Each sample'):
        pp_list = []
        ids = []
        for (name, genotype, patient_id) in sample_list:
            if patient_id not in ids:
                ids.append(patient_id)
                pp_list.append([name, genotype, patient_id])
        sample_list = pp_list
    return sample_list


def get_order_file(species, dataset, locus_order=True, genomic=False):
    file_name = os.path.join(app.config['OUTPUT_PATH'], '%s_%s_%s_order.tsv' % (species, dataset, 'locus' if locus_order else 'alpha'))

    if genomic:
        session = genomic_dbs[species][dataset].session
        gene = GenomicGene
    else:
        session = vdjbase_dbs[species][dataset].session
        gene = Gene

    if locus_order:
        gene_order = session.query(gene.name).order_by(gene.locus_order).all()
    else:
        gene_order = session.query(gene.name).order_by(gene.alpha_order).all()
    gene_order = [x[0] for x in gene_order]

    with open(file_name, 'w') as fo:
        fo.write('\n'.join(gene_order))

    return file_name


# interleave content from multiple files, or return one of them if they are all the same
def get_multiple_order_file(species, rep_datasets, gen_datasets, locus_order=True):
    gene_order = []
    added = False

    for dataset in rep_datasets:
        added, gene_order = merge_order(added, dataset, gene_order, locus_order, species, False)

    for dataset in gen_datasets:
        added, gene_order = merge_order(added, dataset, gene_order, locus_order, species, True)

    if added:
        file_name = make_output_file('tsv')

        with open(file_name, 'w') as fo:
            fo.write('\n'.join(gene_order))
    else:
        if rep_datasets:
            file_name = get_order_file(species, list(rep_datasets)[0], locus_order, False)
        else:
            print(list(gen_datasets))
            file_name = get_order_file(species, list(gen_datasets)[0], locus_order, True)

    return file_name


def merge_order(added, dataset, gene_order, locus_order, species, genomic):
    with open(get_order_file(species, dataset, locus_order, genomic), 'r') as fi:
        if len(gene_order) == 0:
            gene_order = fi.read().split('\n')
        else:
            this_order = fi.read().split('\n')
            prev = ''

            for gene in this_order:
                if gene not in gene_order:
                    if prev in gene_order:
                        gene_order.insert(gene_order.index(prev), gene)
                    else:
                        gene_order.append(gene)
                    added = True
    return added, gene_order






