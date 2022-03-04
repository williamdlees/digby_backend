# Services related to vdjbase repseq-based data sets

from flask import request
from flask_restx import Resource, reqparse, fields, marshal, inputs

from api.reports.report_utils import make_output_file
from api.restx import api
from sqlalchemy import inspect, func, cast, literal, String, select, union_all
from math import ceil
import json
from sqlalchemy_filters import apply_filters
from werkzeug.exceptions import BadRequest
from datetime import datetime
import decimal
import os.path
from os.path import isfile
from db.filter_list import apply_filter_to_list
from api.system.system import digby_protected

from app import vdjbase_dbs, app, genomic_dbs
from db.vdjbase_api_query_filters import sample_info_filters
from db.vdjbase_model import HaplotypesFile, SamplesHaplotype, Allele, AllelesSample, Gene, AlleleConfidenceReport, HaplotypeEvidence
from db.vdjbase_airr_model import GenoDetection, SeqProtocol, Study, TissuePro, Patient, Sample

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
        if '_description' not in name:
            datasets.append({'dataset': name, 'description': vdjbase_dbs[species][name + '_description']})
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
                    if '_description' not in ds_name:
                        session = ds.session
                        novels = session.query(Allele).filter(Allele.novel == 1).all()

                        if novels:
                            if sp not in ret:
                                ret[sp] = {}
                            if ds_name not in ret[sp]:
                                ret[sp][ds_name] = {}
                            for novel in novels:
                                ret[sp][ds_name][novel.study_title] = (novel.seq.replace('.', ''), novel.appears)

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
            .filter(SeqProtocol.read_length == 'Full')\
            .all()

        results = {}
        for novel in novels:
            result = {
                'name': novel.study_title,
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
            best_hap = {'gene_type': '', 'hetero': False, 'count': 0, 'example': novel.samples[0].sample.study_title}

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

                novel_allele = novel.study_title.split('*')[1].upper()
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
                        best_hap = {'gene_type': gene_type, 'hetero': hetero, 'count': novel_count, 'example': haplo.sample.study_title}

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

        stats['description'] = vdjbase_dbs[species][dataset + '_description']
        stats['total_subjects'] = session.query(Patient.id).count()
        stats['total_samples'] = session.query(Sample.id).count()
        stats['sex_count'] = session.query(Patient.sex, func.count(Patient.sex)).group_by(Patient.sex).order_by(func.count(Patient.sex).desc()).all()
        stats['study_count'] = session.query(Study.study_title, func.count(Sample.name)).join(Sample, Sample.study_id == Study.id).group_by(Study.study_title).order_by(func.count(Sample.name).desc()).all()
        stats['condition_count'] = session.query(Patient.disease_diagnosis_label, func.count(Patient.disease_diagnosis_label)).group_by(Patient.disease_diagnosis_label).order_by(func.count(Patient.disease_diagnosis_label).desc()).all()
        stats['celltype_count'] = session.query(TissuePro.sub_cell_type, func.count(TissuePro.sub_cell_type)).group_by(TissuePro.sub_cell_type).order_by(func.count(TissuePro.sub_cell_type).desc()).all()
        stats['tissue_count'] = session.query(TissuePro.tissue_label, func.count(TissuePro.tissue_label)).group_by(TissuePro.tissue_label).order_by(func.count(TissuePro.tissue_label).desc()).all()
        studies = session.query(
            Study.study_title,
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
            row['subjects_in_vdjbase'] = session.query(Study.id).join(Patient, Patient.study_id == Study.id).filter(Study.study_title == row['name']).count()
            row['samples_in_vdjbase'] = session.query(Study.id).join(Sample, Sample.study_id == Study.id).filter(Study.study_title == row['name']).count()
            stats['studies'].append(row)

        return stats



rep_sample_bool_values = {
    'umi': ('UMI', '(blank)')
}

@ns.route('/sample_info/<string:species>/<string:dataset>/<string:sample>')
class SampleInfoApi(Resource):
    @digby_protected()
    def get(self, species, dataset, sample):
        """ Returns information on the selected sample """

        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        attribute_query = []

        for col in sample_info_filters.keys():
            if sample_info_filters[col]['field'] is not None:
                attribute_query.append(sample_info_filters[col]['field'])

        info = session.query(*attribute_query)\
            .join(GenoDetection, GenoDetection.id == Sample.geno_detection_id)\
            .join(Patient, Patient.id == Sample.patient_id)\
            .join(SeqProtocol)\
            .join(TissuePro)\
            .join(Study, Sample.study_id == Study.id)\
            .filter(Sample.name==sample).one_or_none()

        if info:
            info = info._asdict()
            info['date'] = info['date'].isoformat()
            haplotypes = session.query(HaplotypesFile.by_gene_s).join(SamplesHaplotype).join(Sample).filter(Sample.name==sample).order_by(HaplotypesFile.by_gene_s).all()
            info['haplotypes'] = [(h[0]) for h in haplotypes]

        return info



filter_arguments = reqparse.RequestParser()
filter_arguments.add_argument('page_number', type=int)
filter_arguments.add_argument('page_size', type=int)
filter_arguments.add_argument('filter', type=str)
filter_arguments.add_argument('sort_by', type=str)
filter_arguments.add_argument('cols', type=str)


@ns.route('/samples/<string:species>/<string:dataset>')
class SamplesApi(Resource):
    @digby_protected()
    @api.expect(filter_arguments, validate=True)
    def get(self, species, dataset):
        """ Returns the list of samples in the selected dataset """

        if species not in vdjbase_dbs or set(dataset.split(',')).difference(set(vdjbase_dbs[species])):
            return list()

        args = filter_arguments.parse_args(request)

        required_cols = json.loads(args['cols'])
        if 'name' not in required_cols:
            required_cols.append('name')
        if 'dataset' not in required_cols:
            required_cols.append('dataset')

        for col in required_cols:
            if col not in sample_info_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

        if 'genotypes' in required_cols:                # needed to compose paths to files
            for field in ('name', 'patient_name', 'study_name'):
                if field not in required_cols:
                    required_cols.append(field)
            required_cols.append('genotype')
            required_cols.append('genotype_stats')
            required_cols.append('genotype_report')

        attribute_query = [sample_info_filters['id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'id' and sample_info_filters[col]['field'] is not None:
                attribute_query.append(sample_info_filters[col]['field'])

        filter = json.loads(args['filter'])
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
                    if isinstance(el, datetime):
                        el = el.date().isoformat()
                    elif isinstance(el, bool):
                        if f in rep_sample_bool_values:
                            el = rep_sample_bool_values[f][0 if el else 1]
                    if (not isinstance(el, str) or len(el) > 0) and el not in uniques[f]:
                        uniques[f].append(el)
                    if isinstance(el, str) and len(el) == 0 and '(blank)' not in uniques[f]:
                        uniques[f].append('(blank)')

            if filter_applied:
                uniques['names_by_dataset'][s['dataset']].append(s['name'])

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
                if isinstance(v, datetime):
                    rec[k] = v.date().isoformat()
                elif isinstance(v, decimal.Decimal):
                    rec[k] = '%0.2f' % v

        if 'genotypes' in required_cols:
            for r in ret:
                r['genotypes'] = {}
                r['genotypes']['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['name'], 'sort_order': 'Locus'})

                r['genotypes']['path'] = app.config['BACKEND_LINK']
                sp = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset']]) + '/'
                r['genotypes']['tigger'] = sp + r['genotype'].replace('samples', '') if r['genotype'] else ''
                r['genotypes']['ogrdbstats'] = sp + r['genotype_stats'].replace('samples', '') if r['genotype_stats'] else ''
                r['genotypes']['ogrdbplot'] = sp + r['genotype_report'].replace('samples', '') if r['genotype_report'] else ''
                del r['genotype']
                del r['genotype_stats']
                del r['genotype_report']

                session = vdjbase_dbs[species][r['dataset']].session
                igsnper_path = session.query(Sample.igsnper_plot_path, Sample.genotype_report).filter(Sample.name == r['name']).one_or_none()

                if igsnper_path is not None and igsnper_path[0] is not None:
                    r['genotypes']['igsnper'] = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], igsnper_path[0]])
                else:
                    r['genotypes']['igsnper'] = ''


        if 'haplotypes' in required_cols:
            for r in ret:
                session = vdjbase_dbs[species][r['dataset']].session
                haplotypes = session.query(Sample.name, func.group_concat(HaplotypesFile.by_gene_s), func.group_concat(HaplotypesFile.file))
                h = haplotypes.filter(Sample.name == r['name'])\
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
                            r['haplotypes'][hap]['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['name'], 'hap_gene': hap, 'sort_order' : 'Locus'})
                            r['haplotypes'][hap]['rabhit'] = sp
                else:
                    r['haplotypes'] = ''

        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size'])
        }

def find_vdjbase_samples(attribute_query, species, datasets, filter):
    hap_filters = None
    allele_filters = None
    dataset_filters = []
    filter_spec = []

    for f in filter:
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
                    f['value'].append('')
                filter_spec.append(f)
        except:
            raise BadRequest('Bad filter string')
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
            hap_samples = session.query(Sample.name.distinct()).join(SamplesHaplotype).join(HaplotypesFile).filter(
                HaplotypesFile.by_gene_s.in_(hap_filters['value']))
            query = query.filter(Sample.name.in_(hap_samples))

        if allele_filters:
            allele_samples = session.query(Sample.name.distinct()).join(AllelesSample,
                                                                        Sample.id == AllelesSample.sample_id).join(
                Allele, Allele.id == AllelesSample.allele_id).filter(Allele.name.in_(allele_filters['value'])).all()
            if allele_samples is None:
                allele_samples = []
            query = query.filter(Sample.name.in_([s[0] for s in allele_samples]))

        res = query.all()

        for r in res:
            s = r._asdict()
            for k, v in s.items():
                if isinstance(v, datetime):
                    s[k] = v.date().isoformat()
            s['dataset'] = dset
            s['id'] = '%s.%d' % (dset, s['id'])
            ret.append(s)
    return ret


sequence_filters = {
    'name': {'model': 'Allele', 'field': Allele.name, 'sort': 'allele'},
    'gene_name': {'model': 'Gene', 'field': Gene.name.label('gene_name'), 'fieldname': 'name'},
    'igsnper_plot_path': {'model': 'Gene', 'field': Gene.igsnper_plot_path},
    'pipeline_name': {'model': 'Allele', 'field': Allele.pipeline_name},
    'seq': {'model': 'Allele', 'field': Allele.seq, 'no_uniques': True},
    'seq_len': {'model': 'Allele', 'field': Allele.seq_len, 'sort': 'numeric'},
    'similar': {'model': 'Allele', 'field': Allele.similar},
    'appears': {'model': 'Allele', 'field': Allele.appears, 'sort': 'numeric'},
    'is_single_allele': {'model': 'Allele', 'field': Allele.is_single_allele},
    'low_confidence': {'model': 'Allele', 'field': Allele.low_confidence},
    'novel': {'model': 'Allele', 'field': Allele.novel},
    'max_kdiff': {'model': 'Allele', 'field': Allele.max_kdiff, 'sort': 'numeric'},
    'notes': {'model': Allele, 'field': func.group_concat(AlleleConfidenceReport.notes, '\n').label('notes')},
    'notes_count': {'model': Allele, 'field': func.count(AlleleConfidenceReport.id).label('notes_count'), 'sort': 'numeric'},

    'sample_id': {'model': None, 'fieldname': 'sample_id'},
    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset', 'no_uniques': True},
}

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
            return list()

        required_cols = json.loads(args['cols'])

        for col in required_cols:
            if col not in sequence_filters.keys():
                print('bad column in request: %s' % col)
                return list(), 404

        if 'gene_name' in required_cols:
            required_cols.append('igsnper_plot_path')

        datasets = dataset.split(',')
        ret = find_vdjbase_sequences(species, datasets, required_cols, json.loads(args['filter']))
        total_size = len(ret)

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in ret:
            for f in required_cols:
                if sequence_filters[f]['field'] is not None and 'no_uniques' not in sequence_filters[f]:
                    el = s[f]
                    if isinstance(el, datetime):        # convert and add to uniques. Can't convert the returned records until we sort
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
                    uniques[f] = [x[1] for x in uniques[f].sorted(key=lambda x: (x is None or x == '', x))]
            except:
                pass

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        if len(sort_specs) > 0:
            for spec in sort_specs:
                if spec['field'] in sequence_filters.keys():
                    if spec['field'] in ('name'):
                        ret = sorted(ret, key=lambda x : allele_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                    else:
                        ret = sorted(ret, key=lambda x : ((x[spec['field']] is None or x[spec['field']] == ''),  x[spec['field']]), reverse=(spec['order'] == 'desc'))
        else:
            ret = sorted(ret, key=lambda x : allele_sort_key(x['name']))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first : first + args['page_size']]

        for rec in ret:
            for k, v in rec.items():
                if isinstance(v, datetime):
                    rec[k] = v.date().isoformat()
                elif isinstance(v, decimal.Decimal):
                    rec[k] = '%0.2f' % v



        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size'])
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
                            value.append('|' + v.replace(',', '|') + '|')
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
                        f['value'].append('')
                    filter_spec.append(f)
            except:
                raise BadRequest('Bad filter string: %s' % f)

    if 'notes_count' in required_cols and 'notes' not in required_cols:
        required_cols.append('notes')
    ret = []
    if len(dataset_filters) > 0:
        apply_filter_to_list(datasets, dataset_filters)
    for dset in datasets:
        session = vdjbase_dbs[species][dset].session

        attribute_query = []

        for col in required_cols:
            if sequence_filters[col]['field'] is not None:
                attribute_query.append(sequence_filters[col]['field'])

        query = session.query(*attribute_query).join(Gene, Allele.gene_id == Gene.id)
        query = apply_filters(query, filter_spec)

        if 'notes' in required_cols:
            query = query.outerjoin(AlleleConfidenceReport, Allele.id == AlleleConfidenceReport.allele_id).group_by(
                Allele.id)

        res = query.all()

        required_names = []

        if sample_id_filter is not None:
            required_ids = []
            if dset in sample_id_filter['value']:
                for id in sample_id_filter['value'][dset]:
                    required_ids.append(id)
            alleles_with_samples = session.query(Allele.name) \
                .join(AllelesSample) \
                .join(Sample) \
                .filter(Sample.name.in_(required_ids)).all()

            for a in alleles_with_samples:
                required_names.append(a[0])
            # query = query.filter(Allele.name.in_(required_names))

        for r in res:
            if len(required_names) == 0 or r.study_title in required_names:
                s = r._asdict()
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


# Apply filter params to a list of samples in the context of a specific dataset
# wanted_genes is returned in the required search order

def apply_rep_filter_params(params, sample_list, session):
    if 'per_sample' in params:
        sample_list = filter_per_sample(params['per_sample'], sample_list)
    if 'sort_order' in params:
        if params['sort_order'] == 'Alphabetic':
            gq = session.query(Gene).order_by(Gene.alpha_order)
        else:
            gq = session.query(Gene).order_by(Gene.locus_order)
    else:
        gq = session.query(Gene).order_by(Gene.alpha_order)

    if len(params['f_gene_types']) > 0:
        gq = gq.filter(Gene.type.in_(params['f_gene_types']))
    if len(params['f_genes']) > 0:
        gq = gq.filter(Gene.name.in_(params['f_genes']))
    if not params['f_pseudo_genes']:
        gq = gq.filter(Gene.pseudo_gene == 0)
    wanted_genes = gq.all()
    wanted_genes = [gene.study_title for gene in wanted_genes]
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


def get_order_file(species, dataset, locus_order=True):
    file_name = os.path.join(app.config['OUTPUT_PATH'], '%s_%s_%s_order.tsv' % (species, dataset, 'locus' if locus_order else 'alpha'))

    #if not os.path.isfile(file_name):
    session = vdjbase_dbs[species][dataset].session
    if locus_order:
        gene_order = session.query(Gene.name).order_by(Gene.locus_order).all()
    else:
        gene_order = session.query(Gene.name).order_by(Gene.alpha_order).all()
    gene_order = [x[0] for x in gene_order]

    with open(file_name, 'w') as fo:
        fo.write('\n'.join(gene_order))

    return file_name


# interleave content from multiple files, or return one of them if they are all the same
def get_multiple_order_file(species, datasets, locus_order=True):
    gene_order = []
    added = False

    for dataset in datasets:
        with open(get_order_file(species, dataset, locus_order), 'r') as fi:
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

    if added:
        file_name = make_output_file('tsv')

        with open(file_name, 'w') as fo:
            fo.write('\n'.join(gene_order))
    else:
        print(list(datasets))
        file_name = get_order_file(species, list(datasets)[0], locus_order)

    return file_name






