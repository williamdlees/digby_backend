# Services related to vdjbase repseq-based data sets

from flask import request
from flask_restx import Resource, reqparse, fields, marshal, inputs
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

from app import vdjbase_dbs, app
from db.vdjbase_model import Sample, GenoDetection, Patient, SeqProtocol, Study, TissuePro, HaplotypesFile, SamplesHaplotype, Allele, AllelesSample, Gene, AlleleConfidenceReport

VDJBASE_SAMPLE_PATH = os.path.join(app.config['STATIC_PATH'], 'study_data/VDJbase/samples')



# Return SqlAlchemy row as a dict, using correct column names
def object_as_dict(obj):
    return {c.key: getattr(obj, c.key)
            for c in inspect(obj).mapper.column_attrs}


ns = api.namespace('repseq', description='Genes and annotations inferred from RepSeq data')

@ns.route('/species')
@api.response(404, 'No species available!')
class SpeciesApi(Resource):

    def get(self):

        """ Returns the list of species for which information is held """
        sp = list(vdjbase_dbs.keys())

        if sp:
            return sp
        else:
            return None, 404

def find_datasets(species):
    datasets = []
    for name in vdjbase_dbs[species].keys():
        if '_description' not in name:
            datasets.append({'dataset': name, 'description': vdjbase_dbs[species][name + '_description']})
    return datasets


@ns.route('/ref_seqs/<string:species>')
class DataSetAPI(Resource):
    def get(self, species):
        """ Returns the list of datasets available for the selected species """
        if species in vdjbase_dbs:
            return find_datasets(species)
        else:
            return list()

@ns.route('/dataset_info/<string:species>/<string:dataset>')
class DataSetInfoAPI(Resource):

    def get(self, species, dataset):
        """Returns information and statistics on the dataset"""
        if species not in vdjbase_dbs and dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        stats = {}

        stats['description'] = vdjbase_dbs[species][dataset + '_description']
        stats['total_subjects'] = session.query(Patient.id).count()
        stats['total_samples'] = session.query(Sample.id).count()
        stats['sex_count'] = session.query(Patient.sex, func.count(Patient.sex)).group_by(Patient.sex).order_by(func.count(Patient.sex).desc()).all()
        stats['study_count'] = session.query(Study.name, func.count(Sample.name)).join(Sample, Sample.study_id == Study.id).group_by(Study.name).order_by(func.count(Sample.name).desc()).all()
        stats['condition_count'] = session.query(Patient.status, func.count(Patient.status)).group_by(Patient.status).order_by(func.count(Patient.status).desc()).all()
        stats['celltype_count'] = session.query(TissuePro.sub_cell_type, func.count(TissuePro.sub_cell_type)).group_by(TissuePro.sub_cell_type).order_by(func.count(TissuePro.sub_cell_type).desc()).all()
        stats['tissue_count'] = session.query(TissuePro.tissue, func.count(TissuePro.tissue)).group_by(TissuePro.tissue).order_by(func.count(TissuePro.tissue).desc()).all()
        studies = session.query(
            Study.name,
            Study.accession_id,
            Study.institute,
            Study.num_subjects,
            Study.num_samples,
            Study.reference,
            Study.accession_reference
        ).all()
        stats['studies'] = [row._asdict() for row in studies]

        return stats


valid_filters = {
    'name': {'model': 'Sample', 'field': Sample.name},
    'id': {'model': 'Sample', 'field': Sample.id},
    'chain': {'model': 'Sample', 'field': Sample.chain},
    'row_reads': {'model': 'Sample', 'field': Sample.row_reads},
    'date': {'model': 'Sample', 'field': Sample.date},
    'samples_group': {'model': 'Sample', 'field': Sample.samples_group},

    'patient_name': {'model': 'Patient', 'field': Patient.name.label('patient_name'), 'fieldname': 'name'},
    'sex': {'model': 'Patient', 'field': Patient.sex},
    'ethnic': {'model': 'Patient', 'field': Patient.ethnic},
    'status': {'model': 'Patient', 'field': Patient.status},
    'cohort': {'model': 'Patient', 'field': Patient.cohort},
    'age': {'model': 'Patient', 'field': Patient.age},
    'name_in_paper': {'model': 'Patient', 'field': Patient.name_in_paper},
    'country': {'model': 'Patient', 'field': Patient.country},

    'study_name': {'model': 'Study', 'field': Study.name.label('study_name'), 'fieldname': 'name'},
    'institute': {'model': 'Study', 'field': Study.institute},
    'researcher': {'model': 'Study', 'field': Study.researcher},
    'num_subjects': {'model': 'Study', 'field': Study.num_subjects},
    'num_samples': {'model': 'Study', 'field': Study.num_samples},
    'study_reference': {'model': 'Study', 'field': Study.reference.label('study_reference'), 'fieldname': 'reference'},
    'study_contact': {'model': 'Study', 'field': Study.contact.label('study_contact'), 'fieldname': 'contact'},
    'study_acc_id': {'model': 'Study', 'field': Study.accession_id.label('study_acc_id'), 'fieldname': 'accession_id'},
    'study_acc_ref': {'model': 'Study', 'field': Study.accession_reference.label('study_acc_ref'), 'fieldname': 'accession_reference'},

    'tissue_name': {'model': 'TissuePro', 'field': TissuePro.name.label('tissue_name'), 'fieldname': 'name'},
    'species': {'model': 'TissuePro', 'field': TissuePro.species},
    'tissue': {'model': 'TissuePro', 'field': TissuePro.tissue},
    'combined_cell_type': {'model': 'TissuePro', 'field': TissuePro.combined_cell_type},
    'cell_type': {'model': 'TissuePro', 'field': TissuePro.cell_type},
    'sub_cell_type': {'model': 'TissuePro', 'field': TissuePro.sub_cell_type},
    'isotype': {'model': 'TissuePro', 'field': TissuePro.isotype},

    'sequencing_protocol': {'model': 'SeqProtocol', 'field': SeqProtocol.name.label('sequencing_protocol'), 'fieldname': 'name'},
    'umi': {'model': 'SeqProtocol', 'field': SeqProtocol.umi},
    'sequencing_length': {'model': 'SeqProtocol', 'field': SeqProtocol.sequencing_length},
    'primer_3_loc': {'model': 'SeqProtocol', 'field': SeqProtocol.primers_3_location.label('primer_3_loc'), 'fieldname': 'primers_3_location'},
    'primer_5_loc': {'model': 'SeqProtocol', 'field': SeqProtocol.primers_5_location.label('primer_5_loc'), 'fieldname': 'primers_5_location'},
    'seq_platform': {'model': 'SeqProtocol', 'field': SeqProtocol.sequencing_platform.label('seq_platform'), 'fieldname': 'sequencing_platform'},
    'helix': {'model': 'SeqProtocol', 'field': SeqProtocol.helix},

    'pipeline_name': {'model': 'GenoDetection', 'field': GenoDetection.name.label('pipeline_name'), 'fieldname': 'name'},
    'prepro_tool': {'model': 'GenoDetection', 'field': GenoDetection.prepro_tool},
    'aligner_tool': {'model': 'GenoDetection', 'field': GenoDetection.aligner_tool},
    'aligner_ver': {'model': 'GenoDetection', 'field': GenoDetection.aligner_ver},
    'aligner_reference': {'model': 'GenoDetection', 'field': GenoDetection.aligner_reference},
    'geno_tool': {'model': 'GenoDetection', 'field': GenoDetection.geno_tool},
    'geno_ver': {'model': 'GenoDetection', 'field': GenoDetection.geno_ver},
    'haplotype_tool': {'model': 'GenoDetection', 'field': GenoDetection.haplotype_tool},
    'haplotype_ver': {'model': 'GenoDetection', 'field': GenoDetection.haplotype_ver},
    'single_assignment': {'model': 'GenoDetection', 'field': GenoDetection.single_assignment},
    'detection': {'model': 'GenoDetection', 'field': GenoDetection.detection},

    'allele': {'model': None, 'field': None},

    'haplotypes': {'model': None, 'field': None},
    'genotypes': {'model': None, 'field': None},

    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset', 'no_uniques': True},
}



@ns.route('/sample_info/<string:species>/<string:dataset>/<string:sample>')
class SampleInfoApi(Resource):
    def get(self, species, dataset, sample):
        """ Returns information on the selected sample """

        if species not in vdjbase_dbs or dataset not in vdjbase_dbs[species]:
            return None, 404

        session = vdjbase_dbs[species][dataset].session
        attribute_query = []

        for col in valid_filters.keys():
            if valid_filters[col]['field'] is not None:
                attribute_query.append(valid_filters[col]['field'])

        info = session.query(*attribute_query).join(GenoDetection).join(Patient).join(SeqProtocol).join(TissuePro).join(Study, Sample.study_id == Study.id).filter(Sample.name==sample).one_or_none()

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
            if col not in valid_filters.keys():
                raise BadRequest('Bad filter string %s' % args['filter'])

        if 'genotypes' in required_cols:                # needed to compose paths to files
            for field in ('name', 'patient_name', 'study_name'):
                if field not in required_cols:
                    required_cols.append(field)

        attribute_query = [valid_filters['id']['field']]        # the query requires the first field to be from Sample

        for col in required_cols:
            if col != 'id' and valid_filters[col]['field'] is not None:
                attribute_query.append(valid_filters[col]['field'])

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
                if valid_filters[f]['field'] is not None and 'no_uniques' not in valid_filters[f]:
                    el = s[f]
                    if isinstance(el, datetime):
                        s[f] = el.date().isoformat()
                    if (not isinstance(el, str) or len(s[f]) > 0) and s[f] not in uniques[f]:
                        uniques[f].append(s[f])
                    if isinstance(el, str) and len(el) == 0 and '(blank)' not in uniques[f]:
                        uniques[f].append('(blank)')

            if filter_applied:
                uniques['names_by_dataset'][s['dataset']].append(s['name'])

        def name_sort_key(name):
            name = name.split('_')
            for i in range(len(name)):
                name[i] = name[i][1:].zfill(4)
            return name

        for f in required_cols:
            try:
                if f in ('name', 'patient_name'):
                    uniques[f].sort(key=lambda x: name_sort_key(x))
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

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        if len(sort_specs) > 0:
            for spec in sort_specs:
                if spec['field'] in valid_filters.keys():
                    if spec['field'] in ('name', 'patient_name'):
                        ret = sorted(ret, key=lambda x: name_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                    else:
                        ret = sorted(ret, key=lambda x: ((x[spec['field']] is None or x[spec['field']] == ''),  x[spec['field']]), reverse=(spec['order'] == 'desc'))
        else:
            ret = sorted(ret, key=lambda x: name_sort_key(x['name']))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first:first + args['page_size']]

        if 'genotypes' in required_cols:
            for r in ret:
                # TODO: this section should be changed to use database fields to find the paths.
                r['genotypes'] = {}
                r['genotypes']['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['name']})

                r['genotypes']['path'] = app.config['BACKEND_LINK']
                dp = os.path.join(species, r['dataset'], r['study_name'], r['patient_name'])
                if os.path.isdir(os.path.join(VDJBASE_SAMPLE_PATH, dp)):            # old format
                    fp = os.path.join(VDJBASE_SAMPLE_PATH, dp, r['name'] + '_')
                    sp = '/'.join(['static/study_data/VDJbase/samples', dp, r['name'] + '_'])
                    r['genotypes']['tigger'] = sp + 'geno_H_binom.tab' if isfile(fp + 'geno_H_binom.tab') else ''
                    r['genotypes']['ogrdbstats'] = sp + 'genotyped_mut_ogrdb_report.csv' if isfile(fp + 'genotyped_mut_ogrdb_report.csv') else ''
                    r['genotypes']['ogrdbplot'] = sp + '_ogrdb_plots.pdf' if isfile(fp + '_ogrdb_plots.pdf') else ''
                else:                                                               # new format
                    dp = os.path.join(species, r['dataset'], r['study_name'], r['name']).replace('\\', '/')
                    fp = os.path.join(VDJBASE_SAMPLE_PATH, dp, r['name']).replace('\\', '/')
                    sp = '/'.join(['static/study_data/VDJbase/samples', dp, r['name']])
                    r['genotypes']['tigger'] = sp + '.tsv' if isfile(fp + '.tsv') else ''
                    r['genotypes']['ogrdbstats'] = sp + '_ogrdb_report.csv' if isfile(fp + '_ogrdb_report.csv') else ''
                    r['genotypes']['ogrdbplot'] = sp + '_ogrdb_plots.pdf' if isfile(fp + '_ogrdb_plots.pdf') else ''

                session = vdjbase_dbs[species][r['dataset']].session
                igsnper_path = session.query(Sample.igsnper_plot_path).filter(Sample.name == r['name']).one_or_none()

                if igsnper_path is not None:
                    r['genotypes']['igsnper'] = '/'.join(['static/study_data/VDJbase/samples', species, r['dataset'], igsnper_path[0]])
                else:
                    r['genotypes']['igsnper'] = ''



        if 'haplotypes' in required_cols:
            for r in ret:
                session = vdjbase_dbs[species][r['dataset']].session
                haplotypes = session.query(Sample.name, func.group_concat(HaplotypesFile.by_gene_s), func.group_concat(HaplotypesFile.file))
                h = haplotypes.filter(Sample.name == r['name'])\
                    .join(SamplesHaplotype, SamplesHaplotype.samples_id == Sample.id)\
                    .join(HaplotypesFile, HaplotypesFile.id == SamplesHaplotype.haplotypes_files_id)\
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
                            r['haplotypes'][hap]['analysis'] = json.dumps({'species': species, 'repSeqs': [r['dataset']], 'name': r['name'], 'hap_gene': hap})
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
                f['model'] = valid_filters[f['field']]['model']
                if 'fieldname' in valid_filters[f['field']]:
                    f['field'] = valid_filters[f['field']]['fieldname']
                if '(blank)' in f['value']:
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


valid_sequence_cols = {
    'name': {'model': 'Allele', 'field': Allele.name},
    'gene_name': {'model': 'Gene', 'field': Gene.name.label('gene_name'), 'fieldname': 'name'},
    'igsnper_plot_path': {'model': 'Gene', 'field': Gene.igsnper_plot_path},
    'pipeline_name': {'model': 'Allele', 'field': Allele.pipeline_name},
    'seq': {'model': 'Allele', 'field': Allele.seq, 'no_uniques': True},
    'seq_len': {'model': 'Allele', 'field': Allele.seq_len},
    'similar': {'model': 'Allele', 'field': Allele.similar},
    'appears': {'model': 'Allele', 'field': Allele.appears},
    'is_single_allele': {'model': 'Allele', 'field': Allele.is_single_allele},
    'low_confidence': {'model': 'Allele', 'field': Allele.low_confidence},
    'novel': {'model': 'Allele', 'field': Allele.novel},
    'max_kdiff': {'model': 'Allele', 'field': Allele.max_kdiff},
    'notes': {'model': Allele, 'field': func.group_concat(AlleleConfidenceReport.notes, '\n').label('notes')},
    'notes_count': {'model': Allele, 'field': func.count(AlleleConfidenceReport.id).label('notes_count')},

    'sample_id': {'model': None, 'fieldname': 'sample_id'},
    'dataset': {'model': None, 'field': None, 'fieldname': 'dataset', 'no_uniques': True},
}

@ns.route('/sequences/<string:species>/<string:dataset>')
class SequencesApi(Resource):
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
            if col not in valid_sequence_cols.keys():
                print('bad column in request: %s' % col)
                return list(), 404

        if 'gene_name' in required_cols:
            required_cols.append('igsnper_plot_path')

        sample_id_filter = None
        filter_spec = []
        dataset_filters = []

        if args['filter']:
            for f in json.loads(args['filter']):
               # try:
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
                    f['model'] = valid_sequence_cols[f['field']]['model']
                    if 'fieldname' in valid_sequence_cols[f['field']]:
                        f['field'] = valid_sequence_cols[f['field']]['fieldname']
                    if '(blank)' in f['value']:
                        f['value'].append('')
                    filter_spec.append(f)
                # except:
                  #  raise BadRequest('Bad filter string %s' % args['filter'])

        if 'notes_count' in required_cols and 'notes' not in required_cols:
            required_cols.append('notes')

        ret = []

        datasets = dataset.split(',')

        if len(dataset_filters) > 0:
            apply_filter_to_list(datasets, dataset_filters)

        for dset in datasets:
            session = vdjbase_dbs[species][dset].session

            attribute_query = []

            for col in required_cols:
                if valid_sequence_cols[col]['field'] is not None:
                    attribute_query.append(valid_sequence_cols[col]['field'])

            query = session.query(*attribute_query).join(Gene, Allele.gene_id == Gene.id)
            query = apply_filters(query, filter_spec)

            if sample_id_filter is not None:
                required_ids = []
                if dset in sample_id_filter['value']:
                    for id in sample_id_filter['value'][dset]:
                        required_ids.append(id)
                alleles_with_samples = session.query(Allele.name)\
                    .join(AllelesSample)\
                    .join(Sample)\
                    .filter(Sample.name.in_(required_ids)).all()

                required_names = []
                for a in alleles_with_samples:
                    required_names.append(a[0])
                query = query.filter(Allele.name.in_(required_names))

            if 'notes' in required_cols:
                query = query.outerjoin(AlleleConfidenceReport, Allele.id == AlleleConfidenceReport.allele_id).group_by(Allele.id)

            res = query.all()

            for r in res:
                s = r._asdict()
                for k, v in s.items():
                    if isinstance(v, datetime):
                        s[k] = v.date().isoformat()
                    elif isinstance(v, decimal.Decimal):
                        s[k] = '%0.2f' % v

                    if k == 'similar' and v is not None:
                        s[k] = v.replace('|', '')
                s['dataset'] = dset

                if len(s['igsnper_plot_path']) > 0:
                    s['igsnper_plot_path'] = '/'.join([app.config['BACKEND_LINK'], 'static/study_data/VDJbase/samples', species, dset, s['igsnper_plot_path']])

                ret.append(s)


        total_size = len(ret)

        uniques = {}
        for f in required_cols:
            uniques[f] = []

        for s in ret:
            for f in required_cols:
                if valid_sequence_cols[f]['field'] is not None and 'no_uniques' not in valid_sequence_cols[f]:
                    el = s[f]
                    if isinstance(el, datetime):
                        s[f] = el.date().isoformat()
                    elif isinstance(el, decimal.Decimal):
                        s[f] = '%0.2f' % el
                    if (not isinstance(el, str) or len(s[f]) > 0) and s[f] not in uniques[f]:
                        uniques[f].append(s[f])
                    if isinstance(el, str) and len(el) == 0 and '(blank)' not in uniques[f]:
                        uniques[f].append('(blank)')

        uniques['dataset'] = dataset.split(',')

        # For gene order, it would be a good idea if datasets from the same species agreed!

        gene_order = session.query(Gene.name, Gene.alpha_order).all()
        gene_order = {x[0]: x[1] for x in gene_order}

        def allele_sort_key(name):
            if '*' in name:
                gene = name.split('*')
            else:
                gene = (name, '')

            return((gene_order[gene[0]] if gene[0] in gene_order else 999, gene[1]))

        for f in required_cols:
            try:
                if f in ('name'):
                    uniques[f].sort(key=lambda x: allele_sort_key(x))
                else:
                    uniques[f] = [x[1] for x in uniques[f].sorted(key=lambda x: (x is None or x == '', x))]
            except:
                pass

        sort_specs = json.loads(args['sort_by']) if ('sort_by' in args and args['sort_by'] != None)  else [{'field': 'name', 'order': 'asc'}]

        if len(sort_specs) > 0:
            for spec in sort_specs:
                if spec['field'] in valid_sequence_cols.keys():
                    if spec['field'] in ('name'):
                        ret = sorted(ret, key=lambda x : allele_sort_key(x[spec['field']]), reverse=(spec['order'] == 'desc'))
                    else:
                        ret = sorted(ret, key=lambda x : ((x[spec['field']] is None or x[spec['field']] == ''),  x[spec['field']]), reverse=(spec['order'] == 'desc'))
        else:
            ret = sorted(ret, key=lambda x : allele_sort_key(x['name']))

        if args['page_size']:
            first = (args['page_number']) * args['page_size']
            ret = ret[first : first + args['page_size']]

        return {
            'samples': ret,
            'uniques': uniques,
            'total_items': total_size,
            'page_size': args['page_size'],
            'pages': ceil((total_size*1.0)/args['page_size'])
        }


def find_rep_filter_params(species, datasets):
    if species not in vdjbase_dbs or set(datasets).difference(set(vdjbase_dbs[species])):
        raise BadRequest('Unknown AIRR-seq dataset')

    genes = []
    gene_types = []

    for dataset in datasets:
        session = vdjbase_dbs[species][dataset].session
        g_q = session.query(Gene.name, Gene.type).all()
        genes.extend((res[0] for res in g_q))
        gene_types.extend((res[1] for res in g_q))

    genes = sorted(set(genes))
    gene_types = sorted(set(gene_types))

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

    return params


PSEUDO_GENES = [
    "IGHV2-10", "IGHV3-52", "IGHV3-47", "IGHV3-71", "IGHV3-22", "IGHV4-55", "IGHV1-68",
                    "IGHV5-78", "IGHV3-32", "IGHV3-33-2", "IGHV3-38-3", "IGHV3-25", "IGHV3-19", "IGHV7-40", "IGHV3-63",
                    "IGHV3-62", "IGHV3-29", "IGHV3-54", "IGHV1-38-4", "IGHV7-34-1", "IGHV1-38-4", "IGHV3-30-2",
                    "IGHV3-69-1", "IGHV3-30-22", "IGHV1-f", "IGHV3-30-33", "IGHV3-38", "IGHV7-81", "IGHV3-35",
                    "IGHV3-16","IGHV3-30-52","IGHV1-69D", "IGHD1-14", "IGHV3-30-42"
]

# Apply filter params to a list of samples in the context of a specific dataset

def apply_rep_filter_params(params, sample_list, session):
    if 'per_sample' in params:
        sample_list = filter_per_sample(params['per_sample'], sample_list)
    gq = session.query(Gene)
    if len(params['f_gene_types']) > 0:
        gq = gq.filter(Gene.type.in_(params['f_gene_types']))
    if len(params['f_genes']) > 0:
        gq = gq.filter(Gene.name.in_(params['f_genes']))
    if not params['f_pseudo_genes']:
        gq = gq.filter(Gene.name.notin_(PSEUDO_GENES))
    wanted_genes = gq.all()
    wanted_genes = [gene.name for gene in wanted_genes]
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



