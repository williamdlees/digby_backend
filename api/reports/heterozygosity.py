from werkzeug.exceptions import BadRequest
from api.reports.reports import send_report
from api.reports.report_utils import make_output_file, collate_samples, chunk_list
from app import vdjbase_dbs
from db.vdjbase_model import AllelesSample, Gene, Allele, AllelesPattern
from db.vdjbase_airr_model import Patient, Sample
import os
from api.vdjbase.vdjbase import apply_rep_filter_params
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from api.reports.Python_scripts.Heterozygous import create_heterozygous_plot

SAMPLE_CHUNKS = 400

def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format != 'html':
        raise BadRequest('Invalid format requested')

    kdiff = float(params['f_kdiff']) if 'f_kdiff' in params and params['f_kdiff'] != '' else 0
    chain, samples_by_dataset = collate_samples(rep_samples)

    gene_hetrozygous_dis = {}

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        allele_sample_recs = []

        for sample_chunk in chunk_list(samples_by_dataset[dataset], SAMPLE_CHUNKS):
            sample_list = session.query(Sample.sample_name, Sample.genotype, Sample.patient_id).filter(Sample.sample_name.in_(sample_chunk)).all()
            sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
            sample_list = [s[0] for s in sample_list]

            query = session.query(Gene.name, Patient.id, Allele.id, Sample.sample_name, Gene.locus_order, AllelesSample.kdiff, Allele.name) \
                .join(Allele, Gene.id == Allele.gene_id) \
                .join(AllelesSample, Allele.id == AllelesSample.allele_id) \
                .join(Sample, Sample.id == AllelesSample.sample_id) \
                .join(Patient, Patient.id == Sample.patient_id) \
                .filter(Gene.name.in_(wanted_genes)) \
                .filter(Allele.name.notlike('%Del%')) \
                .filter(Allele.name.notlike('%OR%')) \
                .filter(Sample.sample_name.in_(sample_list)) \
                .filter(AllelesSample.kdiff >= kdiff)

            if 'sort_order' in params and params['sort_order'] == 'Locus':
                query = query.order_by(Gene.locus_order, Patient.id, Allele.id)
            else:
                query = query.order_by(Gene.alpha_order, Patient.id, Allele.id)

            if params['ambiguous_alleles'] == 'Exclude':
                query = query.filter(Allele.is_single_allele == True)

            allele_sample_recs.extend(query.all())

        i = 0
        target_gene = ''

        while i < len(allele_sample_recs):
            target_gene = allele_sample_recs[i][0]
            h_counts = [0, 0]

            while i < len(allele_sample_recs):
                if allele_sample_recs[i][0] != target_gene:
                    break

                target_patient = allele_sample_recs[i][1]
                patient_allele_ids = []

                while i < len(allele_sample_recs):
                    if allele_sample_recs[i][0] != target_gene or allele_sample_recs[i][1] != target_patient:
                        break

                    patient_allele_ids.append(allele_sample_recs[i][2])
                    i += 1

                patient_allele_ids = set(patient_allele_ids)

                if params['ambiguous_alleles'] != 'Exclude':
                    patterns = session.query(AllelesPattern.pattern_id)\
                        .filter(AllelesPattern.allele_in_p_id.in_(patient_allele_ids))\
                        .filter(AllelesPattern.pattern_id.in_(patient_allele_ids))\
                        .all()

                    if patterns is not None and len(patterns) > 0:
                        patterns = set([pattern[0] for pattern in patterns])
                        patient_allele_ids = patient_allele_ids - patterns

                if len(patient_allele_ids) > 1:
                    h_counts[1] += 1
                elif len(patient_allele_ids) > 0:
                    h_counts[0] += 1

            if target_gene not in gene_hetrozygous_dis:
                gene_hetrozygous_dis[target_gene] = (target_gene, h_counts[0], h_counts[1])
            else:
                gene_hetrozygous_dis[target_gene] = (target_gene, gene_hetrozygous_dis[target_gene][1] + h_counts[0], gene_hetrozygous_dis[target_gene][2] + h_counts[1])

    haplo_path = make_output_file('tab')
    labels = ['GENE', 'HM', 'HT']
    df = pd.DataFrame(gene_hetrozygous_dis.values(), columns=labels)
    df.to_csv(haplo_path, sep='\t', index=False)
    output_path = make_output_file('html')

    #Plot Creation
    create_heterozygous_plot(haplo_path, output_path, chain)

    if os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format)
    else:
        raise BadRequest('No output from report')


