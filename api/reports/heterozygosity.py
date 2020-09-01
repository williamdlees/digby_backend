# Heterozygosity plot for AIRR-seq samples

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, HaplotypesFile, SamplesHaplotype, AllelesSample, Gene, Allele, Patient, AllelesPattern
import os
from api.vdjbase.vdjbase import VDJBASE_SAMPLE_PATH, apply_rep_filter_params
from sqlalchemy import func
import pandas as pd

HETEROZYGOSITY_SCRIPT = 'Heterozygous.R'

def run(format, species, genomic_samples, rep_samples, params):
    if len(rep_samples) == 0:
        raise BadRequest('No repertoire-derived genotypes were selected.')

    if format != 'html':
        raise BadRequest('Invalid format requested')

    kdiff = params['f_kdiff'] if 'f_kdiff' in params else 0

    samples_by_dataset = {}
    for rep_sample in rep_samples:
        if rep_sample['dataset'] not in samples_by_dataset:
            samples_by_dataset[rep_sample['dataset']] = []
        samples_by_dataset[rep_sample['dataset']].append(rep_sample['name'])

    # Format we need to produce is [(gene_name, hetero count, homo count),...]

    gene_hetrozygous_dis = {}

    for dataset in samples_by_dataset.keys():
        session = vdjbase_dbs[species][dataset].session
        sample_list = session.query(Sample.name, Sample.genotype, Sample.patient_id).filter(Sample.name.in_(samples_by_dataset[dataset])).all()
        sample_list, wanted_genes = apply_rep_filter_params(params, sample_list, session)
        sample_list = [s[0] for s in sample_list]

        query = session.query(Gene.name, Patient.id, Allele.id, Sample.name, Gene.locus_order, AllelesSample.kdiff, Allele.name) \
            .join(Allele) \
            .join(AllelesSample) \
            .join(Sample) \
            .join(Patient) \
            .filter(Gene.name.in_(wanted_genes)) \
            .filter(Allele.name.notlike('%Del%')) \
            .filter(Allele.name.notlike('%OR%')) \
            .filter(Sample.name.in_(sample_list)) \
            .filter(AllelesSample.kdiff >= kdiff) \
            .order_by(Gene.locus_order, Patient.id, Allele.id)

        if(params['ambiguous_alleles'] == 'Exclude'):
            query = query.filter(Allele.is_single_allele == True)

        allele_sample_recs = query.all()

        # As the result is indexed, run over each gene in turn, count the number of alleles found in each patient, update h_counts accordingly

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

                # If we have both an unambiguous allele and an ambiguous allele containing that unambiguous one,
                # drop the unambiguous one because it is already counted

                if(params['ambiguous_alleles'] != 'Exclude'):
                    patterns = session.query(AllelesPattern.pattern_id)\
                        .filter(AllelesPattern.allele_in_p_id.in_(patient_allele_ids))\
                        .filter(AllelesPattern.pattern_id.in_(patient_allele_ids))\
                        .all()

                    if patterns is not None and len(patterns) >0:
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
    attachment_filename = '%s_heterozygosity_plot.pdf' % species

    cmd_line = ["-i", haplo_path,
                "-o", output_path,
                "-s", SYSDATA]

    if run_rscript(HETEROZYGOSITY_SCRIPT, cmd_line) and os.path.isfile(output_path) and os.path.getsize(output_path) != 0:
        return send_report(output_path, format, attachment_filename)
    else:
        raise BadRequest('No output from report')
