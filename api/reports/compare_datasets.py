# Crude comparison of two repseq datasets

from werkzeug.exceptions import BadRequest
from api.reports.reports import SYSDATA, run_rscript, send_report, make_output_file
from app import app, vdjbase_dbs
from db.vdjbase_model import Sample, Allele, AllelesSample, Gene, GenesDistribution, AllelesPattern, \
    AlleleConfidenceReport

HEATMAP_HAPLOTYPE_SCRIPT = "haplotype_heatmap.R"


def run(format, species, genomic_datasets, genomic_samples, rep_datasets, rep_samples, params):
    datasets = []
    for sample in rep_samples:
        if sample['dataset'] not in datasets:
            datasets.append(sample['dataset'])

    if len(datasets) != 2:
        raise BadRequest('Please select exactly two AIRR-seq datasets to compare.')

    if format != 'html':
        raise BadRequest('Invalid format requested')

    output_path = make_output_file('html')

    with open(output_path, 'w') as fo:
        session = []
        alleles = []
        allele_names = []
        allele_similars = []
        allele_lookups = [{}, {}]

        for i in (0, 1):
            session.append(vdjbase_dbs[species][datasets[i]].session)
            alleles.append(session[i].query(Allele).all())
            allele_names.append(set([allele.name for allele in alleles[i]]))

            allele_similars.append({})
            for allele in alleles[i]:
                if allele.similar is not None:
                    sims = [x.replace('|', '') for x in allele.similar.split(', ')]
                    for sim in sims:
                        allele_similars[i][sim] = allele.name

                allele_lookups[i][allele.name] = allele

        common_allele_names = list(allele_names[0] & allele_names[1])

        fo.write('<h2>Comparison of %s and %s</h2>' % (datasets[0], datasets[1]))

        fo.write('<h2>Alleles only in %s</h2>' % datasets[0])
        exc = list(allele_names[0] - allele_names[1])
        exc_com = [('%s (%s in %s)' % (x, allele_similars[1][x], datasets[1]) if x in allele_similars[1] else x) for x in exc]
        fo.write('<br>'.join(exc_com))

        fo.write('<h2>Alleles only in %s</h2>' % datasets[1])
        exc = list(allele_names[1] - allele_names[0])
        exc_com = [('%s (%s in %s)' % (x, allele_similars[0][x], datasets[0]) if x in allele_similars[0] else x) for x in exc]
        fo.write('<br>'.join(exc_com))

        fo.write('<h2>Changed appearance counts</h2>')
        fo.write('<table><tr><th>Allele</th><th>%s</th><th>%s</th></tr>' % (datasets[0], datasets[1]))

        for allele in common_allele_names:
            if allele_lookups[0][allele].appears != allele_lookups[1][allele].appears:
                fo.write('<tr><th>%s</th><th>%d</th><th>%d</th></tr>' % (allele, allele_lookups[0][allele].appears, allele_lookups[1][allele].appears))
        fo.write('</table>')

        fo.write('<h2>Changed max_kdiffs</h2>')
        fo.write('<table><tr><th>Allele</th><th>%s</th><th>%s</th></tr>' % (datasets[0], datasets[1]))

        for allele in common_allele_names:
            if abs(allele_lookups[0][allele].max_kdiff - allele_lookups[1][allele].max_kdiff) > 0.1:
                fo.write('<tr><th>%s</th><th>%.2f</th><th>%.2f</th></tr>' % (allele, allele_lookups[0][allele].max_kdiff, allele_lookups[1][allele].max_kdiff))
        fo.write('</table>')

        fo.write('<h2>Changed confidence levels</h2>')
        fo.write('<table><tr><th>Allele</th><th>%s</th><th>%s</th></tr>' % (datasets[0], datasets[1]))

        for allele in common_allele_names:
            if allele_lookups[0][allele].low_confidence != allele_lookups[1][allele].low_confidence:
                fo.write('<tr><th>%s</th><th>%s</th><th>%s</th></tr>' % (allele, 'low' if allele_lookups[0][allele].low_confidence else 'high', 'low' if allele_lookups[1][allele].low_confidence else 'high'))
        fo.write('</table>')

        fo.write('<h2>Changed number of notes</h2>')
        fo.write('<table><tr><th>Allele</th><th>%s</th><th>%s</th></tr>' % (datasets[0], datasets[1]))

        for allele in common_allele_names:
            c0 = session[0].query(AlleleConfidenceReport).filter(AlleleConfidenceReport.allele_id == allele_lookups[0][allele].id).count()
            c1 = session[1].query(AlleleConfidenceReport).filter(AlleleConfidenceReport.allele_id == allele_lookups[1][allele].id).count()
            if c0 != c1:
                fo.write('<tr><th>%s</th><th>%d</th><th>%d</th></tr>' % (allele, c0, c1))
        fo.write('</table>')

    return send_report(output_path, 'html')


