# Standalone script to run a set of reports

import requests
import time

#reports_url = "http://localhost:5000"
reports_url = "https://vdjbase.org/admin/"


report_requests = [
    # AIRR-seq only

    ("genotype list, human IGH, 7 samples, pdf", "/api/reports/reports/run/rep_genotype?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%7B%22field%22:%22sample_name%22,%22op%22:%22in%22,%22value%22:%5B%22P10_I3_S1%22,%22P10_I4_S2%22,%22P14_I4_S1%22,%22P14_I7_S1%22,%22P1_I40_S1%22,%22P1_I43_S1%22,%22P1_I63_S1%22%5D%7D%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("genotype list, human IGH, 7 samples, V only, locus order, html", "/api/reports/reports/run/rep_genotype?format=html&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%7B%22field%22:%22sample_name%22,%22op%22:%22in%22,%22value%22:%5B%22P10_I3_S1%22,%22P10_I4_S2%22,%22P14_I4_S1%22,%22P14_I7_S1%22,%22P1_I40_S1%22,%22P1_I43_S1%22,%22P1_I63_S1%22%5D%7D%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22sort_order%22:%22Locus%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%5B%22IGHV%22%5D,%22f_genes%22:%22%22%7D"),
    ("genotype heatmap, human IGH, all samples, J only, pdf", "/api/reports/reports/run/rep_genotype_heatmap?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22f_kdiff%22:%22%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%5B%22IGHJ%22%5D,%22f_genes%22:%22%22%7D"),
    ("genotype heatmap, human IGH, all samples, html", "/api/reports/reports/run/rep_genotype_heatmap?format=html&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22f_kdiff%22:%22%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("haplotype report, human IGH, 6 samples, pdf, kdiff >= 1", "/api/reports/reports/run/haplo_heatmap?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%7B%22field%22:%22haplotypes%22,%22op%22:%22in%22,%22value%22:%5B%22J6-02_04%22%5D%7D%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22haplo_gene%22:%22IGHD2-21-01_02%22,%22f_kdiff%22:1,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele appearance, human IGH, all samples, pdf", "/api/reports/reports/run/allele_appearance?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele usage, human IGH, all samples", "/api/reports/reports/run/allele_usage?format=html&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_kdiff%22:%22%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),

    ("allele appearance, human TRB, all samples, pdf", "/api/reports/reports/run/allele_appearance?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=TRB&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele appearance, human TRB, inc pseudo, all samples","/api/reports/reports/run/allele_appearance?format=xls&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=TRB&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:true,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("heterozygosity, human TRB, all samples, html", "/api/reports/reports/run/heterozygosity?format=html&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=TRB&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22sort_order%22:%22Alphabetic%22,%22f_kdiff%22:%22%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),


    ("allele usage, human IGK, all samples, pdf", "/api/reports/reports/run/allele_usage?format=html&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGK&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_kdiff%22:%22%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),

    ("gene frequencies, human IGL, all samples, pdf", "/api/reports/reports/run/gene_frequencies?format=pdf&species=Human&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGL&rep_filters=%5B%5D&params=%7B%22single_sample%22:%22All%20Available%22,%22calculate_by%22:%22Number%20of%20Clones%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),

    ("genotype list, rhesus IGH, 4 samples, pdf", "/api/reports/reports/run/rep_genotype?format=pdf&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%7B%22field%22:%22sample_name%22,%22op%22:%22in%22,%22value%22:%5B%22P22_I12_S1%22,%22P22_I15_S1%22,%22P22_I18_S1%22,%22P22_I10_S1%22%5D%7D%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("genotype list, rhesus IGH, 4 samples, V only, html", "/api/reports/reports/run/rep_genotype?format=html&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%7B%22field%22:%22sample_name%22,%22op%22:%22in%22,%22value%22:%5B%22P22_I12_S1%22,%22P22_I15_S1%22,%22P22_I18_S1%22,%22P22_I10_S1%22%5D%7D%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%5B%22IGHV%22%5D,%22f_genes%22:%22%22%7D"),
    ("genotype heatmap, rhesus IGH, all samples, pdf", "/api/reports/reports/run/rep_genotype_heatmap?format=pdf&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22f_kdiff%22:%22%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("genotype heatmap, rhesus IGH, all samples, locus order, pseudo, html", "/api/reports/reports/run/rep_genotype_heatmap?format=html&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22per_sample%22:%22Each%20sample%22,%22f_kdiff%22:%22%22,%22sort_order%22:%22Locus%22,%22f_pseudo_genes%22:true,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele appearance, rhesus IGH, all samples, pdf", "/api/reports/reports/run/allele_appearance?format=pdf&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele appearance, rhesus IGH, all samples, xls", "/api/reports/reports/run/allele_appearance?format=xls&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("heterozygosity, rhesus IGH, all samples, pdf", "/api/reports/reports/run/heterozygosity?format=html&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22sort_order%22:%22Alphabetic%22,%22f_kdiff%22:%22%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
    ("allele usage, rhesus IGH, all samples, V only, pdf", "/api/reports/reports/run/allele_usage?format=html&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22ambiguous_alleles%22:%22Exclude%22,%22novel_alleles%22:%22Include%22,%22sort_order%22:%22Alphabetic%22,%22f_kdiff%22:%22%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%5B%22IGHV%22%5D,%22f_genes%22:%22%22%7D"),
    ("gene frequencies, rhesus IGH, all samples", "/api/reports/reports/run/gene_frequencies?format=pdf&species=Rhesus%20Macaque&genomic_datasets=&genomic_filters=%5B%5D&rep_datasets=IGH&rep_filters=%5B%5D&params=%7B%22single_sample%22:%22All%20Available%22,%22calculate_by%22:%22Number%20of%20Clones%22,%22sort_order%22:%22Alphabetic%22,%22f_pseudo_genes%22:%22%22,%22f_gene_types%22:%22%22,%22f_genes%22:%22%22%7D"),
]


for report_desc, report_request in report_requests:
    print(f"\n--- {report_desc} ---")
    response = requests.get(reports_url + report_request)

    if response.status_code != 200:
        print("Error: " + str(response.status_code))
        exit()

    response_text = response.json()
    print(response_text)
    report_id = response_text['id']
    status = response_text['status']

    while status == 'queued' or status == 'PENDING':
        time.sleep(1)
        response = requests.get(reports_url + "/api/reports/reports/status/" + report_id)
        response_text = response.json()
        print(response_text)
        status = response_text['status']

print("All reports complete")



