# Check that annotations in genomic directories match BAM contents
# Run from samples/Genomic
# Must be run under linux, requires pysam

import glob
import os
import csv
import pysam


for bamfile in glob.glob('**/*.bam', recursive=True):
    annot_file = bamfile.replace('.bam', '.csv')

    if not os.path.exists(annot_file):
        print(f"Annotation file {annot_file} not found")
        continue

    with open(annot_file) as f:
        recs = {}
        annot_reader = list(csv.DictReader(f))
        annot_contigs = [x['contig'] for x in annot_reader]

        samfile = pysam.AlignmentFile(bamfile, "rb")
        samfile_contigs = [x.query_name for x in samfile.fetch()]

        if set(annot_contigs) - set(samfile_contigs):
            print(f"Annotations in {annot_file} not found in {bamfile}")

        annot_haps = [x['haplotype'] for x in annot_reader]
        if len(annot_haps) < 3:
            print(f"Annotation file {annot_file} has less than 3 haplotypes")
            continue
