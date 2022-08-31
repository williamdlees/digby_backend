# Format definitions for VDJbase pipeline files

new_format = True

if new_format:
    GENE_COLUMN = "gene"
    COUNTS_COLUMN = "counts"
    TOTAL_COLUMN = "total"
    GENOTYPE_KDIFF_COLUMN = "k_diff"
    ALLELES_COLUMN = "alleles"
    INT_SEP = ";"
else:
    GENE_COLUMN = "GENE"
    COUNTS_COLUMN = "COUNTS"
    TOTAL_COLUMN = "TOTAL"
    GENOTYPE_KDIFF_COLUMN = "K_DIFF"
    ALLELES_COLUMN = "ALLELES"
    INT_SEP = ","


# unvarying between formats
GENOTYPED_ALLELES_COLUMN = "GENOTYPED_ALLELES"

FREQ_BY_SEQ = "Freq_by_Seq"
FREQ_BY_CLONE = "Freq_by_Clone"

HAPLO_KDIFF_COLUMN = "K"
