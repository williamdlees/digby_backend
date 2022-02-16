library(rabhit)
library(dplyr)
library(optparse)

pdf(NULL)  # stop spurious Rplots.pdf being produced

########## VDJbase server ##############

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name", metavar="character"),
  make_option(c("-k", "--Kdiff"), type="character", default=NULL,
              help="The minimal kdiff", metavar="character"),
  make_option(c("-c", "--chain"), type="character", default="IGH",
              help="chain: IGH, IGK, IGL, TRB, TRA"),
  make_option(c("-g", "--gene_order_file"), type="character", default=NULL,
              help="genes listed in desired order (tsv file)")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$output_file)){
  stop("output reference file must be supplied", call.=FALSE)
}

if (is.null(opt$Kdiff)){
  stop("the minimal kdiff vaue must be supplied", call.=FALSE)
}

######### loading data(use melt function) #############

# read genotype table
haplotypes_path <- opt$input_file
output_file <- opt$output_file
haplotypes <- readHaplotypeDb(file=haplotypes_path)
kdiff <- opt$Kdiff

if (!is.null(opt$gene_order_file)){
    gene_order = read.delim(file=opt$gene_order_file, header=FALSE, sep="\t", stringsAsFactors = F)
    gene_order = gene_order$V1
} else {
    gene_order = NULL
}

hapHeatmap(haplotypes, lk_cutoff = kdiff, file = output_file, genes_order=gene_order, chain=opt$chain)
