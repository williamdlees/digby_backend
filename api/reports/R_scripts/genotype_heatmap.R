library(optparse)
library(vdjbaseVis)

########## VDJbase server ##############

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name"),
  make_option(c("-s", "--sysdata_file"), type="character", default=NULL,
              help="sysdata file name"),
  make_option(c("-k", "--Kdiff"), type="numeric", default=NULL,
              help="The minimal kdiff"),
  make_option(c("-t", "--is_html"), type="character", default=NULL,
              help="type of file F - pdf T - html")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$output_file)){
  stop("output reference file must be supplied", call.=FALSE)
}

if (is.null(opt$sysdata_file)){
  stop("sysdata file must be supplied", call.=FALSE)
}

if (is.null(opt$Kdiff)){
  stop("the minimal kdiff vaue must be supplied", call.=FALSE)
}

if (is.null(opt$is_html)){
  stop("type of file must be supplied", call.=FALSE)
}
######### loading data(use melt function) #############
# read genotype table
genotype_path<-opt$input_file
output_file<-opt$output_file
genotypes <- read.delim(file= genotype_path ,header=TRUE,sep="\t",stringsAsFactors = F)

kdiff <- as.numeric(opt$Kdiff)
html_output <- as.logical(opt$is_html) # for pdf set "F"

if (html_output) {
  #num_of_genes <- length(unique(genotypes$GENE))
  #width <- num_of_genes * 0.24 + 1.5
  vdjbaseVis::genoHeatmap_html(genotypes, lk_cutoff = kdiff,
                               file = file.path(normalizePath(dirname(output_file)),basename(output_file)))
} else {
  vdjbaseVis::genoHeatmap(genotypes, lk_cutoff = kdiff,
                          file = file.path(normalizePath(dirname(output_file)),basename(output_file)))
}
