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
  make_option(c("-s", "--sysdata_file"), type="character", default=NULL,
              help="sysdata file name", metavar="character"),
  make_option(c("-k", "--Kdiff"), type="character", default=NULL,
              help="The minimal kdiff", metavar="character")
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

######### loading data(use melt function) #############

# read genotype table
haplotypes_path<-opt$input_file
output_file<-opt$output_file
haplotypes <- read.delim(file=haplotypes_path ,header=TRUE,sep="\t",stringsAsFactors = F)
names(haplotypes)[c(1,2,5:ncol(haplotypes))] <- tolower(names(haplotypes)[c(1,2,5:ncol(haplotypes))])
kdiff <- opt$Kdiff

# load the "sysdata"
load(opt$sysdata_file)

hapHeatmap(haplotypes, lk_cutoff = kdiff, file = output_file)
