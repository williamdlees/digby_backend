# Gene Usage Graph
# html new colors

# required libraries
library("ggplot2")
library("tigger")
library("dplyr")
library("mltools")  #need to install for igguest
library('plotly')
#### load variables: ####
library(optparse)
library("vdjbaseVis")

########################

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name", metavar="character"),
  make_option(c("-s", "--sysdata_file"), type="character", default=NULL,
              help="sysdata file name", metavar="character")
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
  stop("sys file name must be supplied", call.=FALSE)
}

######### loading data to data frame (use melt function) #############

# read genotype table
input_file<-opt$input_file
output_file<-opt$output_file
load(opt$sysdata_file)

gene_segment <- read.delim(input_file, header=TRUE, sep="\t",stringsAsFactors = T)
heterozygous_graph <- vdjbaseVis::heterozygousBar_html(gene_segment)

htmlwidgets::saveWidget(heterozygous_graph , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)
