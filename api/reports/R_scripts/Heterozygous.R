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
library("vdjbasevis")

pdf(NULL)  # stop spurious Rplots.pdf being produced


########################

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name", metavar="character"),
  make_option(c("-c", "--chain"), type="character", default="IGH",
              help="chain: IGH, IGK, IGL, TRB, TRA")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$output_file)){
  stop("output reference file must be supplied", call.=FALSE)
}

######### loading data to data frame (use melt function) #############

# read genotype table
input_file<-opt$input_file
output_file<-opt$output_file

gene_segment <- read.delim(input_file, header=TRUE, sep="\t",stringsAsFactors = T)
#gene_segment <- data.frame(GENE = c("V3-3",'V1-2','D2-8','D3-16','J4','J6'), HM = c(20,60,55,7,30,0) , HT = c(80,40,45,93,0,45))
heterozygous_graph <- vdjbasevis::heterozygousBar_html(gene_segment, chain=opt$chain)

htmlwidgets::saveWidget(heterozygous_graph , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)
