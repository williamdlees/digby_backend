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
              help="sysdata file name", metavar="character"),
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
  stop("sys file name must be supplied", call.=FALSE)
}

if (is.null(opt$is_html)){
  stop("type of file must be supplied", call.=FALSE)
}
######### loading data to data frame (use melt function) #############

# read genotype table
input_file<-opt$input_file
output_file<-opt$output_file
load(opt$sysdata_file)

html_output <- opt$is_html  # for pdf set "F"
html_output <- ifelse(html_output == "T", TRUE, FALSE)

frequencies <- read.delim(input_file, header=TRUE, sep="\t",stringsAsFactors = T)
gene_usage_graph <- geneUsage(frequencies, plot_style="ggplot")

if (html_output) {
  htmlwidgets::saveWidget(gene_usage_graph , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)
} else {
  pdf(output_file,onefile = F, width = 12.5, height = 15, family = "serif")
  plot(gene_usage_graph)
  dev.off()
  embedFonts(output_file)
}



