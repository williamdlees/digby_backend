# Gene Usage Graph
# html new colors

# required libraries
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

alleles_appearance <- read.delim(input_file, header=TRUE, sep="\t",stringsAsFactors = T)
alleles_appearance_graph <- vdjbasevis::alleleUsageBar_html(alleles_appearance, chain=opt$chain)

htmlwidgets::saveWidget(alleles_appearance_graph , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)