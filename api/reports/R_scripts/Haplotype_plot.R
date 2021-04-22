# Haplotype Graph
# html new colors

# required libraries
library('plotly')
#### load variables: ####
library(optparse)
library('rabhit')

########################

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name", metavar="character"),
  make_option(c("-t", "--is_html"), type="character", default=NULL,
              help="type of file F - pdf T - html"),
  make_option("--samp", type="character", default=NULL,
              help="Sample name"),
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

if (is.null(opt$is_html)){
  stop("type of file must be supplied", call.=FALSE)
}

######### loading data to data frame (use melt function) #############

# read genotype table
haplotype_path<-opt$input_file
output_file<-opt$output_file
sample_name <- opt$sample_name

html_output <- opt$is_html  # for pdf set "F"
html_output <- ifelse(html_output == "T", TRUE, FALSE)

haplo_db_J6 <- read.delim(haplotype_path, header=TRUE, sep="\t",stringsAsFactors = F, colClasses = "character")
if (!is.null(opt$samp)) {
  haplo_db_J6$subject <- opt$samp
}
names(haplo_db_J6)[c(1,2,5:ncol(haplo_db_J6))] <- tolower(names(haplo_db_J6)[c(1,2,5:ncol(haplo_db_J6))])

haplotype_graph <- plotHaplotype(haplo_db_J6, html_output = html_output, text_size = 11)

if (html_output) {
  htmlwidgets::saveWidget(haplotype_graph , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)
} else {
  pdf(output_file,onefile = F, width = 12.5, height = 15, family = "serif")
  print(haplotype_graph)
  dev.off()
  embedFonts(output_file)
}



