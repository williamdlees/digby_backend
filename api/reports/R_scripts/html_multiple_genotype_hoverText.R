# Or - Create Multiple Genotype Graphs html/pdf("T"/"F").

library(optparse)
library(vdjbasevis)
library(stringr)
library(plotly)

pdf(NULL)  # stop spurious Rplots.pdf being produced


########## VDJbase server ##############

#pdf(NULL)     # This prevents blank Rscript.pdf files being written

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name"),
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


######### VDJbase server- loading data(use melt function) #############

# read genotype table
genotype_path<-opt$input_file
gene_order_file = opt$gene_order_file
output_file<-opt$output_file
chain = opt$chain
data<- read.delim(file= genotype_path ,header=TRUE,sep="\t",stringsAsFactors = F)

if (!is.null(opt$samp)) {
  data$subject <- opt$samp
}

if (!is.null(opt$gene_order_file)){
    gene_order = read.delim(file=opt$gene_order_file, header=FALSE, sep="\t", stringsAsFactors = F)
    gene_order = gene_order$V1
} else {
    gene_order = NULL
}

html_output <- as.logical(opt$is_html)  # for pdf set "F"


######################### Run multiGenotype fuction ##########################################
num_of_subjects <- length(unique(data$subject))
genotype_graph <- vdjbasevis::multipleGenoytpe(geno_table=data, chain=chain, ordered_genes=gene_order, html=html_output)

save.image(file='foo.rdata')

if (html_output) {
  htmlwidgets::saveWidget(genotype_graph[[1]] , file.path(normalizePath(dirname(output_file)),basename(output_file)), background = "white", selfcontained = F)
} else {
  pdf(output_file,onefile = F, width = genotype_graph[[3]], height = genotype_graph[[2]], family = "serif")
  print(genotype_graph)
  dev.off()
  embedFonts(output_file)
}
