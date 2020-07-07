# html new colors

# required libraries
library("ggplot2")
library("tigger")
library("dplyr")
library("mltools")  #need to install for igguest
library("cowplot")
library("gridExtra")
library("grid")
library('plotly')
#### load variables: ####
library(optparse)

########################


option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, 
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="graph.pdf file name", metavar="character"),
  make_option(c("-s", "--sysdata_file"), type="character", default=NULL,
              help="sysdata file name", metavar="character"),
  make_option(c("-p", "--sample_name"), type="character", default=NULL,
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

if (is.null(opt$sample_name)){
  stop("sample name must be supplied", call.=FALSE)
}

if (is.null(opt$sysdata_file)){
  stop("sysdata file must be supplied", call.=FALSE)
}

######### loading data to data frame (use melt function) #############


# read genotype table
genotype_path<-opt$input_file
output_file<-opt$output_file
sys_file <- opt$sysdata_file
sample_name <- opt$sample_name
# # for testing 
# genotype_path<-"/home/aviv/Data/P1_Celiac_Data/genotypes/P1_I100_S1_geno_H.tab"
# output_file<-"/home/aviv/test.html"


# load the "sysdata"
# load("/home/aviv/Rscripts/sysdata.rda")
load(sys_file)


removeIGH = TRUE
text_size = 12 # or 14
plotYaxis=TRUE
html_output=FALSE
chain='IGH'
gene_sort='position'

########################

######## read table #################  



# # read genotype table
data<- read.delim(file= genotype_path ,header=TRUE,sep="\t",stringsAsFactors = F)
gen_table <- data
gen_table <- data%>% select(GENE,ALLELES,COUNTS,TOTAL,K_DIFF,GENOTYPED_ALLELES)


######## end: read table #################  



#~~~~~~~~~~~~~~~~~~~Creat genotype plot (remove pseudo genes + insert tranper color for new allele) : ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

GENE.loc.tmp <- GENE.loc[[chain]]  # take just IGH genes
genotype = data
SAMP=sample_name
# for testing
# SAMP <- "P1_I100_S1"
alleles = strsplit(genotype$GENOTYPED_ALLELES, ",")

# split all multi alleles to number of raws

genotype <- genotype[, -which(names(genotype) == "NOTE")]
geno2 <- genotype

r = 1
for (g in 1:nrow(genotype)) {
  for (a in 1:length(alleles[[g]])) {
    geno2[r, ] = genotype[g, ]
    geno2[r, ]$ALLELES = alleles[[g]][a]
    r = r + 1
  }
}

kval.df = subset(genotype,select=c(GENE,GENOTYPED_ALLELES,K_DIFF)) # Kval per gene

# sort GENE coulmn 
if(gene_sort=='name'){
  geno2$GENE = factor(geno2$GENE, levels = rev(sortAlleles(unique(geno2$GENE), method = gene_sort)))
  kval.df$GENE = factor(kval.df$GENE, levels = rev(sortAlleles(unique(kval.df$GENE), method = gene_sort)))
} else {
  
  names(GENE.loc.tmp) <- GENE.loc.tmp
  geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp),ordered = TRUE)
  kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))
  
  #remove rows with NA
  geno2 <- na.omit(geno2)
  kval.df <- na.omit(kval.df)
  
  if(removeIGH){
    GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
    names(GENE.loc.tmp) <- GENE.loc.tmp
    geno2$GENE <- gsub('IG[H|K|L]','',geno2$GENE)
    kval.df$GENE <- gsub('IG[H|K|L]','',kval.df$GENE)
    geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
    kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))
  } else {
    names(GENE.loc.tmp) <- GENE.loc.tmp
    
    geno2$GENE = factor(geno2$GENE, levels = rev(GENE.loc.tmp))
    kval.df$GENE = factor(kval.df$GENE, levels = rev(GENE.loc.tmp))
    
  }
  
}


AlleleCol <- c(sort(grep('[012]',unique(geno2$ALLELES),value = T,perl = T)),'Unk')


## Alleles plot 
p = ggplot(geno2, aes(x = GENE, fill = factor(ALLELES,levels=AlleleCol))) + theme_bw() + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = text_size), strip.background = element_blank(), 
        strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black"),
        panel.spacing = unit(0, "cm"),
        strip.switch.pad.grid = unit(0, "cm"),
        plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill") + 
  coord_flip() + xlab("Gene") + ylab("") 




# No need? 
if(!plotYaxis){
    p=p+theme(axis.text.y=element_blank())
  }
  
  
  if(is.null(ALLELE_PALETTE)){
    AlleleCol <- c(sort(grep('[012]',unique(geno2$ALLELES),value = T,perl = T)),'Unk')
    tmp.col <- scale_fill_hue(h.start = 10, h=c(10, 270))
    len.col <- length(AlleleCol)-2
    names(AlleleCol) <-  c(tmp.col$palette(len.col),'#dedede','#6d6d6d')
    p = p +scale_fill_manual(values=(names(AlleleCol)),name='ALLELES')
  } else {
    
    AlleleCol <- grep('[012]',unique(geno2$ALLELES),value = T,perl = T)
    AlleleCol.tmp <- sort(unique(sapply(strsplit(AlleleCol,'_'),'[',1)))
    tmp.col <- ALLELE_PALETTE[AlleleCol.tmp]
    
    novels <- grep('_',AlleleCol,value = T)
    if(length(novels) > 0){
      novels.col <- ALLELE_PALETTE[sapply(strsplit(novels,'_'),'[',1)]
      names(novels.col) <- novels 
      alleles.comb <- c(tmp.col,novels.col)[order(names(c(tmp.col,novels.col)))]
    } else {
      alleles.comb <- c(tmp.col)[order(names(c(tmp.col)))]
      
    }
    
    AlleleCol<- names(c(alleles.comb,Unk='#dedede'))
    names(AlleleCol) <- c(alleles.comb,Unk='#dedede')
    
    
    
    transper <- sapply(AlleleCol,function(x){if(grepl('_',x)){mom_allele <- strsplit(x,'_')[[1]][1];
    all_novel <- grep(paste0(mom_allele,'_'),AlleleCol,value=T);
    if(length(all_novel)==1){return(0.5)};
    if(length(all_novel)==2){m=which(all_novel==x);return(ifelse(m==1,0.6,0.3))}
    if(length(all_novel)==3){m=which(all_novel==x);if(m==1){return(0.6)} ; return(ifelse(m==2,0.4,0.2))}
    } else (1)})
    names(transper) <- AlleleCol
    
    #? remove 'mother' allele if added (when there is no germline allele but there is a novel)
    AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep('[012]',unique(geno2$ALLELES),value = T,perl = T)),'Unk')]
    transper <- transper[names(transper) %in% AlleleCol ]
    p = p +scale_fill_manual(values=alpha(names(AlleleCol),transper),name='Alleles')
    
  }
  
  ## plot K values   
  kval.df$K_GROUPED <- bin_data(kval.df$K, bins=c(0, 1,2,3,4,5,10,20,50,Inf), binType = "explicit")
  #? blue_values <-c("#9ECAE1","#B0C4DE","#4292C6","#2171B5", "#08519C" ,"#00009B","#08306B",)
  pk <- ggplot(kval.df, aes(x = GENE, fill = K_GROUPED))+theme_bw() + 
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          text = element_text(size = text_size), strip.background = element_blank(), 
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0, "cm"),
          strip.switch.pad.grid = unit(0, "cm"),
          plot.margin = unit(c(0.25, 0, 0.2, 0), "cm")) + geom_bar(position = "fill")  +
    coord_flip() +xlab("")+ ylab("") 
  
#------------------------------------------create html ------------------------------------------------------------    
  
  blue_values <-c("#4682B433","#B0C4DE","#9ECAE1","#4292C6","#2171B5", "#08519C" ,"#00008B","#00009B","#08306B")
  #blue_values <-c("#4682B433","#B0C4DE","#9ECAE1","#4292C6","#00CCFF", "#08519C" ,"#0000FF","#000080","#08306B")
  #blue_values <- rainbow(16,s=c(0.5:0.6))
  pk <- pk + scale_fill_manual(values =blue_values)  # +ggtitle(paste0("Genotype Graph for: ",SAMP))
  
  pk.l <- ggplotly(pk,height = 1000,width = 700) %>% plotly::layout(showlegend=TRUE)
  
  p2.l <- ggplotly(p,height = 1300,width = 1000) %>% plotly::layout(margin=list(b=50,t=50),
                                                                    yaxis = list(title = paste0(c(rep("&nbsp;", 3),
                                                                                                  "Gene",
                                                                                                  rep("&nbsp;", 3),
                                                                                                  rep("\n&nbsp;", 1)),
                                                                                                collapse = "")),
                                                                    showlegend=TRUE)
  
  # title of legends      
  p.l.c <- suppressWarnings(subplot(p2.l,pk.l,widths = c(0.4,0.1),shareY = T,titleX = TRUE,margin = 0.01,which_layout = 1))
  
  p.l.c$x$layout$annotations[[2]]$text = "log<sub>10</sub>(K)-------------"
  p.l.c$x$layout$annotations[[2]]$xanchor="center"
  p.l.c$x$layout$annotations[[2]]$y=0.96-0.0234*(length(AlleleCol)-0.5)#0.52
  p.l.c$x$layout$annotations[[2]]$x=1.02
  
  p.l.c$x$layout$annotations[[1]] <- p.l.c$x$layout$annotations[[2]]
  p.l.c$x$layout$annotations[[1]]$text = "Alleles-------------"
  p.l.c$x$layout$annotations[[1]]$y=0.96
  p.l.c$x$layout$annotations[[1]]$x=1.02
  p.l.c$x$layout$annotations[[1]]$legendTitle=FALSE
  
  # title for graph [name of sample]
  p.l.c$x$layout$title= paste0("Sample Name : ",SAMP)
  print(GENE.loc.tmp)
  
  # save html 
  #htmlwidgets::saveWidget(p.l.c , paste0(getwd(),"/interactive/Genotype_html_",SAMP,".html"),selfcontained = T)
  htmlwidgets::saveWidget(p.l.c , file.path(normalizePath(dirname(output_file)),basename(output_file)),selfcontained = F)
  
  
############# combine Html plotly ############### 
  # install.packages('shiny')
  # require(shiny)
  # install.packages('manipulateWidget')
  # require(manipulateWidget)
  # list_plots = list(p.l.c,p.l.c)
  # combineWidgets(list=list_plots,nrow=1)  
############# end: combine Html plotly ############### 
  