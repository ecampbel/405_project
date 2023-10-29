#/bin/envr/R!
library(tidyr)
library(dplyr)

#Goal: add a column to the read data with gene function
gene_function_file <- "uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t") 
gene_function_info <- gene_function_info[order(gene_function_info$Gene.Names),]

read_files <- c("10d1ReadsPerGene.out.tab","10d2ReadsPerGene.out.tab","10d3ReadsPerGene.out.tab",
                "apo1ReadsPerGene.out.tab","apo2ReadsPerGene.out.tab","apo3ReadsPerGene.out.tab")

file <- read.delim("10d1ReadsPerGene.out.tab",header=FALSE,sep="\t") %>% tail(-4) 
colnames(file) <-  c("Gene.Names","exp1","exp2","exp3")

for (rf in read_files) {main(r)}

main <- function(rf){
  file <- read.delim(rf,header=FALSE,sep="\t") %>% tail(-4) 
  colnames(file) <-  c("Gene.Names","Rep1","Rep2","Rep3")
  file_annotated <- left_join(file,gene_function_info,by=c("Gene.Name"))
  write_it(file,rf)
}

write_it <- function(annotated,rf){
  filename <- paste(rf,".annotated",sep="")
  write(annotated, file=filename,sep="\t")
}
