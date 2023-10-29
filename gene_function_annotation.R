#/bin/envr/R!
library(tidyr)
library(dplyr)

#Goal: add a column to the read data with gene function
gene_function_file <- "uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t") 
gene_function_info <- gene_function_info[order(gene_function_info$Gene.Names),]
gene_function_info <- gene_function_info["EGW08"%in%gene_function_info$Gene.Names==TRUE,]

read_files <- c("10d1ReadsPerGene.out.tab","10d2ReadsPerGene.out.tab","10d3ReadsPerGene.out.tab",
                "apo1ReadsPerGene.out.tab","apo2ReadsPerGene.out.tab","apo3ReadsPerGene.out.tab")

file <- read.delim("10d1ReadsPerGene.out.tab",header=FALSE,sep="\t") %>% tail(-4) 
colnames(file) <-  c("Gene.Names","exp1","exp2","exp3")

for (rf in read_files) {main(r)}

main <- function(read_file){
  file <- read.delim(read_file,header=FALSE,sep="\t") %>% tail(-4) 
  colnames(file) <-  c("Gene.Names","exp1","exp2","exp3")
  file_annotated <- annotate(file)
  write_it(file)
}

annotate <- function(file){
  
}

write_it <- function(file){
  
}