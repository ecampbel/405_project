#/bin/envr/R!
library(tidyr)
library(dplyr)

#Goal: add a column to the read data with gene function
gene_function_file <- "uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t")

read_files <- c("10d1ReadsPerGene.out.tab","10d2ReadsPerGene.out.tab","10d3ReadsPerGene.out.tab",
                "apo1ReadsPerGene.out.tab","apo2ReadsPerGene.out.tab","apo3ReadsPerGene.out.tab")

file <- read.delim("10d1ReadsPerGene.out.tab",header=FALSE,sep="\t")
file <- file %>% tail(-4), colnames()

for (rf in read_files) {annotate(r)}

annotate <- function(read_file){
  file <- read.delim(read_file)
  
}