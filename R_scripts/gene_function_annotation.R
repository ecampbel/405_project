#/bin/envr/R!
library(tidyr)
library(dplyr)

#Read gene function data
gene_function_file <- "uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t") 
gene_function_info <- gene_function_info[order(gene_function_info$Gene.Names),]

rf10d1 <- "10d1ReadsPerGene.out.tab"
rf10d2 <- "10d2ReadsPerGene.out.tab"
rf10d3 <- "10d3ReadsPerGene.out.tab"
apod1 <- "apo1ReadsPerGene.out.tab"
apod2 <- "apo2ReadsPerGene.out.tab"
apod3 <- "apo3ReadsPerGene.out.tab"

#Read and join our ReadsPerOut files
merged_10d <- merge(read_it(rf10d1),read_it(rf10d2),by=c("geneID")) %>% merge(read_it(rf10d3),by=c("geneID"))
merged_apo <- merge(read_it(apod1),read_it(apod2),by=c("geneID")) %>% merge(read_it(apod3),by=c("geneID"))

#For merged annotation files and their associated name (character), joins to gene functional info and writes to csv file
main <- function(merged_file,name){
  file_annotated <- left_join(merged_file,gene_function_info,by=c("geneID"))
  write_it(file_annotated,name)
}

### HELPER FUNCTIONS #####

#Reads a ReadsPerGene.out.tab file, renames its columns.
read_it <- function(rf){
  file <- read.delim(rf,header=FALSE,sep="\t") %>% tail(-4) 
  colnames(file) <-  c("geneID",paste(substring(rf,1,4),"Rep1",sep="_"),paste(substring(rf,1,4),"Rep2",sep="_"),paste(substring(rf,1,4),"Rep3",sep="_"))
  return(file)
  }

#Writes annotated files (with a name you provide)
write_it <- function(annotated,name){
  filename <- paste(name,".annotated.tsv",sep="")
  print(filename)
  write.table(annotated, file=filename,sep="\t",row.names=FALSE,col.names=TRUE)
  }

##### ANNOTATION OF NORMALISED FILE
normalised <- read.csv("echlorotica_TPM.csv") 
names(gene_function_info)[names(gene_function_info) == 'Gene.Names'] <- 'geneID'
normalsied_annotated <- left_join(normalised,gene_function_info,by=c("geneID"))
