#/bin/envr/R!
# install.packages("reshape")
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

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
normalised_annotated <- left_join(normalised,gene_function_info,by=c("geneID"))
write_it(normalised_annotated,"normalised_annotated")

###### Imogen messes around ################################################################################################
analysis <- matrix(nrow=nrow(normalised_annotated),ncol=10)
for (n in 1:6){analysis[,n] <- normalised_annotated[,n+2]}
for (n in 1:nrow(analysis)){analysis[n,7] <- mean(analysis[n,1:3])}  #d10 means
for (n in 1:nrow(analysis)){analysis[n,8] <- sd(analysis[n,1:3])}   #d10 stdev
for (n in 1:nrow(analysis)){analysis[n,9] <- mean(analysis[n,4:6])}  #apo means
for (n in 1:nrow(analysis)){analysis[n,10] <- sd(analysis[n,4:6])}  #apo stdev
analysis_df <- as.data.frame(analysis) 
colnames(analysis_df) <- c("t10R1","t10R2","t10R3","apoR1","apoR2","apoR3","t10_av","t10_sd","apo_av","apo_sd")

analysis_df2 <- as.data.frame(matrix(nrow=nrow(normalised_annotated),ncol=6)) 
for (n in 1:6){analysis_df2[,n] <- normalised_annotated[,n+2]}
colnames(analysis_df2) <- c("t10R1","t10R2","t10R3","apoR1","apoR2","apoR3")
analysis_df2 <- melt(analysis_df2)
ggplot(analysis_df2,mapping=aes(variable,value))+geom_boxplot()+scale_y_log10() #recreates gene count 
ggplot(analysis_df2,mapping=aes(value))+geom_histogram()+scale_x_log10() #as histogram, logarithmic
