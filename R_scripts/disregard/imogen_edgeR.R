### EdgeR running! Trial
library(edgeR)
library(tibble)
# setwd("../")
source("/R_scripts/gene_function_annotation.R")

gene_function_file <- "../uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t") 
gene_function_info <- gene_function_info[order(gene_function_info$Gene.Names),]

#Helpers
#Reads a ReadsPerGene.out.tab file, renames its columns.
read_it <- function(rf){
  file <- read.delim(rf,header=FALSE,sep="\t") %>% tail(-4) 
  colnames(file) <-  c("geneID",paste(substring(rf,1,4),"Rep1",sep="_"),paste(substring(rf,1,4),"Rep2",sep="_"),paste(substring(rf,1,4),"Rep3",sep="_"))
  return(file)}
#Writes annotated files (with a name you provide)
write_it <- function(annotated,name){
  filename <- paste(name,".annotated.tsv",sep="")
  print(filename)
  write.table(annotated, file=filename,sep="\t",row.names=FALSE,col.names=TRUE)}

rf10d1 <- "../rawCounts/10d1ReadsPerGene.out.tab"
rf10d2 <- "../rawCounts/10d2ReadsPerGene.out.tab"
rf10d3 <- "../rawCounts/10d3ReadsPerGene.out.tab"
apod1 <- "../rawCounts/apo1ReadsPerGene.out.tab"
apod2 <- "../rawCounts/apo2ReadsPerGene.out.tab"
apod3 <- "../rawCounts/apo3ReadsPerGene.out.tab"

#Read and join our ReadsPerOut files
merged_10d <- merge(read_it(rf10d1)[,1:2],read_it(rf10d2)[,1:2],by=c("geneID")) %>% merge(read_it(rf10d3)[,1:2],by=c("geneID"))
merged_apo <- merge(read_it(apod1)[,1:2],read_it(apod2)[,1:2],by=c("geneID")) %>% merge(read_it(apod3)[,1:2],by=c("geneID"))
colnames(merged_10d) <- c("geneID","1","2","3")
colnames(merged_apo) <- c("geneID","1","2","3")
x <- merge(merged_10d,merged_apo,by="geneID") %>% column_to_rownames("geneID") 
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
keep<-filterByExpr(y)
y<-y[keep,,keep.lib.sizes=FALSE]
y<-normLibSizes(y)
design<-model.matrix(~group)
y<-estimateDisp(y,design)

## Quasi-likehood F-tests
fit<-glmQLFit(y,design)
qlf<-glmQLFTest(fit,coef=2)
topTags(qlf)



