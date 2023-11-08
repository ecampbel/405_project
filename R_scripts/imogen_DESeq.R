#/bin/envr/R!
library(DESeq2)
library(apeglm)
library(tibble)
library(EnhancedVolcano)

source("gene_function_annotation.R")
sampleTable <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"),
  Condition = c("d10", "d10", "d10", "apo","apo","apo"))

#Read gene function data
gene_function_file <- "../uniprotkb_taxonomy_id_188477_2023_10_29.tsv"
gene_function_info <- read.delim(gene_function_file,header=TRUE,sep="\t") 
gene_function_info <- gene_function_info[order(gene_function_info$Gene.Names),]
# normalised <- read.csv("../echlorotica_TPM.csv") 
# names(gene_function_info)[names(gene_function_info) == 'Gene.Names'] <- 'geneID'
# normalised_annotated <- left_join(normalised,gene_function_info,by=c("geneID"))

# Create a count matrix with genes as rows and samples as columns
merged_10d <- merge(read_it(rf10d1)[,1:2],read_it(rf10d2)[,1:2],by=c("geneID")) %>% merge(read_it(rf10d3)[,1:2],by=c("geneID"))
merged_apo <- merge(read_it(apod1)[,1:2],read_it(apod2)[,1:2],by=c("geneID")) %>% merge(read_it(apod3)[,1:2],by=c("geneID"))
colnames(merged_10d) <- c("geneID","d10","d10","d10")
colnames(merged_apo) <- c("geneID","apo","apo","apo")
countMatrix <- as.matrix(merge(merged_10d,merged_apo,by="geneID") %>% column_to_rownames("geneID") )

# Create a count matrix with genes as rows and samples as columns
# countMatrix <- as.matrix(normalised[,3:8])
rownames(countMatrix) <- normalised[,1]
colnames(countMatrix) <- c("d10","d10","d10","apo","apo","apo")

### Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(countMatrix), colData = sampleTable, design = ~ Condition) 
dds <- DESeq(dds)
res <- results(dds, name="Condition_d10_vs_apo")
resLFC <- lfcShrink(dds, coef="Condition_d10_vs_apo", type="apeglm")
resOrdered <- na.omit(res[order(res$pvalue),])

### Get summary data
# summary(res)

### TEST: removal of batch effects
# mat <- assay(vsd)
# mm <- model.matrix(~Condition, colData(vsd))
# mat <- limma::removeBatchEffect(mat, batch=vsd$Sample, design=mm)
# assay(vsd) <- mat
# plotPCA(vsd,intgroup=c("Condition"))

### Create MA-plot (log2 fold changes over mean of normalised counts)
DESeq2::plotMA(resLFC, main="title", colNonSig = "gray60",colSig = "blue",  colLine = "grey40") 
idx <- identify(res$baseMean, res$log2FoldChange) #run to activate identifier <-  click on points to get gene ID
rownames(res)[idx]

### Create PCA plots
# test<-rlogTransformation(dds,blind=TRUE, fitType='local')
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
DESeq2::plotPCA(vsd, intgroup=c("Condition"))

### Create Volcano plots using EnhancedVolcano
EnhancedVolcano(resLFC, lab=rownames(res), x='log2FoldChange', y='pvalue',
                pCutoff=10e-32, FCcutoff=0.5, pointSize=2.0, labSize = 4.0)
                #can select labels of interest if we wish
############################################################################################################
### Filter for significance -- not necessary.
# sig <- na.omit(res[na.omit(res$padj)<0.05,])
# sig_up <- as.data.frame(sig[sig$log2FoldChange>0,])
# sig_down <- as.data.frame(sig[sig$log2FoldChange<0,])
# go_data <- gene_function_info[,c(4,8)]
# 
# sig_up_annotated <- left_join(as.data.frame(rownames_to_column(sig_up)),go_data,by=c("rowname"="Gene.Names"))
# sig_up_annotated <- sig_up_annotated[sig_up_annotated$Gene.Ontology.IDs!="",] 
# sig_down_annotated <- left_join(as.data.frame(rownames_to_column(sig_down)),go_data,by=c("rowname"="Gene.Names"))
# sig_down_annotated <- sig_down_annotated[sig_down_annotated$Gene.Ontology.IDs!="",] 

