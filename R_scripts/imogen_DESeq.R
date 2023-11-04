#/bin/envr/R!
library(DESeq2)
library(apeglm)
library(tibble)
sampleTable <- data.frame(Sample = c("Sample1", "Sample2", "Sample3"),
  Condition = c("d10", "d10", "d10", "apo","apo","apo"))

# Create a count matrix with genes as rows and samples as columns
merged_10d <- merge(read_it(rf10d1)[,1:2],read_it(rf10d2)[,1:2],by=c("geneID")) %>% merge(read_it(rf10d3)[,1:2],by=c("geneID"))
merged_apo <- merge(read_it(apod1)[,1:2],read_it(apod2)[,1:2],by=c("geneID")) %>% merge(read_it(apod3)[,1:2],by=c("geneID"))
colnames(merged_10d) <- c("geneID","d10","d10","d10")
colnames(merged_apo) <- c("geneID","apo","apo","apo")
countMatrix <- as.matrix(merge(merged_10d,merged_apo,by="geneID") %>% column_to_rownames("geneID") )

### Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleTable, design = ~ Sample + Condition) %>% DESeq()
res <- results(dds, name="Condition_d10_vs_apo")
resLFC <- lfcShrink(dds, coef="Condition_d10_vs_apo", type="apeglm")
resOrdered <- na.omit(res[order(res$pvalue),])

### Get summary data
summary(res)

### Create MA-plot (log2 fold changes over mean of normalised counts)
plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange) #run to activate identifier <-  click on points to get gene ID
rownames(res)[idx]

### Create PCA plots
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Condition","Sample"))
