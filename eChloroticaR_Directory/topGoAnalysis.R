# Installing topGO
# BiocManager::install("topGO")

# load topGO
library("topGO")

#filter DESeq for significance
sig <- na.omit(res)
sig <- sig[sig$padj<0.05,]
sig_up <- as.data.frame(sig[sig$log2FoldChange>=1,])
sig_up <- tibble::rownames_to_column(sig_up, "gene_id")
sig_down <- as.data.frame(sig[sig$log2FoldChange<=1,])
sig_down <- tibble::rownames_to_column(sig_down, "gene_id")

# read in the "gene universe" file
go_data <- normalised_annotated[,c(1,19)]
write_tsv(go_data,"echlor_GOIDs.tsv", col_names = FALSE)
geneID2GO <- readMappings("echlor_GOIDs.tsv")
geneUniverse <- names(geneID2GO)

# read in the genes of interest

# if feeding directly from DESeq2, then use the following instead
upregulated_genes <- as.character(sig_up$gene_id)
downregulated_genes <- as.character(sig_down$gene_id)

# factor the names
up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_genes))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_genes))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

# build the GOdata object in topGO for upregulated
up_GO_data <- new("topGOdata",
                  description = "EChlor_apo_d10",
                  ontology = "MF",
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

# build the GOdata object in topGO for downregulated
down_GO_data <- new("topGOdata",
                    description = "EChlor_apo_d10",
                    ontology = "MF",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# run the Fisher's exact tests with the weight01 algorithm
# resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")
down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# summarize the results
up_summary <- GenTable(up_GO_data,
                       weight01 = up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 50)
down_summary <- GenTable(down_GO_data,
                         weight01 = down_result,
                         orderBy = "down_result",
                         ranksOf = "down_result",
                         topNodes = 50)

up_GO_filtered <- up_summary %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

up_GO_filtered %>% 
  ggplot(aes(x = Term, y = GeneRatio)) +
  geom_col(width = 0.05) + 
  geom_point(size = 3) +
  coord_flip() # Flip the axes so the x-axis labels are readable

down_GO_filtered <- down_summary %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

down_GO_filtered %>% 
  ggplot(aes(x = Term, y = GeneRatio)) +
  geom_col(width = 0.05) + 
  geom_point(size = 3) +
  coord_flip() # Flip the axes so the x-axis labels are readable

down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
order_term <- down_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector


down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, colour = weight01)) +
  geom_col(width = 0.05) +
  geom_point(size = 3) +
  coord_flip() +
  scale_x_discrete(limits = order_term) + 
  scale_colour_gradient(low = "red", high = "blue")

down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  
  # Add these to make our plot prettier
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # this changes the scale of the axes
