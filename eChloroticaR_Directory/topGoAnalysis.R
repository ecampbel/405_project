# Installing topGO
# BiocManager::install("topGO")

# load topGO
library("topGO")

#filter DESeq for significance
sig <- na.omit(res)
sig <- sig[sig$padj<0.05,]
sigdf <- as.data.frame(sig)
sigdf <- tibble::rownames_to_column(sigdf, "gene_id")
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

up_GO_filtered_arranged <- up_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
up_order_term <- up_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

up_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  ggtitle("Significant Upregulated GO Terms with Ontology MF") +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = up_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  
  # Add these to make our plot prettier
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # this changes the scale of the axes

down_GO_filtered <- down_summary %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  filter(weight01 <= 0.05) %>%
  head(n = 20)

down_GO_filtered_arranged <- down_GO_filtered %>% 
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term))

# Now let's extract the order of the term column
down_order_term <- down_GO_filtered_arranged %>% 
  pull(Term) # pull() extracts a column as a vector

down_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  ggtitle("Significant Downregulated GO Terms with Ontology MF") +
  geom_col(width = 0.05) +
  geom_point(aes(size = Significant)) + 
  coord_flip() +
  scale_x_discrete(limits = down_order_term) +
  scale_color_gradient(low = "red", high = "blue") +
  
  # Add these to make our plot prettier
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) # this changes the scale of the axes

# Add labels to upregulated and downregulated dataframes
up_GO <- up_summary %>% 
  mutate(up_down = "UP")

down_GO <- down_summary %>% 
  mutate(up_down = "DOWN")

# Make a joined dataframe
joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
  filter(weight01 <= 0.05) %>%
  mutate(GeneRatio = Significant/Annotated, weight01 = as.numeric(weight01)) %>%
  arrange(GeneRatio) %>%
  mutate(Term = factor(Term)) %>%
  head(n = 40)

# Extract the column order
order_term_joined <- joined_GO_filtered_arranged %>% 
  pull(Term)

joined_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = weight01)) +
  geom_point(aes(size= Significant)) +
  coord_flip() +
  scale_x_discrete(limits = order_term_joined) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") +
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96"), 
        strip.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  facet_grid(.~ up_down)

echlor_GOIDs <- read_tsv("echlor_GOIDs.tsv", col_names = c("gene_id", "GO_id"))
GO_ID <- down_GO_filtered_arranged$GO.ID[2]
filtered_GO_ID <- echlor_GOIDs %>%
  filter(grepl(pattern = GO_ID, GO_id))

gene_id_vector <- filtered_GO_ID %>% 
  pull(gene_id)

# Filter the data set
dat_GO_filtered <- sigdf %>%
  filter(gene_id %in% gene_id_vector) %>%
  filter(padj <= 0.05 & log2FoldChange != 0)

dat_GO_filtered %>% 
  ggplot(aes(x = gene_id, y = log2FoldChange)) +
  ggtitle("Genes Associated with GO:0004104 - Cholinesterase Activity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax =log2FoldChange + lfcSE), width = 0.4) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  coord_flip()+ 
  labs(y = "log2(FC)", x= "Gene ID", fill = "Regulation") +
  geom_hline(yintercept = 0) +
  theme_light() +
  theme(panel.border = element_rect(color = "black"), panel.grid = element_line(colour = "white"))