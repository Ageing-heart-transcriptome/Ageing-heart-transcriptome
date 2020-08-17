library(DEGreport)
library(ggplot2)
library(clusterProfiler)
library(factoextra)

setwd("~/UROP 2019/R")
source("get_top_genes.R")

## Gene clustering ##

# pseudocounts for log2 transformation
LRT_age_genes_cluster <- top_genes_LRT_age+1
LRT_age_genes_cluster <- log2(LRT_age_genes_cluster)
sig_genes_LRT <- rownames(results_LRT_age[which(results_LRT_age$padj < 0.05),])

# run clustering
clustering_age <- degPatterns(LRT_age_genes_cluster[sig_genes_LRT,], conditions, time = "group", reduce = TRUE)

# modify gene names so we can grab their p values
gene_names <- ifelse(grepl('^AC|^AL|^CR', clustering_age$df$genes), clustering_age$df$genes, str_replace_all(clustering_age$df$genes, '\\.', '-'))
gene_names <- ifelse(grepl('^X\\d', gene_names), str_replace_all(gene_names, '^X', ''), gene_names)
clustering_age$df$padj <- results_LRT_age$padj[match(gene_names, rownames(results_LRT_age))]

# get genes that were excluded from the clustering
unclustered_genes <- rownames(top_genes_LRT_age)[!(rownames(top_genes_LRT_age) %in% gene_names)]
unclustered_gene_pvalues <- results_LRT_age$padj[match(unclustered_genes, rownames(results_LRT_age))]
unclustered_genes_df <- data.frame(genes=unclustered_genes, cluster=NA, padj=unclustered_gene_pvalues)

# append unclustered genes to the clustering df, to save as supplementary file
full_clustering_df <- rbind(clustering_age$df, unclustered_genes_df)

# write csv
write.csv(full_clustering_df, 'supplementary_files/Supplementary file 4.csv')

# make a plot of the clustering groups
cluster_plot_age <- degPlotCluster(table = clustering_age$normalized, time = "group") + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.position = "none", text = element_text(size = 16), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

ggsave('figures/Figure 3.png', cluster_plot_age)

# make list with member genes listed under each cluster
genes_per_cluster <- lapply(clustering_age$df$cluster, function(cluster){
  clustering_age$df$genes[clustering_age$df$cluster == cluster]
})
names(genes_per_cluster) <- clustering_age$df$cluster


# use hypergeometric test to scan clusters for overrepresentation (compareCluster function used)
Profile_clusters <- function(cluster_df){
  # get genes in each cluster
  genes <- lapply(unique(cluster_df$cluster), function(x){
    return(cluster_df$genes[which(cluster_df$cluster == x)])})
  # get entrez ids from gene ids
  entrez <- lapply(genes, function(x){
    entrez_ids <- Symbol_to_entrez(x)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    return(entrez_ids)
  })
  # add X cluster numbers
  names(entrez) <- paste("X", levels(clustering_age$df$cluster), sep="")
  result_GO <- compareCluster(entrez, fun="enrichGO", OrgDb = "org.Mm.eg.db", qvalueCutoff=0.1)
  result_KEGG <- compareCluster(entrez, fun="enrichKEGG", organism = "mmu", qvalueCutoff=0.1)
  all_results <- list(GO=result_GO, KEGG=result_KEGG)
  return(all_results)
}

# cluster profiling results
age_cluster_results <- Profile_clusters(clustering_age$df)

# get gene names associated with gene ontology for Ub transferase-like activity
Ub_transferase_genes <- age_cluster_results$GO@compareClusterResult$geneID[1] %>% str_split("/") %>% 
  unlist %>% Entrez_to_symbol()

# get gene names associated with KEGG pathway for protein processing in the ER
Protein_processing_genes <- age_cluster_results$KEGG@compareClusterResult$geneID[1] %>% str_split("/") %>% 
  unlist %>% Entrez_to_symbol()

# compress df into mean count values for each time point 
# can be used to double-check clustering done by degPatterns
ages <- c("4wk", "15wk", "8mo", "22mo")
compressed_df <- lapply(ages, function(age){
  df <- annotated[,grep(age, colnames(annotated))]
  compressed <- apply(df, 1, mean) %>% as.data.frame
  rownames(compressed) <- rownames(df)
  colnames(compressed) <- age
  return(compressed)
}) %>% do.call(what=cbind)
rownames(compressed_df) <- annotated$gene_id
