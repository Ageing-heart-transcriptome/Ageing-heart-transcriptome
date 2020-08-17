library(dplyr)
library(org.Mm.eg.db)
library(stringr)

## This script runs DESeq analysis scripts and processes the results ##
## It also initializes the variables used by scripts for downstream analyses ##

source("DE_analysis_LRT.R")
source("DE_analysis_pairwise.R")

## ID mapping functions ##
Mouse_IDs <- org.Mm.eg.db

# map entrez gene IDs to gene symbols
Entrez_to_symbol <- function(gene_names){
  symbols <- mapIds(Mouse_IDs, keys = gene_names, column = "SYMBOL", keytype = "ENTREZID") 
  return(symbols[!is.na(symbols)])
}

# map gene symbols to entrez gene IDs
Symbol_to_entrez <- function(gene_names){
  ids <- mapIds(Mouse_IDs, keys = gene_names, column = "ENTREZID", keytype = "SYMBOL") 
  return(ids[!is.na(ids)])
}

# map gene symbols to gene names
Symbol_to_name <- function(gene_names){
  ids <- mapIds(Mouse_IDs, keys = gene_names, column = "GENENAME", keytype = "SYMBOL") 
  return(ids[!is.na(ids)])
}


# function for getting transcript counts for genes below a p value threshold
get_counts_for_sig_genes <- function(results, pvalue_threshold=0.05){
  results_df <- as.data.frame(results, row.names=rownames(results))
  significant_subset <- subset(results_df, padj <= pvalue_threshold)
  significant_genes <- rownames(significant_subset)
  return(count_matrix_anno[significant_genes,])
}

## Isolate all genes under padj = 0.05 ##

# LRTs
top_genes_LRT_age <- get_counts_for_sig_genes(results_LRT_age) %>% as.matrix

top_genes_LRT_int <- get_counts_for_sig_genes(results_LRT_int) %>% as.matrix

# Wald tests (pairwise comparisons)
top_genes_t1t2 <- get_counts_for_sig_genes(t1_t2) %>% as.matrix

top_genes_t2t3 <- get_counts_for_sig_genes(t2_t3) %>% as.matrix

top_genes_t3t4 <- get_counts_for_sig_genes(t3_t4) %>% as.matrix

top_genes_t1t3 <- get_counts_for_sig_genes(t1_t3) %>% as.matrix

top_genes_t1t4 <- get_counts_for_sig_genes(t1_t4) %>% as.matrix

top_genes_t2t4 <- get_counts_for_sig_genes(t2_t4) %>% as.matrix

# get up and downregulated genes - returns subset of DESeq results data frame, including on the top n up or downregulated genes
get_top_n_genes <- function(DESeq_results, direction, n){
  dynamic <- subset(DESeq_results, if(direction == "Up"){
    log2FoldChange > 0
    }
    else if(direction == "Down"){
      log2FoldChange < 0
      }
    else{
      log2FoldChange != 0
      })
  
  # delete NA p-values and remove unannotated genes (Gm/Rik)
  dynamic <- dynamic[-c(which(is.na(dynamic$padj)), grep("Gm|Rik", rownames(dynamic))),]
  
  # order by adjusted p-value 
  ordered <- reorder(rownames(dynamic), dynamic$padj) 
  
  # get top n genes
  top_n <- levels(ordered)[1:n] %>% as.character
  
  # make a new column for top n genes
  output <- cbind(Gene = top_n, DESeq_results[top_n,] %>% apply(2, format, digits = 3, nsmall = 2)) %>% as.data.frame
  
  # reformat pvalue strings
  output$padj <- sapply(output$padj, function(x){
    if(grepl("e", x)){
      split <- str_split(x, "e")
      base <- split[[1]][1] %>% format(digits = 2)
      exp <- split[[1]][2] %>% str_remove("\\+") %>% format(digits = 1)
      pvalue_formatted <- paste(base, "*10^", exp, sep = "")
      return(pvalue_formatted)}
    else{return(x)}
  })
  if(direction == "LRT"){
    output$LFC_4wk_15wk<- t1_t2[top_n,"log2FoldChange"] %>% round(digits = 2) %>% format(nsmall = 2)
    output$LFC_15wk_8mo<- t2_t3[top_n,"log2FoldChange"] %>% round(digits = 2) %>% format(nsmall = 2)
    output$LFC_8mo_22mo<- t3_t4[top_n,"log2FoldChange"] %>% round(digits = 2) %>% format(nsmall = 2)
  }
  return(output)
}

# get top 20 genes for both directions in each pairwise comparison 
t1t2_top_downreg <- get_top_n_genes(t1_t2, "Down", 20)
t1t2_top_upreg<- get_top_n_genes(t1_t2, "Up", 20)
t2t3_top_downreg <- get_top_n_genes(t2_t3, "Down", 20)
t2t3_top_upreg <- get_top_n_genes(t2_t3, "Up", 20)
t3t4_top_downreg <- get_top_n_genes(t3_t4, "Down", 20)
t3t4_top_upreg <- get_top_n_genes(t3_t4, "Up", 20)

# top genes for LRTs
LRT_top_age <- get_top_n_genes(results_LRT_age, "LRT", 20)

# get DE genes for all pairwise comparisons, all directions
DE_t1t2_all <- subset(t1_t2, padj <= 0.05) %>% rownames
DE_t1t2_up <- subset(t1_t2, padj <= 0.05 & log2FoldChange > 0) %>% rownames
DE_t1t2_down <- subset(t1_t2, padj <= 0.05 & log2FoldChange < 0) %>% rownames
DE_t2t3_all <- subset(t2_t3, padj <= 0.05) %>% rownames
DE_t2t3_up <- subset(t2_t3, padj <= 0.05 & log2FoldChange > 0) %>% rownames
DE_t2t3_down <- subset(t2_t3, padj <= 0.05 & log2FoldChange < 0) %>% rownames
DE_t3t4_all <- subset(t3_t4, padj <= 0.05) %>% rownames
DE_t3t4_up <- subset(t3_t4, padj <= 0.05 & log2FoldChange > 0) %>% rownames
DE_t3t4_down <- subset(t3_t4, padj <= 0.05 & log2FoldChange < 0) %>% rownames
total_gene_set <- rownames(count_matrix_anno)

# map gene names to entrez IDs for SPIA
Entrez_IDs_t1t2 <- Symbol_to_entrez(DE_t1t2_all)
Entrez_IDs_t1t2_up <- Symbol_to_entrez(DE_t1t2_up)
Entrez_IDs_t1t2_down <- Symbol_to_entrez(DE_t1t2_down)
Entrez_IDs_t2t3 <- Symbol_to_entrez(DE_t2t3_all)
Entrez_IDs_t2t3_up <- Symbol_to_entrez(DE_t2t3_up)
Entrez_IDs_t2t3_down <- Symbol_to_entrez(DE_t2t3_down)
Entrez_IDs_t3t4 <- Symbol_to_entrez(DE_t3t4_all)
Entrez_IDs_t3t4_up <- Symbol_to_entrez(DE_t3t4_up)
Entrez_IDs_t3t4_down <- Symbol_to_entrez(DE_t3t4_down)
Entrez_IDs_all <- Symbol_to_entrez(total_gene_set)

# get rid of NAs
Entrez_IDs_t1t2 <- Entrez_IDs_t1t2[which(!is.na(Entrez_IDs_t1t2))]
Entrez_IDs_t1t2_up <- Entrez_IDs_t1t2_up[which(!is.na(Entrez_IDs_t1t2_up))]
Entrez_IDs_t1t2_down <- Entrez_IDs_t1t2_down[which(!is.na(Entrez_IDs_t1t2_down))]
Entrez_IDs_t2t3 <- Entrez_IDs_t2t3[which(!is.na(Entrez_IDs_t2t3))]
Entrez_IDs_t2t3_up <- Entrez_IDs_t2t3_up[which(!is.na(Entrez_IDs_t2t3_up))]
Entrez_IDs_t2t3_down <- Entrez_IDs_t2t3_down[which(!is.na(Entrez_IDs_t2t3_down))]
Entrez_IDs_t3t4 <- Entrez_IDs_t3t4[which(!is.na(Entrez_IDs_t3t4))]
Entrez_IDs_t3t4_up <- Entrez_IDs_t3t4_up[which(!is.na(Entrez_IDs_t3t4_up))]
Entrez_IDs_t3t4_down <- Entrez_IDs_t3t4_down[which(!is.na(Entrez_IDs_t3t4_down))]
Entrez_IDs_all <- Entrez_IDs_all[which(!is.na(Entrez_IDs_all))]

# set gene symbols to a new object - these are the subset of DE genes with entrez IDs
t1t2_genes <- names(Entrez_IDs_t1t2)
t1t2_genes_up <- names(Entrez_IDs_t1t2_up)
t1t2_genes_down <- names(Entrez_IDs_t1t2_down)
t2t3_genes <- names(Entrez_IDs_t2t3)
t2t3_genes_up <- names(Entrez_IDs_t2t3_up)
t2t3_genes_down <- names(Entrez_IDs_t2t3_down)
t3t4_genes <- names(Entrez_IDs_t3t4)
t3t4_genes_up <- names(Entrez_IDs_t3t4_up)
t3t4_genes_down <- names(Entrez_IDs_t3t4_down)
all_genes <- names(Entrez_IDs_all)