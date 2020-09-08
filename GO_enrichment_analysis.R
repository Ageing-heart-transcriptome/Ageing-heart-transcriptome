library(clusterProfiler)
library(stringr)

setwd("~/UROP 2019/R")
source("get_top_genes.R")

# clean results of GO enrichment
clean_enrichment_results <- function(enrichment_results){
  enrichment_results %>% 
    as.data.frame %>% 
    dplyr::select(Ontology = ONTOLOGY, Description, p.adjust, geneID) %>% 
    # get gene symbols by splitting entrez IDs, converting to symbols, and pasting back together
    mutate(Symbol = str_split(geneID, "/") %>% lapply(Entrez_to_symbol) %>% lapply(paste, collapse = ", ") %>% unlist) %>%
    dplyr::select(Ontology, Description, p.adjust, Symbol) 
}


# get top n terms for each ontology type
get_top_terms <- function(GO_df, ontology_type, n=5){
  rows <- GO_df[which(GO_df$Ontology == ontology_type),]
  ordered <- arrange(rows, p.adjust)
  n <- ifelse(nrow(rows) <= n, nrow(rows), n)
  return(ordered[1:n,]) 
}

# for each time point, conduct gene ontology enrichment analyses, separating genes that are upregulated and downregulated between each group
GO_t1t2_up <- enrichGO(gene = Entrez_IDs_t1t2_up, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t1t2_up_short <- lapply(levels(GO_t1t2_up$Ontology), get_top_terms, GO_df = GO_t1t2_up) %>% do.call(what=rbind)
GO_t1t2_down <- enrichGO(gene = Entrez_IDs_t1t2_down, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t1t2_down_short <- lapply(levels(GO_t1t2_down$Ontology), get_top_terms, GO_df = GO_t1t2_down) %>% do.call(what=rbind)

GO_t2t3_up <- enrichGO(Entrez_IDs_t2t3_up, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t2t3_up_short <- lapply(levels(GO_t2t3_up$Ontology), get_top_terms, GO_df = GO_t2t3_up) %>% do.call(what=rbind)
GO_t2t3_down <- enrichGO(Entrez_IDs_t2t3_down, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t2t3_down_short <- lapply(levels(GO_t2t3_down$Ontology), get_top_terms, GO_df = GO_t2t3_down) %>% do.call(what=rbind)

GO_t3t4_up <- enrichGO(Entrez_IDs_t3t4_up, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t3t4_up_short <- lapply(levels(GO_t3t4_up$Ontology), get_top_terms, GO_df = GO_t3t4_up) %>% do.call(what=rbind)
GO_t3t4_down <- enrichGO(Entrez_IDs_t3t4_down, OrgDb = "org.Mm.eg.db", ont = "ALL")@result %>% clean_enrichment_results
GO_t3t4_down_short <- lapply(levels(GO_t3t4_down$Ontology), get_top_terms, GO_df = GO_t3t4_down) %>% do.call(what=rbind)

# round p values
# only done for GO tables corresponding to downregulated genes, as those corresponding to upregulated genes are empty
GO_t1t2_down_short$p.adjust <- round(GO_t1t2_down_short$p.adjust, digits = 4)
GO_t2t3_down_short$p.adjust <- round(GO_t2t3_down_short$p.adjust, digits = 4)
GO_t3t4_down_short$p.adjust <- round(GO_t3t4_down_short$p.adjust, digits = 4)

# save full tables as CSV files
write.csv(GO_t1t2_down, 'supplementary_files/GO_t1t2_down.csv')
write.csv(GO_t2t3_down, 'supplementary_files/GO_t2t3_down.csv')
write.csv(GO_t3t4_down, 'supplementary_files/GO_t3t4_down.csv')