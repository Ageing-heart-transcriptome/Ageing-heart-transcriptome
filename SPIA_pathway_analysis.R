library(pathview)
library(SPIA)
library(graphite)

source("get_top_genes.R")

## Pathway analysis with SPIA ## 

# get KEGG pathways
pwys <- pathways(species="mmusculus", database="kegg")
prepareSPIA(pwys, "pwys_for_spia")

# isolate DESeq2's log2 fold changes
LFCs_t1t2 <- t1_t2[t1t2_genes,]$log2FoldChange
names(LFCs_t1t2) <- Entrez_IDs_t1t2

LFCs_t2t3 <- t2_t3[t2t3_genes,]$log2FoldChange
names(LFCs_t2t3) <- Entrez_IDs_t2t3

LFCs_t3t4 <- t3_t4[t3t4_genes,]$log2FoldChange  
names(LFCs_t3t4) <- Entrez_IDs_t3t4

# get pathway regulation score (PRS)
Pathway_results_t1t2 <- spia(de=LFCs_t1t2, all=Entrez_IDs_all, organism="mmu", combine="fisher")
Pathway_results_t2t3 <- spia(de=LFCs_t2t3, all=Entrez_IDs_all, organism="mmu", combine="fisher")
Pathway_results_t3t4 <- spia(de=LFCs_t3t4, all=Entrez_IDs_all, organism="mmu", combine="fisher")

# isolate pathways with p < 0.05
Get_DE_pathways <- function(pathway_results){
  DE_pathways_idx <- which(pathway_results[,"pG"] < 0.05)
  return(pathway_results[DE_pathways_idx,"Name"])
}

# get DE pathways for pairwise comparisons
DE_pathways_t1t2 <- Get_DE_pathways(Pathway_results_t1t2)
DE_pathways_t2t3 <- Get_DE_pathways(Pathway_results_t2t3)
DE_pathways_t3t4 <- Get_DE_pathways(Pathway_results_t3t4)

# find genes corresponding to different pathways
Pathways_to_genes <- function(Pathway_results){
  DE_pathways <- which(Pathway_results$pG < 0.05)
  genes <- Pathway_results[DE_pathways,"KEGGLINK"]
  genes_remove_url <- str_split(genes, "\\?") %>% lapply(function(x){x[-1]})
  genes_sep <- str_split(genes_remove_url, "\\+") %>% lapply("[", -1)
  names(genes_sep) <- Pathway_results[DE_pathways,"Name"]
  return(genes_sep)
}

# isolate genes for all pairwise comparisons
Pathway_genes_t1t2 <- Pathways_to_genes(Pathway_results_t1t2)
Pathway_genes_t2t3 <- Pathways_to_genes(Pathway_results_t2t3)
Pathway_genes_t3t4 <- Pathways_to_genes(Pathway_results_t3t4)


## Pathway visualisation with pathview ##

# function to get KEGG ID from pathway name
Names_to_IDs <- function(pathway_names){
  IDs <- sapply(pathway_names, function(x){
    pathway_obj <- pwys[[x]]
    if(!is.null(pathway_obj)){
      full_id <- pathwayId(pathway_obj)
      id_no_species <- unlist(str_split(full_id, ":"))[2] 
      return(id_no_species)
    }else{
      return(NA)
    }
  })
  return(IDs[!is.na(IDs)])
}

# function for running pathview to visualize the pathway data from SPIA
Run_pathview <- function(DESeq_results, pathway, file_name = NULL){
  
  # isolate all non-NA value LFCs 
  LFCs <- DESeq_results[,"log2FoldChange"]
  names(LFCs) <- rownames(DESeq_results)
  LFCs <- LFCs[which(!is.na(LFCs))]
  
  # get genes that translate properly to entrez (no NAs)
  gene_ids <- Symbol_to_entrez(names(LFCs))
  LFCs_final <- LFCs[names(gene_ids)]
  names(LFCs_final) <- sapply(names(LFCs_final), function(name){gene_ids[name]})
  
  # get info from pathway name
  pathway_id <- Names_to_IDs(pathway)
  if(is.null(file_name)){
    file_name <- str_replace_all(pathway, ' ', '_') 
  }
  
  setwd('supplementary_files')
  
  # use pathview to visualize pathways and output using file_names
  pathway_map <- pathview(LFCs_final, pathway.id = pathway_id, species = "mmu", out.suffix = file_name, 
                          high = "yellow", low = "red", limit = c(-4,4), node.sum="sum")
  
  setwd('..')
  return(pathway_map)
}

# pathview for pathways identified by spia
spia_pathview_t1t2 <- Run_pathview(t1_t2, "Regulation of actin cytoskeleton", file_name = 'Figure 4')
spia_pathview_t1t2 <- Run_pathview(t1_t2, "Wnt signaling pathway")
spia_pathview_t2t3 <- Run_pathview(t2_t3, "VEGF signaling pathway")
spia_pathview_t3t4 <- Run_pathview(t3_t4, "T cell receptor signaling pathway")