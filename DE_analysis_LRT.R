library(dplyr)
library(stringr)
library(DESeq2)

# run import data script to import count data and experiment conditions
source("import_count_data.R")

# initialise DESeq matrix using count matrix for annotated genes only
deseq_matrix <- count_matrix_anno

## DESEq analysis ###

# construct DESeq model with interaction term included
dds <- DESeqDataSetFromMatrix(countData = deseq_matrix, colData = conditions, design=~group*sex)
dds <- estimateSizeFactors(dds)
normalised_counts_full <- counts(dds, normalized=TRUE)

# likelihood ratio test for age
dds_LRT_age <- DESeq(dds, test="LRT", reduced=~sex)
results_LRT_age <- results(dds_LRT_age)
rownames(results_LRT_age) <- rownames(count_matrix_anno)

# likelihood ratio test for age/sex interaction
dds_LRT_int <- DESeq(dds, test="LRT", reduced=~group+sex)
results_LRT_int <- results(dds_LRT_int, contrast = c("sex", "M", "F"))
rownames(results_LRT_int) <- rownames(count_matrix_anno)

# save LRT for age/sex interaction to supplementary figures (LRT for age results saved after clustering)
LRT_int_significant <- subset(results_LRT_int, padj < 0.05) %>% as.data.frame %>% dplyr::select(log2FoldChange, padj)
write.csv(t1_t2_significant, 'supplementary_files/Supplementary file 5.csv')