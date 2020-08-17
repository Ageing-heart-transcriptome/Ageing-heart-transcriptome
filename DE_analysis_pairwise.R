library(dplyr)
library(stringr)
library(DESeq2)

# run import data script to import count data and experiment conditions
source("import_count_data.R")

# initialise DESeq matrix using count matrix for annotated genes only
deseq_matrix <- count_matrix_anno

## DESEq code ##

# construct DESeq model with interaction term included
dds <- DESeqDataSetFromMatrix(countData = deseq_matrix, colData = conditions, design=~group)
dds <- estimateSizeFactors(dds)
normalised_counts <- counts(dds, normalized=TRUE)

dds <- DESeq(dds)

t1_t2 <- results(dds, contrast=c("group", "t2", "t1"))
rownames(t1_t2) <- rownames(count_matrix_anno)

t1_t2_significant <- subset(t1_t2, padj < 0.05) %>% as.data.frame %>% dplyr::select(log2FoldChange, padj)
write.csv(t1_t2_significant, 'supplementary_files/Supplementary file 1.csv')

t2_t3 <- results(dds, contrast=c("group", "t3", "t2"))
rownames(t2_t3) <- rownames(count_matrix_anno)

t2_t3_significant <- subset(t2_t3, padj < 0.05) %>% as.data.frame %>% dplyr::select(log2FoldChange, padj)
write.csv(t2_t3_significant, 'supplementary_files/Supplementary file 2.csv')

t3_t4 <- results(dds, contrast=c("group", "t4", "t3"))
rownames(t3_t4) <- rownames(count_matrix_anno)

t3_t4_significant <- subset(t3_t4, padj < 0.05) %>% as.data.frame %>% dplyr::select(log2FoldChange, padj)
write.csv(t3_t4_significant, 'supplementary_files/Supplementary file 3.csv')

t1_t3 <- results(dds, contrast=c("group", "t3", "t1"))
rownames(t1_t3) <- rownames(count_matrix_anno)

t1_t4 <- results(dds, contrast=c("group", "t4", "t1"))
rownames(t1_t4) <- rownames(count_matrix_anno)

t2_t4 <- results(dds, contrast=c("group", "t4", "t2"))
rownames(t2_t4) <- rownames(count_matrix_anno)