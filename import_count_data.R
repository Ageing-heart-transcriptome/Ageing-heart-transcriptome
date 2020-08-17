library(dplyr)

# import raw data set of unnormalised counts
gene_data <- read.csv("data/gene_count_matrix.csv")

# names of columns corresponding to transcript counts
count_columns <- colnames(gene_data)[grep('mouse', colnames(gene_data))]

# set transcript count columns to numeric type
gene_data[,count_columns] <- apply(gene_data[,count_columns], MARGIN = 2, function(col){
  col[is.na(col)] <- 0
  return(as.numeric(col))
  }) %>% as.data.frame

# remove transcripts (rows) with zero counts across all mice
non_zero <- apply(gene_data[,count_columns], MARGIN = 1, function(x){any(x > 0)}) 
gene_data <- gene_data[which(non_zero),]

# masks for annotated and novel transcripts
is_annotated <- -grep("STRG.", gene_data$gene_id)
is_novel <- grep("STRG.", gene_data$gene_id)

# make new dfs for annotated and novel transcripts
annotated <- gene_data[is_annotated,]
novel <- gene_data[is_novel,]

# plot hierarchical clustering of logged transcript counts
log_counts <- log2(annotated[,count_columns] + 1)
scaled_log_counts <- scale(log_counts)
dist_matrix <- dist(t(scaled_log_counts))
clustering <- hclust(dist_matrix, method = 'average')
plot(clustering)

# remove outlier
gene_data <- gene_data[,-grep("mouse11_1", colnames(gene_data))]
annotated <- annotated[,-grep("mouse11_1", colnames(annotated))]
novel <- novel[,-grep("mouse11_1", colnames(novel))]

# make count matrix by averaging technical replicates
combine_technical_replicates <- function(df){
  
  # apply summation over each pair of columns
  combined <- sapply(colnames(df[,seq(1, ncol(df), by = 2)]), function(x){
    # if a column corresponds to 1 of 2 technical replicate, take the sum with the next column
    if(grepl("_1_", x)){
      idx <- match(x, colnames(df))
      sum <- df[,idx]+df[,idx+1]
      return(sum)
    # otherwise just take the column, which corresponds to either a single tech rep or a non-count column
    }else{
      return(df[,x])}
  }) %>% data.frame
  
  # edit column names
  colnames(combined) <- lapply(colnames(df[,seq(1, ncol(df), by = 2)]), function(x){
    if(grepl("wk|mo", x)){
      split <- str_split(x, "_")[[1]][-2]
      newname <- paste(split, collapse = "_")
      return(newname)
    }else{
      return(x)}
  }) 
  
  # subset the columns which correspond to transcript counts
  combined <- combined[,grep("wk|mo", colnames(combined))]
  
  return(combined)
  
}

count_matrix_all <- combine_technical_replicates(gene_data) %>% as.matrix
count_matrix_anno <- combine_technical_replicates(annotated) %>% as.matrix
count_matrix_novel <- combine_technical_replicates(novel) %>% as.matrix


# assign row names to be gene names
rownames(count_matrix_all) <- gene_data[,1]
rownames(count_matrix_anno) <- annotated[,1]
rownames(count_matrix_novel) <- novel[,1]


# reorder columns to be in chronological order
reorder_timepoints <- function(count_matrix){
  
  t1 <- which(grepl("4wk", colnames(count_matrix)))
  t2 <- which(grepl("15wk", colnames(count_matrix)))
  t3 <- which(grepl("8mo", colnames(count_matrix)))
  t4 <- which(grepl("22mo", colnames(count_matrix)))
  
  time_points <- c(t1, t2, t3, t4)
  
  return(count_matrix[,time_points])
  
}

# set new count matrix objects with reordered timepoints
count_matrix_all <- reorder_timepoints(count_matrix_all)
count_matrix_anno <- reorder_timepoints(count_matrix_anno)
count_matrix_novel <- reorder_timepoints(count_matrix_novel)

# make metadata frame of experimental conditions, including both age and sex
conditions <- data.frame(group = c(rep('t1', 3), rep('t2', 4), rep('t3', 4), rep('t4', 3)))
conditions$group <- factor(conditions$group)
conditions$sex <- ifelse(grepl("_F_", colnames(count_matrix_anno)), "F", "M") %>% as.factor
rownames(conditions) <- colnames(count_matrix_anno)
