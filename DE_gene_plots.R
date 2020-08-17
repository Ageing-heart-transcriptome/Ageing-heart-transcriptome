library(EnhancedVolcano)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)

setwd("~/UROP 2019/R")
source("get_top_genes.R")

## Volcano plots (figure 1) ##

LRT_plot <- EnhancedVolcano(results_LRT_age, lab = rownames(results_LRT_age), 
                              x = "log2FoldChange", y="padj",
                              title = "Likelihood ratio test", transcriptLabSize = 6, 
                              subtitle = "Significance of age's effect on gene \nexpression",
                              caption = "x-axis indicates LFC between 4 week and 22 month groups",
                              captionLabSize = 14, titleLabSize = 24,
                              subtitleLabSize = 22,
                              pCutoff = 0.001) + xlim(-40, 40)
t1_t2_plot <- EnhancedVolcano(t1_t2, lab = rownames(t1_t2), 
                                x = "log2FoldChange", y="padj",
                                title = "Wald test", transcriptLabSize = 6, 
                                subtitle = "4 weeks | 15 weeks", 
                                caption = "", titleLabSize = 24,
                                subtitleLabSize = 22,
                                pCutoff = 0.001) + xlim(-40, 40)
t2_t3_plot <- EnhancedVolcano(t2_t3, lab = rownames(t2_t3), 
                                x = "log2FoldChange", y="padj",
                                title = "Wald test", transcriptLabSize = 6, 
                                subtitle = "15 weeks | 8 months",
                                caption = "", titleLabSize = 24,
                                subtitleLabSize = 22,
                                pCutoff = 0.001) + xlim(-50, 50)
is_akt1 <- rownames(t3_t4) == 'Akt1'
t3_t4 <- rbind(t3_t4[is_akt1,], t3_t4[!is_akt1,])
t3_t4_plot <- EnhancedVolcano(t3_t4, lab = rownames(t3_t4), 
                                x = "log2FoldChange", y="padj",
                                title = "Wald test", transcriptLabSize = 6, 
                                subtitle = "8 months | 22 months",
                                caption = "", titleLabSize = 24,
                                subtitleLabSize = 22,
                                pCutoff = 0.001) + xlim(-40, 40)
 
grid.arrange(LRT_plot, t1_t2_plot, t2_t3_plot, t3_t4_plot)

ggsave("figures/Figure 1.png")

## Heatmap (not included in manuscript) ##

# apply DESeq's rlog transform to normalise counts
rlog_dds_anno <- rlog(count_matrix_anno)
rownames(rlog_dds_anno) <- rownames(count_matrix_anno)

# get top annotated genes
top_genes <- get_top_n_genes(rlog_dds_anno, results_LRT_age[-grep("Gm|Rik", rownames(results_LRT_age)),], "padj", 40)

# shorten column names for heatmap
indexes <- sapply(gregexpr("_", colnames(top_genes)), "[", 1)
colnames(top_genes) <- substr(colnames(top_genes), start=indexes+1, stop=sapply(colnames(top_genes), nchar))

# distance matrices
sampledistmatrix <- rlog_dds_anno %>% t %>% dist 

# colours for plot
reds <- colorRampPalette(brewer.pal(3, "Reds"))(255)

sample_heatmap <- heatmap.2(as.matrix(sampledistmatrix), col=reds, trace = "none", key.title = NA)
gene_heatmap <- heatmap.2(as.matrix(top_genes), col=reds, trace = "none", key.title = NA, key.xlab = "Normalised expression", 
                          main = "Top 40 age-dependent genes")

## Plotting transcript counts over time (figure 2) ##
 
# make plotting function
do_ggplot <- function(df){
  gene <- as.character(df$gene)
  gene_name <- tryCatch(expr = {Symbol_to_name(gene)}, error = function(e){gene}, warning = function(w){gene}) %>% str_to_sentence
  plot <- ggplot(df, aes(x = group, y = count, color = sex, group = interaction(gene, sex))) + geom_point() + 
    stat_summary(fun.y = mean, geom = "line") + scale_y_log10() + labs(y = "Count", color = "Sex", title = gene_name) + 
    theme(plot.title = element_text(hjust = 0.5, size=14), text=element_text(size=14))
  return(plot + scale_x_discrete(name = "Age group", breaks = c("t1", "t2", "t3", "t4"), labels = c("4 wk", "15 wk", "8 mo", "22 mo")))
}

# plot counts over time, arrange in a grid for top n genes
make_plots <- function(DESeq_dataset, DESeq_results, group, title, top_n=12){
  
  # get rid of unannotated genes
  unannotated <- grep("Gm|Rik", rownames(DESeq_results))
  results <- DESeq_results[-unannotated,]
  dataset <- DESeq_dataset[-unannotated,]
  
  # get row indices of top n genes in terms of padj
  top_pval_idx <- order(results$padj)[1:top_n]
  
  # apply plotCounts function to the gene corresponding to every index of that list
  plot_data_int <- lapply(top_pval_idx, plotCounts, dds = DESeq_dataset, intgroup = c("group", "sex"), returnData = T)
  
  # add gene names to each df
  plot_data_dfs <- mapply(function(df, gene){
    cbind(df, gene = rep(rownames(results)[gene], nrow(df)))
    }, plot_data_int, top_pval_idx) %>% 
    apply(2, as.data.frame)
  
  # apply plotting function to our list of dfs
  ggplots <- lapply(plot_data_dfs, do_ggplot)
  
  # make grob object with title
  title <- bquote(.(title))
  
  #Make text element for big title in plot grid
  big_title <- textGrob(title, gp=gpar(fontsize=28, fontface="bold"))
  
  #Put grid together by applying arrange to our list of plots
  grid.arrange(grobs = ggplots, top = big_title)
}

# figure 2 (age/sex interaction plots)
int_plots <- make_plots(dds_LRT_int, results_LRT_int,
                        title = "Genes with expression patterns associated with age/sex interaction\n")

ggsave("figures/Figure 2.png")


# count plots for LRT for age top 12 genes, not used in manuscript
age_plots <- make_plots(dds_LRT_age, results_LRT_age, 
                        title = "Genes with significant differential expression associated with age\n")
