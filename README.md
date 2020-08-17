# Cross-sectional transcriptional analysis of the murine cardiac lifespan

Code for the transcriptomic analysis [published by Greenig et al](https://www.frontiersin.org/articles/10.3389/fmolb.2020.565530) in Frontiers in Molecular Biosciences. We analysed RNA-sequencing data collected from cardiac tissue of mice of four different age groups, identifying components and processes associated with ageing in wild type mice.

## Table of contents

The workflow can be roughly divided into three parts:
1. Preprocessing
2. Differential expression analysis
3. Downstream analysis

### Preprocessing

We provide a single script - `import_count_data.R` to import and clean the transcript count data. The data is sourced from the data/ directory, as follows:

#### **`import_count_data.R`**
```
gene_data <- read.csv('data/gene_count_matrix.csv')
```

In this script we perform an exploratory hierarchical clustering and identified an outlier technical replicate in the 4 week group.

#### **`import_count_data.R`**
```
log_counts <- log2(annotated[,count_columns] + 1)
scaled_log_counts <- scale(log_counts)
dist_matrix <- dist(t(scaled_log_counts))
clustering <- hclust(dist_matrix, method = 'average')
plot(clustering)
```

The script removes the outlier replicate, combines the other technical replicates, and generates the metadata for the experimental design. 

### Differential expression analysis

We provide two scripts for differential expression analysis: `DE_analysis_LRT.R` and `DE_analysis_pairwise.R`.

The base model contains a separate coefficient for each age group, a single coefficient for sex, and coefficients for age/sex combinations (interaction terms).

```
dds <- DESeqDataSetFromMatrix(countData = deseq_matrix, colData = conditions, design=~group*sex)
dds <- estimateSizeFactors(dds)
```

`DE_analysis_LRT.R` conducts a likelihood ratio test by comparing the likelihood of the age-included model to the likelihood of the model in which all age coefficients are removed:

#### **`DE_analysis_LRT.R`**
``` 
dds_LRT_age <- DESeq(dds, test="LRT", reduced=~sex)
```

This identifies genes whose expression levels are well-explained by the age variables in the model. In other words, it identifies genes that are broadly associated with the ageing process.

We also test a reduced model in which only the age/sex interaction coefficients were removed:

#### **`DE_analysis_LRT.R`**
``` 
dds_LRT_int <- DESeq(dds, test="LRT", reduced=~group+sex)
```
This identifies genes whose expression are well-explained by the age/sex interaction variables in the model. These are genes whose expression varies based on an interaction between sex and age - they are male-expressed at certain ages, and female-expressed at others.


`DE_analysis_pairwise.R` conducts Wald tests for GLM contrast coefficients (see [DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/)). This measures pairwise differential expression, identifying genes whose expression varies significantly between two different age groups.

For example, to run a Wald test on the contrast coefficient calculated between the 4 week and the 15 week coefficients, we use:

#### **`DE_analysis_pairwise.R`**
``` 
t1_t2 <- results(dds, contrast=c("group", "t2", "t1"))
```

whereas for the 15 week-8 month comparison we use:

#### **`DE_analysis_pairwise.R`**
``` 
t2_t3 <- results(dds, contrast=c("group", "t3", "t2"))
```

### Downstream analysis

The downstream analysis begins with a script processing the results of the differential expression analysis: `get_top_genes.R`. This scripts loads objects containing information about the genes identified as statistically-significant by DESeq2.

The four other scripts perform further analyses on the objects created by `get_top_genes.R`. These are:
- `DE_gene_plots.R` -- creates volcano plots and transcript count plots for all statistically significant genes 
- `gene_clustering.R` -- performs hierarchical clustering on the genes identified as statistically significant by the likelihood ratio test
- `SPIA_pathway_analysis.R` -- uses SPIA to perform pathway enrichment analysis on genes identified as statistically significant by pairwise tests between age groups
- `GO_enrichment_analysis.R` -- conducts gene ontology (GO) enrichment analysis on statistically significant genes identified in pairwise tests between age groups
