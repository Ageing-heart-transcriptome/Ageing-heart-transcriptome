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

We performed an exploratory hierarchical clustering and identified an outlier technical replicate in the 4 week group.

#### **`import_count_data.R`**
```
dist_matrix <- dist(t(scaled_log_counts))
clustering <- hclust(dist_matrix, method = 'average')
plot(clustering)
```

The script removes the outlier replicate, combines the other technical replicates, and generates the metadata for the experimental design. 

### Differential expression analysis

We provide two scripts for differential expression analysis: `DE_analysis_LRT.R` and `DE_analysis_pairwise.R`.

`DE_analysis_LRT.R` conducts a likelihood ratio test by comparing the likelihoods of the age-included and age-excluded models. This identifies genes whose expression levels are well-explained by the age variables in the model. In other words, it identifies genes that are broadly associated with the ageing process.

`DE_analysis_pairwise.R` conducts Wald tests for GLM contrast coefficients (see [DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/)). This measures differential expression between pairs of age groups, identifying genes whose expression varies significantly between two different age groups.

### Downstream analysis

The downstream analysis begins with a script processing the results of the differential expression analysis: `get_top_genes.R`. This scripts loads objects containing information about the genes identified as statistically-significant by DESeq2.

The four other scripts perform further analyses on the objects created by `get_top_genes.R`. These are:
- `DE_gene_plots.R` -- creates volcano plots and transcript count plots for all statistically significant genes 
- `gene_clustering.R` -- performs hierarchical clustering on the genes identified as statistically significant by the likelihood ratio test
- `SPIA_pathway_analysis.R` -- uses SPIA to perform pathway enrichment analysis on genes identified as statistically significant by pairwise tests between age groups
- `GO_enrichment_analysis.R` -- conducts gene ontology (GO) enrichment analysis on statistically significant genes identified in pairwise tests between age groups
