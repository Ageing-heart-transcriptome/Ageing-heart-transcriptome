# Ageing-heart_transcriptome

Code for the transcriptomic analysis [published](https://www.frontiersin.org/articles/10.3389/fmolb.2020.565530) by Greenig et al in Frontiers in Molecular Biosciences. We analysed RNA-sequencing data collected from cardiac tissue of mice of four different age groups, identifying components and processs associated with the ageing process in wild type mice.

## Table of contents

The workflow can be roughly divided into three parts:
- Preprocessing
- Differential expression analysis
- Downstream analysis

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

`DE_analysis_pairwise.R` conducts Wald tests for GLM contrast coefficients ([see DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/)). This measures differential expression between pairs of age groups, identifying genes whose expression varies significantly between two different age groups.

### Downstream analysis

