# Knight Lab Rotation Winter 2019
## Gibraan Rahman

The purpose of this project is to identify relationships between immune cell abundances and newly-discovered microbial communities in TCGA cancers.

### Features:

* Perform Similarity Network Fusion and Consensus Clustering on OTUs/immune cells
* Test cluster results for associations with race, sex, age
* Perform survival analysis on resulting clusters comparing overall survival
* PERMANOVA differential composition analysis of clusters
* (TODO) Differential gene expression using GSEA

### SNF-CC

Create consensus clusters of patients from fused network of microbial profiles and CIBERSORT-inferred immune cell abundances. Performed 9 times with number of clusters = 2:10. Silhouette widths of each of these clusters numbers is saved along with the grouping of each sample.

Input the CIBERSORT file to use and output directory.

Currently trying different pre-processing procedures on the CIBERSORT data to improve clinical relevant of clustered.

Usage:

`qsub -v CIBFILE=<path to cibersort file>,OUTDIR=<output directory> SNF_CC_all.pbs`

### Covariate Testing

Post-hoc analysis of clustering results to ensure that samples are not simply being separated by race, age, or sex. Inputs are cancer type and number of clusters.

TODO: Create submittable script to perform on group file rather than input cancer type and number of clusters.

Usage:

`Rscript covariate_stats.R <cancer_type> <num_clusters>`

### Survival Analysis

Create survival analysis plots for all cancer types (all cluster numbers). This analysis uses updated metadata from TCGA (days to death OR right-censored days to last follow-up) to create survival curves.

TODO: Parameterize this procedure such that a group file can be input instead of a cancer type. This will allow for more flexibility in analysis.

Usage:

`qsub survival_analysis.pbs`

### PERMANOVA

Perform PERMANOVA on cluster results to analyze differential abundances among clusters.

Usage:

`qsub -v GROUP_FILE=<group file> permanova.pbs`

### Feature ANOVA

Perform ANOVA on all features and association with identified group. Right now only supports CIBERSORT features.

TODO: Implement OTU as well.

Usage:

`Rscript feature_anova.R <path to group file> <path to CIBERSORT file> <path to output file>`
