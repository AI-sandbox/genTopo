# Toplogical Data Analysis in Population Genomics

This repo contains code relating to Harry Emeric's
rotation in the Bustamante lab in Spring 2020. The files contain
exploratory work applying TDA to the Himalayan data. The project
was supervised and coordinated by Alex Ioannidis.

## Files

### ripser.ipynb

Explores the python ripser package and charts birth-death plots 
of the Himalayan data. Due to compute, dimensionality reduction is
key, and various subsets of the PCs as well as precomuted distance
matrices on various metrics are considered.

### Plot_PCA2.ipynb

This notebook explores some PCA2 plots for various datasets along 
with a pipeline for getting PCAs using plink2.

### plink_workflows.md

Contains linux commands for converting between .bed .bim .fam and .vcf,
and how to run pca.

### distance_matrix.ipynb

Has some experiments realting to distance matrix and inclusion rips.
Adapted froa original by Brad Nelson.
