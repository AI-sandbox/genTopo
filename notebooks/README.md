# Notebooks

This directory contains various notebooks used for exploratory analysis and results for this project.

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
Adapted from original by Brad Nelson.

### visualize_cocylces.ipynb

Contains code to create a plotly graph object with both a small dataset and
one with a sample of 5000 points (approximately the number of PCs in the
Himalayan data) to see what cocycles look like.


## PCA Pipeline in plink2

To get eigenvectors and eigenvalues, run this from the directory you want to keep the .eigenvec and .eigenval files


`plink2 --vcf <filename> --pca`

For example:

phase data was used with filled in missing data
`plink2 --vcf /home/projects/HimalGenAsia/HimalGen.phase.vcf.gz --pca`

`plink2 --vcf /home/projects/HimalGenAsia/HimalGen.Final.Asia.cln.vcf --pca`


Convert .bed .bim .fam to .vcf with 

```bash
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
```

ref: https://www.biostars.org/p/108499/