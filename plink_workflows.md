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