# Toplogical Data Analysis in Population Genomics

This repo contains code relating to Harry Emeric's
rotation in the Bustamante lab in Spring 2020. The files contain
exploratory work applying TDA to the Himalayan data. The project
was supervised and coordinated by Alex Ioannidis.


## Files

### notebooks/

This directory contains various notebooks used for exploratory analysis and results for this project.

### results/

This directory stores any results used such as flat files and images.

### cocycleIndividualPlot.py

Module with code to create a plot comparing PCA2, Birth - Death plot from the rips complex, and 
for each of these point on the BD plot, called representative cocycyles, how they break down by
population of the samples.

Example usage:

1. From .vcf file:

```python
cocycles_ind_plot_vcf = cocycleIndividualPlot(vcf_file='/home/projects/HimalGenAsia/HimalGen.phase.vcf.gz',
                                              popinfo_path='~/../projects/HimalGenAsia/HimalGen.popinfo.csv')
                                              
fig_vcf = cocycles_ind_plot_vcf.display_cocycle_charts(
    cocycle_number_list=[0,1,2,3,5],
    cocycle_individuals_file='results/cocycle_individuals.txt',
    birth_death_coordinates_file='results/birth_death_coordinates_file.txt')

fig_vcf.suptitle('Population Breakdown of Most Persistant Cocycles for All Principal Components from vcf', fontsize=13)

fig_vcf.show()                                
```

2. From precalculated genotype matrix and ripser object

```python
cocycles_ind_plot_gt_pcs = cocycleIndividualPlot(popinfo_path='~/../projects/HimalGenAsia/HimalGen.popinfo.csv',
                                                 gt_matrix_PCs=gt_matrix_PCs,
                                                 ripser_result=result_gt_pcs)
fig = cocycles_ind_plot_gt_pcs.display_cocycle_charts(cocycle_number_list=[0,1,2,3,5])
fig.suptitle('Population Breakdown of Most Persistant Cocycles for All Principal Components', fontsize=13)

fig.show()

```

### ripser_individuals.ipynb

Contains example usage of `cocycleIndividualPlot` on the Himalayan dataset.

### TDA_in_Population_Genetics.pdf

The final report

### plotly_test.py

Code to test deployment to turn cocycleIndividualPlot.py into a live web app.