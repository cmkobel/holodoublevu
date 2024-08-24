# Explanation of plots in holodoublevu/results/*

This file aims to help navigate and understand the output plots in the results/ directory.

Sorted by order of creation in the pipeline.


**results/filtered/filter_f1.pdf**

Scatterplot with number of proteins and the average intensity per protein. On the number of proteins, some samples have been excluded from downstream analysis.

There are three imputation groups: [aberdeen, luing, both]. The following results exists for each of these groups.


## Imputation

**no figs yet.**

## Module calling

**results/ig/{imputation_group}/wgcna/pst_*.pdf**

Before the modules are called, the proteomic intensities are transformed with an exponent. This exponent is also known as a soft threshold and is set such that scale freeness of the protein network is achieved. Digesta layers tend to get a threshold around 12-15. Liver and well 6-8.
  
**results/ig/{imputation_group}/wgcna/hclust_*.pdf**

Hierarchical clustering of _all_ proteins and heatmap showing phenotypic traits
Can be used to see if any traits are obviously linked to a group of samples.


## Inspection of modules
  

**results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_digestaXliver_*.pdf**
**results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_digestaXwall_*.pdf**
**results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_liverXwall_*.pdf**

These files represent correlating modules pairwisely across the digesta-wall-liver axis. Significant correlations are marked.
  
**results/ig/{imputation_group}/wgcna/inspected/dendro_*.pdf**
  
Proteins are hierarchically clustered to form modules. This is a diagnostic plot which has no biological interpretative quality.
  
**results/ig/{imputation_group}/wgcna/inspected/count_*.pdf**

Every protein is allocated to a single module. If it isn't it ends up in "module 0". Module names (0...N) are arranged by the number of proteins.

**results/ig/{imputation_group}/wgcna/inspected/pheno_*.pdf**
  
Modules are correlated (across proteins) to phenotypic traits.

**results/ig/{imputation_group}/wgcna/inspected/me_{imputation_group}_*.pdf**
  
Modules are correlated to themselves. Useful for inspecting merging of modules.


## Module membership and trait significance

**results/ig/{imputation_group}/wgcna/pathway_enrichment/pathway_*.pdf**

The proteins in each module are subjected to pathway enrichment analysis. This heatmap shows which pathways are significantly enriched for in each module. In the bottom, another heatmap shows which modules are significantly correlated to phenotypic traits.
  
**results/ig/{imputation_group}/wgcna/species_table/tax_binomial_*.pdf**

The proteins in each module come from a specific species. Thes heatmap shows the count of proteins from each species in each module. In the bottom, another heatmap shows which modules are significantly correlated to phenotypic traits.


**results/ig/{imputation_group}/wgcna/inspected/module_membership_trait_significance/moduletrait_\<trait>_*_\<module>*.pdf**

The correlation between protein and module (across samples) is also referred to as module membership or module eigengene connectivity. Proteins with high module membership are regarded as hub nodes because they're well connected within the module. These plots show the module-trait correlation as a function of module membership.
