# Explanation of plots in holodoublevu/results/*

This file aims to help navigate and understand the output plots in the results/ directory.

Sorted by order of creation in the pipeline.


  - results/filtered/filter_f1.pdf
    Scatterplot with number of proteins and the average intensity per protein. On the number of proteins, some samples have been excluded from downstream analysis.
    
There are three imputation groups: [aberdeen, luing, both]. The following results exists for each of these groups.


## Imputation

  - no figs yet.

## Module calling

  - results/ig/{imputation_group}/wgcna/pst_fig5.pdf
    Before the modules are called, the tra
  
  - results/ig/{imputation_group}/wgcna/hclust_*.pdf
    Hierarchical clustering of _all_ proteins and heatmap showing phenotypic traits
    Can be used to see if any traits are obviously linked to a group of samples.
    
    
## Inspection of modules
  
    
  - results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_digestaXliver_fig14.pdf
  - results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_digestaXwall_fig13.pdf
  - results/ig/{imputation_group}/wgcna/inspected/axis_{imputation_group}_liverXwall_fig15.pdf
    These files represent correlating modules pairwisely across the digesta-wall-liver axis. Significant correlations are marked.
    
  - results/ig/{imputation_group}/wgcna/inspected/count_*.pdf
    Every protein is allocated to a single module. If it isn't it ends up in "module 0". Module names (0...N) are arranged by the number of proteins.
    
  - results/ig/{imputation_group}/wgcna/inspected/dendro_*.pdf

  - results/ig/{imputation_group}/wgcna/inspected/pheno_*.pdf
    
    
  - results/ig/{imputation_group}/wgcna/inspected/me_{imputation_group}_g1_D_slaughter_fig9.pdf
    
    
## Module membership and trait significance
  - results/ig/{imputation_group}/wgcna/inspected/module_membership_trait_significance/moduletrait_<trait>_*_<module>*.pdf
    
    
    

    
  - results/ig/{imputation_group}/wgcna/pathway_enrichment/pathway_*.pdf
    
  
    
  - results/ig/{imputation_group}/wgcna/species_table/tax_binomial_*.pdf
    
