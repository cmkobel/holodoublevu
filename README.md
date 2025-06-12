# holodoublevu 🪢👀



[![DOI](https://zenodo.org/badge/840225002.svg)](https://doi.org/10.5281/zenodo.15646870)



_"HOlistic LOw-level Denoising & OUtputting for BLEnded multi-View Unification"_ - (Sorry that backronym was made with chatgpt. But Holodoublevu is mine!)


> I'm a pipeline kinda guy.

This is a snakemake workflow that contains (nearly) all steps to perform the proteomics network analysis and interpretation. 

Hopefully this workflow is very _flexible_. But not so much that it false apart when being held up against the sun.


## IO

### Inputs

  - Sample metadata. (key `animal`)
  - Proteome amino acid sequences. (key `protein`)
  - Proteome intensities. (key `protein`)
  - Proteome-to-genome mapping. (key `protein`)
  

### Outputs

  - Filtered and transformed intensities. (key `protein`)
  - Imputed intensities. (key `protein`)
  - Function/intensity main overview table. (key `protein`)
  - Sample pathway enrichment. (key `sample`x`pathway`)
  - Sample species overview. (key `sample`x`pathway`)
  - Ordination plots (PCA, NMDS, UMAP). 
  - Network analysis (WGCNA). (key `module eigengene`)
  - WGCNA module pathway enrichment. (key `module eigengene`x`pathway`)
  - Other?


## Roadmap

[Busy](https://www.youtube.com/watch?v=oPQ3o14ksaM)!

## Run

Have conda/[mamba](https://github.com/conda-forge/miniforge#install). And set the principal inputs in config/config.yaml. Then do:

```bash
mamba env create -f environment.yaml
conda activate holodoublevu
snakemake
```


---
[<img width="666" alt="Screenshot 2024-08-14 at 23 32 17" src="https://github.com/user-attachments/assets/a9bc126d-f788-43f6-9d1c-ceff7d2ca498">](https://www.youtube.com/watch?v=OQSNhk5ICTI)
