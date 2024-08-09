# holodoublevu 🪢👀

> I'm a pipeline kinda guy.

## IO

### Inputs

  - Sample metadata. (key `sample`)
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
# conda env create -f environment.yaml
# conda activate holodoublevu
snakemake
```
