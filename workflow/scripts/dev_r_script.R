#!/usr/bin/env Rscript

library(tidyverse)

paste(snakemake@input)
paste(snakemake@output["thing"]) # named IOs

file.create(snakemake@output["thing"] %>% as.character())
