#!/usr/bin/env Rscript

## In debian R needs the packages r-base and r-base-dev. Also
## libcurl4-openssl-dev is needed for cancereffectsizeR

## Create the .Rlibs directory before running this file
.libPaths(c("./.Rlibs", .libPaths()))

options(timeout = 600)

install.packages("remotes", repos = "http://cran.us.r-project.org")

remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@*release",
  dependencies = TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17") # this is so that it works with R version 4.3

library(BiocManager)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19")