# At first we set working directory

setwd('C:/minorBio/final_project_git/hse21_H3K4me3_G4_human/src')

# libraries

# at first, we install all libraries, which are needed to make the project

# method to check, if library is already installed and if not -> install
package.check <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}

# Installing BiocManager
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ChIPseeker")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force=TRUE)
#BiocManager::install("clusterProfiler", force=TRUE)
#BiocManager::install("GenomicFeatures", force=TRUE)
#BiocManager::install("org.Hs.eg.db", force=TRUE)

# Libraries from BiocManager
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(GenomicFeatures)
library(org.Hs.eg.db)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# constants

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'