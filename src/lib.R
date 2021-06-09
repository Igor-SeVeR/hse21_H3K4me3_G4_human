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
#BiocManager::install("ChIPpeakAnno", force=TRUE)



# Libraries from BiocManager
package.check("ChIPseeker")
package.check("TxDb.Hsapiens.UCSC.hg19.knownGene")
package.check("clusterProfiler")
package.check("GenomicFeatures")
package.check("org.Hs.eg.db")
package.check("ChIPpeakAnno")

package.check("ggplot2")
package.check("dplyr")
package.check("tidyr")
package.check("tibble")

# constants

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'