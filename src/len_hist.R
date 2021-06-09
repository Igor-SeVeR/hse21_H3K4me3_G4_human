# we use sources from lib.R, where contain all libs ans constants
source('lib.R')

# file name, with which we work
#NAME <- 'H3K4me3_H1.ENCFF254ACI.hg19'
#NAME <- 'H3K4me3_H1.ENCFF668YOE.hg19'
#NAME <- 'H3K4me3_H1.ENCFF254ACI.hg38'
#NAME <- 'H3K4me3_H1.ENCFF668YOE.hg38'
#NAME <- 'H3K4me3_H1.ENCFF254ACI.hg19.filtered'
#NAME <- 'H3K4me3_H1.ENCFF668YOE.hg19.filtered'
#NAME <- 'GSM3003539_Homo'
#NAME <- 'H3K4me3_H1.intersect_with_GSM3003539_Homo'

# code to build histograms on reads length
bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end')#, 'name', 'score') # if we build GSM3003539_Homo we need only 'chrom', 'start', 'end'. 
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

# hist(bed_df$len)

# building plot and saving it

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)