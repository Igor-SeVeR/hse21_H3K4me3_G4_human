# we use sources from lib.R, where contain all libs ans constants
source('lib.R')

# file name, with which we work
#NAME <- 'H3K4me3_H1.ENCFF254ACI.hg19'
#NAME <- 'H3K4me3_H1.ENCFF668YOE.hg19'

# code to read reads
bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

### sort reads by their length and create filtered reads

# checking row number before filtering
nrow(bed_df)

# filtering
bed_df <- bed_df%>%
  arrange(-len) %>%
  #filter(len < 50000) #if ENCFF254ACI
  filter(len < 9000) #if ENCFF668YOE

# checking row number after filtering
nrow(bed_df)

# printing head of reads
head(bed_df, 10)

# writing filtered file
bed_df %>%
  select(-len) %>%
  write.table(file=paste0(DATA_DIR, NAME,'.filtered.bed'),
              col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)

