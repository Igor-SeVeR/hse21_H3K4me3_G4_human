# we use sources from lib.R, where contain all libs ans constants
source('lib.R')

# file name, with which we work
#NAME <- 'H3K4me3_H1.ENCFF254ACI.hg19.filtered'
#NAME <- 'H3K4me3_H1.ENCFF668YOE.hg19.filtered'
#NAME <- 'GSM3003539_Homo'

BED_FN <- paste0(DATA_DIR, NAME, '.bed')

###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(BED_FN, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.pdf'))
pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.pdf'))
plotAnnoPie(peakAnno)
dev.off()

# if working with GSM3003539_Homo -> comment this part
peak <- readPeakFile(BED_FN)
pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.covplot.pdf'))
covplot(peak)
dev.off()
 