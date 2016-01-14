library(PopSV)

## Load some annotation to play with
library(AnnotationHub)
ah = AnnotationHub()
ah2 = query(ah, c("Homo sapiens", "hg19", "Genecode"))
genes = ah2[[1]]
ah3 = query(ah, c("Homo sapiens", "hg19", "peak", "narrow"))
peaks = ah3[[4]]

## Build controls regions that fit the peaks size and overlap with genes
peaks.cont = draw.controls(peaks, list(gene=genes), chr.prefix="chr")
