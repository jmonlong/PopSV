library(PopSV)
library(ggplot2)

##
## Constructing control regions for enrichment analysis
## controlling for overlap with specific features.
##

## Load some annotation to play with
library(AnnotationHub)
ah = AnnotationHub()
genes = ah[["AH49010"]] ## Genes
dgv = ah[["AH5120"]] ## SVs from DGV
dgv = dgv[sample.int(length(dgv), 1e4)] ## Reduce to 10K random SVs

## Build controls regions that fit the SV size and overlap with genes
dgv.cont = draw.controls(dgv, list(gene=genes), chr.prefix="chr")

## Same size distribution ?
size.df = rbind(data.frame(reg="dgv", size=width(dgv)),
                 data.frame(reg="control", size=width(dgv.cont)))
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge")
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge") + scale_x_log10()

## Same overlap with genes ?
mean(overlapsAny(dgv, genes))
mean(overlapsAny(dgv.cont, genes))

##
## Now with more features
gap = ah[["AH6444"]]
segdups = ah[["AH5121"]]

## Build controls regions that fit the SV size and overlap with genes, gap and segmental duplication
dgv.cont2 = draw.controls(dgv, list(gene=genes, gap=gap, sd=segdups), chr.prefix="chr")

## Same size distribution ?
size.df = rbind(data.frame(reg="dgv", size=width(dgv)),
                 data.frame(reg="control", size=width(dgv.cont2)))
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge")
ggplot(size.df, aes(x=size, fill=reg)) + geom_histogram(position="dodge") + scale_x_log10()

## Same overlap with features ?
mean(overlapsAny(dgv, genes))
mean(overlapsAny(dgv.cont2, genes))
mean(overlapsAny(dgv, gap))
mean(overlapsAny(dgv.cont, gap)) ## No
mean(overlapsAny(dgv.cont2, gap)) ## Yes
mean(overlapsAny(dgv, segdups))
mean(overlapsAny(dgv.cont, segdups)) ## No
mean(overlapsAny(dgv.cont2, segdups)) ## Yes
