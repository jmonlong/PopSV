##' The median and deviation to the median is normalized across samples. 
##' @title Median-variance normalization of bin counts
##' @param bc a matrix or data.frame with the bin counts (bin x sample).
##' @param ref.samples a vector with the names of the samples to be used as reference.
##' @param bc.support TODO
##' @param z.poisson Should the Z-score be computed as an normal-poisson hybrid (see
##' Details). Default is FALSE.
##' @return a matrix with the normalized bin counts (bin x sample).
##' @author Jean Monlong
##' @export
medvar.norm <- function(bc, ref.samples, bc.support = NULL, z.poisson = FALSE) {
    all.samples = setdiff(colnames(bc), c("chr", "start", "end"))
    ref.samples.ii = which(all.samples %in% ref.samples)
    rownames(bc) = bins = paste(bc$chr, as.integer(bc$start), sep = "-")
    norm.stats = createEmptyDF(c("character", rep("integer", 2), rep("numeric", 3)), 
        length(bins))
    colnames(norm.stats) = c("chr", "start", "end", "m", "sd", "nb.remove")
    bc.norm = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), 
        length(bins))
    z = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), 
        length(bins))
    fc = createEmptyDF(c("character", rep("integer", 2), rep("numeric", length(all.samples))), 
        length(bins))
    colnames(bc.norm) = colnames(z) = colnames(fc) = c("chr", "start", "end", all.samples)
    norm.stats$chr = bc.norm$chr = z$chr = fc$chr = bc[bins, "chr"]
    norm.stats$start = bc.norm$start = z$start = fc$start = bc[bins, "start"]
    norm.stats$end = bc.norm$end = z$end = fc$end = bc[bins, "end"]
    if (z.poisson) {
        z.comp <- function(x, mean.c, sd.c) {
            z.n = (x - mean.c)/sd.c
            z.p = qnorm(ppois(x, mean.c))
            n.ii = abs(z.n) < abs(z.p)
            z.p[n.ii] = z.n[n.ii]
            z.p
        }
    } else {
        z.comp <- function(x, mean.c, sd.c) {
            (x - mean.c)/sd.c
        }
    }
    
    if (is.null(bc.support)) {
        bc.support = bc[, all.samples]
    } else {
        bc.support = bc.support[, all.samples]
    }
    med = apply(bc.support, 2, median, na.rm = TRUE)
    med.c = mean(med)
    bc.support = t(t(bc.support) * med.c/med)
    bc.support = bc.support - med.c
    md = apply(bc.support, 2, function(x) median(abs(x), na.rm = TRUE))
    md.c = median(abs(bc.support), na.rm = TRUE)
    bc = bc[, all.samples] - med.c
    bc.norm[, -(1:3)] = t(t(bc) * md.c/md) + med.c
    if (any(bc.norm[, -(1:3)] < 0)) {
        bc.norm[, -(1:3)][bc.norm[, -(1:3)] < 0] = 0
    }
    msd = apply(bc.norm[, ref.samples], 1, function(ee) unlist(mean.sd.outlierR(ee)))
    norm.stats[, 4:6] = cbind(msd[1, ], msd[2, ], msd[3, ])
    z[, -(1:3)] = apply(bc.norm[, -(1:3)], 2, z.comp, mean.c = norm.stats$m, sd.c = norm.stats$sd)
    fc[, -(1:3)] = bc.norm[, -(1:3)]/norm.stats$m
    
    list(norm.stats = norm.stats, bc.norm = bc.norm, z = z, fc = fc, z.poisson = z.poisson)
} 
