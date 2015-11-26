##' Chromosomes with clear aneuploidy are detected using a simple threshold-based approach. It is used to flag chromosome with complete aneuploidy that could handicap later analysis. However it might not detect partial chromosomal aberration or aneuploidy in samples with noisy read coverage.
##' @title Flag chromosomal aneuploidy
##' @param samp the name of the sample to analyze.
##' @param files.df a data.frame with information about samples such as the path to the count files.
##' @param col.file the name of the column in 'files.df' with the path information to use. Default is 'bc.gc.gz', i.e. GC corrected bin counts.
##' @param nb.bins the number of bins to use when subsampling each chromosome. Default is 1000.
##' @param prop.aneu the minimum proportion of normal bins at which a chromosome is considered normal.
##' @param plot should the read coverage per chromosome be displayed. Default is FALSE.
##' @return a list with:
##' \item{aneu.chrs}{vector with the names of the flagged(aneuploid) chromosomes.}
##' \item{aneu.chrs.fc}{the fold change in read count in the aneuploid chromosomes.}
##' @author Jean Monlong
##' @export
aneuploidy.flag <- function(samp, files.df, col.file = "bc.gc.gz", nb.bins = 1000, prop.aneu = 0.1, plot = FALSE) {
    bc = aneu.flag = chr = . = NULL  ## Appease R checks

    localMax <- function(x, min.max.prop = 0.1) {
        d = density(x, na.rm = TRUE)
        im = 1 + which(diff(sign(diff(d$y))) == -2)
        my = max(d$y)
        max.id = im[which(d$y[im] >= min.max.prop * my)]
        max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
        return(list(lM = d$x[max.id.o], h = d$y[max.id.o]/my))
    }

    df = read.table(files.df[files.df$sample == samp, col.file], header = TRUE, as.is = TRUE)
    chrs = unique(df$chr)
    df = df[which(df$bc > 0), ]

    df.sub = df[unlist(tapply(1:nrow(df), df$chr, function(ii)sample(ii,nb.bins))),]
    lm.o = localMax(df.sub$bc)$lM[1]
    core.chrs = df.sub$chr[order(abs(lm.o - df.sub$bc))[1:(nrow(df.sub) * 0.3)]]
    cchrs.t = table(core.chrs)
    aneu.chrs = NULL
    if (length(cchrs.t) < length(chrs)) {
        aneu.chrs = setdiff(chrs, names(cchrs.t))
    }
    if (any(cchrs.t < prop.aneu * nb.bins * 0.3)) {
        aneu.chrs = c(aneu.chrs, names(cchrs.t)[which(cchrs.t < prop.aneu * nb.bins *
            0.3)])
    }

    aneu.chrs.fc = NULL
    if (!is.null(aneu.chrs)) {
        aneu.chrs.fc = sapply(aneu.chrs, function(chr.i) localMax(df.sub$bc[which(df.sub$chr ==
            chr.i)])$lM[1]/lm.o)
    }

    if (plot) {
        df$aneu.flag = df$chr %in% aneu.chrs
        df$bc = winsor(df$bc, u = quantile(df$bc, probs = 0.995, na.rm = TRUE))
        p1 = ggplot2::ggplot(df, ggplot2::aes(x = bc, fill = aneu.flag)) + ggplot2::geom_density(alpha = 0.4) +
            ggplot2::theme_bw() + ggplot2::facet_grid(chr ~ ., scales = "free") +
            ggplot2::xlab("read coverage") + ggplot2::ylab("bin density") + ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
            ggplot2::theme(legend.position = "bottom")
        df$all = "all"
        p2 = ggplot2::ggplot(df, ggplot2::aes(x = bc)) + ggplot2::geom_density(fill = "grey70",
            alpha = 0.4) + ggplot2::theme_bw() + ggplot2::ylab("bin density") + ggplot2::xlab("read coverage") +
            ggplot2::theme(axis.text.y = ggplot2::element_blank()) + ggplot2::facet_grid(all ~
            ., scales = "free") + ggplot2::ggtitle(samp)
        grid::grid.newpage()
        print(p1, vp = grid::viewport(1, 0.8, x = 0.5, y = 0.4))
        print(p2, vp = grid::viewport(1, 0.2, x = 0.5, y = 0.9))
    }

    return(list(aneu.chrs = aneu.chrs, aneu.chrs.fc = aneu.chrs.fc))
}
