##' Compute weight that represent how similar samples are to the reference samples, from the PCA. These weights are used to adjust the stringency of the threshold for the calls. Indeed, some subtle 'abnormal' signal detected might be due to fact that a sample is globally different from the reference samples, hence a more stringent threshold is advised. 
##' @title Weight sample difference from a set of reference samples
##' @param pca.mat a matrix with the principal components as columns and the sample names as row names. The first two components are used.
##' @param ref.samples the names of the samples used as reference.
##' @param plot should some graphs be displayed. Default is FALSE.
##' @return a vector with the weights for each sample.
##' @author Jean Monlong
##' @export
weight.ref.pca <- function(pca.mat, ref.samples, plot=FALSE){
    if(ncol(pca.mat)>2){
        pca.mat = pca.mat[,1:2]
    }

    weight.f <- function(x, min.x, max.x){
        w = (x-min.x)/(max.x-min.x)
        return(sqrt(max(c(min(c(w,1)),0))))
    }

    ref.samples = intersect(ref.samples, rownames(pca.mat))
    if(length(ref.samples)==0){
        stop("Inconsistent samples in 'ref.samples' and 'pca.mat' row names.")
    }
    
    centroid.ref = apply(pca.mat[ref.samples,],2,median, na.rm=TRUE)
    d.cent = sqrt(rowSums((pca.mat-matrix(centroid.ref, nrow(pca.mat),ncol=2, byrow=TRUE))^2))
    d.cent.med = median(d.cent[ref.samples])
    w.pca = sapply(d.cent, weight.f, min.x=d.cent.med, max.x=2*d.cent.med)

    if(plot){
        pc.df = as.data.frame(pca.mat)
        pc.df$sample = rownames(pca.mat)
        pc.df$ref = pc.df$sample %in% ref.samples
        pc.df$d.cent = d.cent
        print(ggplot2::ggplot(pc.df, ggplot2::aes(x=PC1, y=PC2, shape=ref, colour=factor(ceiling(d.cent/d.cent.med)))) + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::theme(legend.position="bottom") + ggplot2::annotate(geom="point",x=centroid.ref[1],y=centroid.ref[2], size=3, shape=8) + ggplot2::scale_colour_hue(name="distance to the reference centroid relative to median distance"))
        print(ggplot2::ggplot(subset(pc.df, ref), ggplot2::aes(x=PC1, y=PC2, shape=ref, colour=factor(ceiling(d.cent/d.cent.med)))) + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::theme(legend.position="bottom") + ggplot2::annotate(geom="point",x=centroid.ref[1],y=centroid.ref[2], size=3, shape=8) + ggplot2::scale_colour_hue(name="distance to the reference centroid relative to median distance"))
    }

    return(w.pca)
}
