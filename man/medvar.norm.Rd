% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medvar.norm.R
\name{medvar.norm}
\alias{medvar.norm}
\title{Median-variance normalization of bin counts}
\usage{
medvar.norm(bc, ref.samples, bc.support = NULL, z.poisson = FALSE)
}
\arguments{
\item{bc}{a matrix or data.frame with the bin counts (bin x sample).}

\item{ref.samples}{a vector with the names of the samples to be used as reference.}

\item{bc.support}{TODO}

\item{z.poisson}{Should the Z-score be computed as an normal-poisson hybrid (see
Details). Default is FALSE.}
}
\value{
a matrix with the normalized bin counts (bin x sample).
}
\description{
The median and deviation to the median is normalized across samples.
}
\author{
Jean Monlong
}
