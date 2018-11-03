##' Compute the null variances. Here the null distribution is modeled as two Normal distrution centered in 0. The variance are fitted to the empirical distribution. P-values and Q-values are derived from the Z-score and this fitted null distribution. 
##' @title P-values estimation from mixture of 2 centered normal
##' @param z a vector with the Z-scores
##' @return a list with
##' \item{sigma.est.dup}{the estimated null distribution variance for positive Z-scores}
##' \item{sigma.est.del}{the estimated null distribution variance for negative Z-scores}
##' @author Jean Monlong
##' @keywords internal
fdrtool.quantile.2N <- function(z) {

  localMax <- function(x, min.max.prop = 0.1) {
    d = stats::density(x, na.rm = TRUE)
    im = 1 + which(diff(sign(diff(d$y))) == -2)
    my = max(d$y)
    max.id = im[which(d$y[im] >= min.max.prop * my)]
    max.id.o = max.id[order(d$y[max.id], decreasing = TRUE)]
    return(list(lM = d$x[max.id.o], h = d$y[max.id.o]/my))
  }
  fit2norm.sd.cens <- function(z, p0 = c(p = 0.75, s1 = 1, s2 = 2), z0) {
    z = z[abs(z) < z0]
    mix.obj <- function(p, x) {
      e <- p[1] * stats::dnorm(x/p[2])/((stats::pnorm(z0, 0, p[2]) - stats::pnorm(-z0, 0, p[2])) *
                                 p[2]) + (1 - p[1]) * stats::dnorm(x/p[3])/((stats::pnorm(z0, 0, p[3]) - stats::pnorm(-z0,
                                                                                                 0, p[3])) * p[3])
      if (any(e <= 0, na.rm = TRUE) | p[1] < 0 | p[1] > 1 | p[2]>2*z0 | p[3]>2*z0)
      Inf else -sum(log(e))
    }
    lmix2a <- stats::deriv(~ -log(p * dnorm(x/s1)/((pnorm(z0, 0, s1) - pnorm(-z0, 0, s1)) * s1) + (1 - p) * dnorm(x/s2)/((pnorm(z0, 0, s2) - pnorm(-z0, 0, s2)) * s2)), c("p", "s1", "s2"), function(x, p, s1, s2) NULL)
    mix.gr <- function(pa, x) {
      p <- pa[1]
      s1 <- pa[2]
      s2 <- pa[3]
      colSums(attr(lmix2a(x, p, s1, s2), "gradient"))
    }
    results = stats::optim(p0, mix.obj, mix.gr, x = z)
    if (results$par[1] < 0.5) {
      results$par[1] = 1 - results$par[1]
      results$par[2:3] = results$par[3:2]
    }
    results
  }
  sim2norm.sd <- function(pars, nb.sims = 1e+06) {
    if (pars["p"] < 0 | pars["p"] > 1)
    return(stats::rnorm(nb.sims, 0, pars["s1"]))
    c(stats::rnorm(pars["p"] * nb.sims, 0, pars["s1"]), stats::rnorm((1 - pars["p"]) * nb.sims,
                                                       0, pars["s2"]))
  }
  find.par <- function(z) {
    findPar <- function(z0){
      do.call(rbind, lapply(z0,function(z0){
        p = fit2norm.sd.cens(z, z0 = z0)
        zsim = sim2norm.sd(p$par)
        dz = stats::density(z[abs(z) < z0], from = -z0, to = z0, n = 512)
        dzsim = stats::density(zsim[abs(zsim) < z0], from = -z0, to = z0, n = 512)
        data.frame(z0 = z0, dens.diff = sum(abs(dz$y - dzsim$y))/sum(dzsim$y),
                   p = p$par[1], s1 = p$par[2], s2 = p$par[3])
      }))
    }
    step = c(.2,.5,1)
    continue = TRUE
    cpt = 0
    first.scan = findPar(2:8)
    z0 = first.scan$z0[which.min(first.scan$dens.diff)]
    p.c = first.scan[which.min(first.scan$dens.diff),]
    while(continue & cpt < 8){
      p.l = findPar(z0-step)
      p.u = findPar(z0+step)
      if(sum(p.l$dens.diff<p.c$dens.diff)>sum(p.u$dens.diff<p.c$dens.diff)){
        z0 = max(z0-step[sum(p.l$dens.diff<p.c$dens.diff)],2)
        p.c = p.l[sum(p.l$dens.diff<p.c$dens.diff),]
      } else if(sum(p.l$dens.diff<p.c$dens.diff)<sum(p.u$dens.diff<p.c$dens.diff)){
        z0 = z0+step[sum(p.u$dens.diff<p.c$dens.diff)]
        p.c = p.u[sum(p.u$dens.diff<p.c$dens.diff),]
      } else {
        continue = FALSE
      }
      cpt = cpt + 1
    }
    list(par = unlist(p.c))
  }
  p2norm <- function(z, pars) {
    pars["p"] * stats::pnorm(z, 0, pars["s1"]) + (1 - pars["p"]) * stats::pnorm(z, 0, pars["s2"])
  }

  res = list(sigma.est.dup = NA, sigma.est.del = NA)
  z = z[which(!is.infinite(z) & !is.na(z) & z != 0)]

  sup.ss = 50000
  ## Duplication
  z.dup = z[z > 0]
  z.dup = sample(c(-1, 1), length(z.dup), replace = TRUE) * z.dup
  if (length(z.dup) > sup.ss) {
    p = find.par(sample(z.dup, sup.ss))
  } else {
    p = find.par(z.dup)
  }
  res$sigma.est.dup = p$par["s1"]
  ## Deletion
  z.del = z[z < 0]
  z.del = sample(c(-1, 1), length(z.del), replace = TRUE) * z.del
  if (length(z.del) > sup.ss) {
    p = find.par(sample(z.del, sup.ss))
  } else {
    p = find.par(z.del)
  }
  res$sigma.est.del = p$par["s1"]

  return(res)
}
