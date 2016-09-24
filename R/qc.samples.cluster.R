##' Interactive representation of the bin counts with potential clustering of the samples. This app can be used to check if the reference samples are homogeneous, and define outliers or batches.
##'
##' There is no need to use all the bin count, a few thousand bins is enough to detect potential problems. 'quick.count' is useful to retrieve bin counts on a subset of bins.
##' @title QC sample : clustering
##' @param bc.df a data.frame with the bin counts.
##' @param nb.rand.bins if non-NULL, the number of bins to randomly choose. Default is NULL.
##' @param shiny.app If TRUE (default), a shiny app is launched in a web browser. If FALSE, the function simply returns the PCA results.
##' @return a shiny app in the web browser, or a data.frame with the results of the PCA (if shiny.app=FALSE).
##' @author Jean Monlong
##' @export
qc.samples.cluster <- function(bc.df, nb.rand.bins=NULL, shiny.app=TRUE){
  if(!is.null(nb.rand.bins)){
    bins.df = bins.df[sample.int(nrow(bins.df), min(nrow(bins.df), nb.rand.bins)),]
  }
  samples = setdiff(colnames(bc.df), c("chr","start","end"))

  cat("Computing Principal Components...\n")
  bc.mv = medvar.norm.internal(bc.df[, samples])
  ## PCA
  pc = stats::prcomp(t(stats::na.exclude(bc.mv)))
  pc.df = data.frame(pc$x[, 1:3])
  pc.df$sample = samples
  cat("Done.\n")

  if(shiny.app){
      PC1 = PC2 = x= y = xend = yend = sample = NULL ## Uglily appease R checks
      shiny::runApp(list(

          ui = shiny::fluidPage(
              shiny::titlePanel("PopSV - QC sample clustering"),
              shiny::sidebarPanel(
                  shiny::numericInput("nb.clust", "Number of clusters",1, 1, 10, 1),
                  shiny::selectInput("cl.meth", "Clustering method", c("complete","average","ward")),
                  shiny::numericInput("cl.sel", "Cluster selected",1, 1, 10, 1),
                  shiny::hr(),
                  shiny::actionButton("exp","Done"),
                  shiny::hr(),
                  shiny::textOutput("export")
                  ),
              shiny::mainPanel(shiny::plotOutput("pca", height=800))
              ),

          server = function(input, output) {
              samples.ref <- shiny::reactive({
                  d.pca = stats::dist(pc.df[, c("PC1","PC2")])
                  hc.o = stats::hclust(d.pca, method=input$cl.meth)
                  samp.cl = stats::cutree(hc.o, input$nb.clust)
                  data.frame(pc.df,
                             cluster=samp.cl,
                             selected = samp.cl == input$cl.sel)
              })
              output$pca = shiny::renderPlot({
                  samp.df = samples.ref()
                  pdf = merge(pc.df, samp.df)
                  pdf$cluster = factor(pdf$cluster)

                  cluster = selected = NULL ## Uglily appeases R checks
                  ggplot2::ggplot(pdf, ggplot2::aes(x=PC1, y=PC2, colour=cluster, size=selected)) +
                      ggplot2::geom_point(alpha=.6) + ggplot2::theme_bw() +
                          ggplot2::scale_size_manual(values=2:3)
              })
              output$export = shiny::renderText({
                  samp.df = shiny::isolate(samples.ref())
                  if(input$exp){
                      shiny::stopApp(samp.df)
                  }
                  return(paste(sum(samp.df$selected),"samples selected."))
              })
          }))
  } else {
      return(pc.df)
  }
}
