##' Interactive representation of the bin counts with potential clustering of the samples. This app can be used to check if the reference samples are homogeneous, and define outliers or batches.
##'
##' There is no need to use all the bin count, a few thousand bins is enough to detect potential problems. 'quick.count' is useful to retrieve bin counts on a subset of bins.
##' @title QC sample : clustering
##' @param bc.df a data.frame with the bin counts.
##' @param nb.rand.bins if non-NULL, the number of bins to randomly choose. Default is NULL.
##' @return a shiny app in the web browser
##' @author Jean Monlong
##' @export
qc.samples.cluster <- function(bc.df, nb.rand.bins=NULL){
  if(!is.null(nb.rand.bins)){
    bins.df = bins.df[sample.int(nrow(bins.df), min(nrow(bins.df), nb.rand.bins)),]
  }
  samples = setdiff(colnames(bc.df), c("chr","start","end"))

  cat("Computing Principal Components...\n")
  bc.mv = medvar.norm.internal(bc.df[, samples])
  ## PCA
  pc = prcomp(t(na.exclude(bc.mv)))
  pc.df = data.frame(pc$x[, 1:3])
  pc.df$sample = samples
  cat("Done.\n")

  PC1 = PC2 = x= y = xend = yend = sample = NULL ## Uglily appease R checks
  shiny::runApp(list(

      ui = shiny::fluidPage(
          shiny::titlePanel("PopSV - QC sample clustering"),
          shiny::sidebarPanel(
              shiny::numericInput("nb.clust", "Number of clusters",1, 1, 10, 1),
              shiny::selectInput("cl.meth", "Clustering method", c("complete","average","ward")),
              shiny::numericInput("cl.sel", "Cluster selected",1, 1, 10, 1),
              shiny::hr(),
              shiny::actionButton("exp","Export selected cluster to 'ref-samples.RData'"),
              shiny::hr(),
              shiny::textOutput("export"),
              shiny::hr(),
              shiny::actionButton("expAll","Export all clusters to 'samples-clusters.RData'"),
              shiny::hr(),
              shiny::textOutput("exportAll")
              ),
          shiny::mainPanel(shiny::plotOutput("pca", height=800))
          ),

      server = function(input, output) {
        samples.ref <- shiny::reactive({
            d.pca = dist(pc.df[, c("PC1","PC2")])
            hc.o = hclust(d.pca, method=input$cl.meth)
            samp.cl = cutree(hc.o, input$nb.clust)
            data.frame(sample=pc.df$sample,
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
            input$exp
            samp.df = shiny::isolate(samples.ref())
            ref.samples = as.character(samp.df$sample[which(samp.df$selected)])
            save(ref.samples, file="ref-samples.RData")
            return(paste(length(ref.samples)," samples saved in 'ref-samples.RData'."))
          })
        output$exportAll = shiny::renderText({
            input$expAll
            samp.df = shiny::isolate(samples.ref())
            samp.df$selected = NULL
            save(samp.df, file="samples-clusters.RData")
            return(paste(length(unique(samp.df$cluster))," clusters saved in 'samples-clusters.RData'."))
          })
      }))

}
