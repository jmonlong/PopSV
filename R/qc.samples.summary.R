##' Interactive representation of the QC metrics on the reference samples. 
##' @title QC sample summary
##' @param qc.res a list D statistic, PCA and correlation information, produced
##' by 'qc.samples'.
##' @return a shiny app in the web browser
##' @author Jean Monlong
##' @export
qc.samples.summary <- function(qc.res){
        shiny::runApp(list(

        ui = shiny::fluidPage(
            shiny::headerPanel("PopSV - Sample QC"),
            shiny::sidebarPanel(
                shiny::textInput("d", "D statistic threshold", "0.8"),
                shiny::conditionalPanel(condition = "input.conditionPanels == 'Clustering'",
                                 shiny::selectInput("cl.meth", "Linkage method: ",
                                                    c("Complete"="complete","Average"="average","Ward"="ward")))
                ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("D distribution", shiny::plotOutput("d.hist")),
                    shiny::tabPanel("Clustering", shiny::plotOutput("clust")),
                    shiny::tabPanel("PCA", shiny::plotOutput("pca")), id="conditionPanels"
                    )
                )
            ),

        server = function(input, output) {
            samples.ref <- shiny::reactive({
                subset(qc.res$dstat, Dstat >= as.numeric(input$d))
            })
            output$d.hist = shiny::renderPlot({
                plot.df = qc.res$dstat
                plot.df$reference = plot.df$sample %in% samples.ref()
                return(ggplot2::ggplot(plot.df, ggplot2::aes(x=Dstat)) + 
                       ggplot2::geom_bar(ggplot2::aes(fill=reference)) + 
                       ggplot2::theme(legend.position=c(0,1),legend.justification=c(0,1)) + 
                       ggplot2::theme_bw() +
                       ggplot2::geom_vline(xintercept=input$d,linetype=2) +
                       ggplot2::ylab("number of samples")
                       )
            })
            output$clust = shiny::renderPlot({
                hc.o = hclust(as.dist(1-qc.res$cor.pw), method=input$cl.meth)
                dd <- ggdendro::dendro_data(as.dendrogram(hc.o))
                l.df = dd$labels
                l.df$reference = l.df$label %in% samples.ref()
                return(ggplot2::ggplot(dd$segments) +
                       ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
                       ggplot2::geom_point(ggplot2::aes(x=x, y=y, colour=reference),shape=18, size=5,data=l.df) +
                       ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom") +
                       ggplot2::scale_x_discrete(labels=l.df$label) +
                       ggplot2::xlab("sample") + ggplot2::ylab("") + 
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,hjust = 1))
                       )
            })
            output$pca = shiny::renderPlot({
                plot.df = data.frame(rownames(qc.res$pc.1.3), qc.res$pc.1.3, stringsAsFactors=FALSE)
                plot.df$reference = plot.df$sample %in% samples.ref()
                return(ggplot2::ggplot(plot.df, ggplot2::aes(x=PC1, y=PC2)) + 
                       ggplot2::geom_point(ggplot2::aes(colour=reference, size=reference)) + 
                       ggplot2::theme(legend.position="bottom") + 
                       ggplot2::theme_bw()
                       )                
            })
        }))

}
