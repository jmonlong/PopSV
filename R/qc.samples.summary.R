##' Interactive representation of the QC metrics on the reference samples. 
##' @title QC sample summary
##' @param qc.res a list D statistic, PCA and correlation information, produced
##' by 'qc.samples'.
##' @return a shiny app in the web browser
##' @author Jean Monlong
##' @export
qc.samples.summary <- function(qc.res){
  PC1 = PC2 = Dstat = reference = x= y = xend = yend = sample = NULL ## Uglily appease R checks

  shiny::runApp(list(
        
        ui = shiny::fluidPage(
            shiny::titlePanel("PopSV - Sample QC"),
            shiny::sidebarPanel(
                shiny::textInput("d", "D statistic threshold", "0.8"),
                shiny::textInput("ols", "Outlier samples (',' separated)", ""),
                shiny::conditionalPanel(condition = "input.conditionPanels == 'PCA'",
                                        shiny::radioButtons("samp.lab", "Point : ",
                                                            c("Point","Sample name"))),
                shiny::conditionalPanel(condition = "input.conditionPanels != 'D distribution'",
                                        shiny::radioButtons("samp.set", "Samples : ",
                                                            c("All","Reference"))),
                shiny::conditionalPanel(condition = "input.conditionPanels == 'Clustering'",
                                        shiny::selectInput("cl.meth", "Linkage method: ",
                                                           c("Average"="average","Complete"="complete","Ward"="ward.D"))),
                shiny::hr(),
                shiny::textOutput("nbsamp")
                ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("D distribution", shiny::plotOutput("d.hist")),
                    shiny::tabPanel("Clustering", shiny::plotOutput("clust")),
                    shiny::tabPanel("PCA", shiny::plotOutput("pca")),
                    shiny::tabPanel("Export", shiny::actionButton("exp","Export to 'ref-samples.RData'"),shiny::hr(),shiny::textOutput("export")),id="conditionPanels"
                    )
                )),
            
        server = function(input, output) {
            samples.ref <- shiny::reactive({
                sr = as.character(qc.res$dstat$sample[which(qc.res$dstat$Dstat >= as.numeric(input$d))])
                if(input$ols != ""){
                    ols = unlist(strsplit(input$ols,","))
                    return(setdiff(sr, ols))
                } else {
                    return(sr)
                }
            })
            output$d.hist = shiny::renderPlot({
                plot.df = qc.res$dstat
                plot.df$reference = factor(plot.df$sample %in% samples.ref(), levels=c("TRUE","FALSE"))
                return(ggplot2::ggplot(plot.df, ggplot2::aes(x=Dstat)) + 
                       ggplot2::geom_histogram(ggplot2::aes(fill=reference), binwidth=.005) + 
                       ggplot2::scale_fill_brewer(palette="Set1") +
                       ggplot2::theme(legend.position=c(0,1),legend.justification=c(0,1)) + 
                       ggplot2::theme_bw() +
                       ggplot2::geom_vline(xintercept=as.numeric(input$d),linetype=2) +
                       ggplot2::ylab("number of samples"))
            })
            output$clust = shiny::renderPlot({
                refs = samples.ref()
                if(input$samp.set=="Reference"){
                    hc.o = hclust(as.dist(1-qc.res$cor.pw[refs,refs]), method=input$cl.meth)
                } else {
                    hc.o = hclust(as.dist(1-qc.res$cor.pw), method=input$cl.meth)
                }
                dd <- ggdendro::dendro_data.hclust(hc.o)
                l.df = dd$labels
                l.df$reference = factor(l.df$label %in% refs, levels=c("TRUE","FALSE"))
                return(ggplot2::ggplot(dd$segments) +
                       ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
                       ggplot2::geom_point(ggplot2::aes(x=x, y=y, colour=reference),shape=18, size=5,data=l.df) +
                       ggplot2::scale_colour_brewer(palette="Set1") +
                       ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom") +
                       ggplot2::scale_x_discrete(labels=l.df$label) +
                       ggplot2::xlab("sample") + ggplot2::ylab("") + 
                       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,hjust = 1))
                       )
            })
            output$pca = shiny::renderPlot({
                plot.df = data.frame(sample=rownames(qc.res$pc.1.3), qc.res$pc.1.3, stringsAsFactors=FALSE)
                plot.df$reference = factor(plot.df$sample %in% samples.ref(), levels=c("TRUE","FALSE"))
                if(input$samp.set=="Reference") plot.df = plot.df[which(plot.df$reference=="TRUE"),]
                xlims = c(min(plot.df$PC1,na.rm=TRUE),max(plot.df$PC1,na.rm=TRUE))
                xlims = c(xlims[1]-.07*diff(xlims),xlims[2]+.07*diff(xlims))
                gpl = ggplot2::ggplot(plot.df, ggplot2::aes(x=PC1, y=PC2, colour=reference)) + 
                       ggplot2::theme(legend.position="bottom") + 
                           ggplot2::scale_colour_brewer(palette="Set1") +
                               ggplot2::theme_bw() + ggplot2::xlim(xlims[1],xlims[2])
                if(input$samp.lab=="Point") return(gpl + ggplot2::geom_point(alpha=.7))
                else return(gpl + ggplot2::geom_text(ggplot2::aes(label=sample)))
            })
            output$export = shiny::renderText({
                input$exp
                ref.samples = shiny::isolate(samples.ref())
                save(ref.samples, file="ref-samples.RData")
                return(paste(length(ref.samples)," samples saved in 'ref-samples.RData'."))
            })
            output$nbsamp = shiny::renderText({
                ref.samples = samples.ref()
                return(paste(length(ref.samples)," reference samples selected."))
            })
        }))

}
