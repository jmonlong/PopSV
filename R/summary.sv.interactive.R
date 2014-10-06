##' Interactive summary of the individual bins with abnormal coverage as detected per 'call.abnomal.cov'. 
##' @title Interactive summary
##' @param res.df the data.frame with the results.
##' @return a shiny app in the web browser. 
##' @author Jean Monlong
##' @export
summary.sv.interactive <- function(res.df){
    shiny::runApp(list(
        ui = shiny::fluidPage(
            shiny::headerPanel("PopSV - Results summary"),
            shiny::sidebarPanel(
                shiny::textInput("fdr", "False Discovery rate", "0.05")
                ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("Number of calls", shiny::plotOutput("nb.calls")),
                    shiny::tabPanel("Copy number estimates", shiny::plotOutput("cn")),
                    shiny::tabPanel("Frequency distribution", shiny::plotOutput("freq"))
                    )
                )
            ),
        server = function(input, output) {
            output$nb.calls = shiny::renderPlot({
                fdr = as.numeric(input$fdr)
                ggplot2::ggplot(subset(res.df,qv<fdr), ggplot2::aes(x=sample)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw()
            })
            output$cn = shiny::renderPlot({
                fdr = as.numeric(input$fdr)
                res.m = mergeConsBin.simple(subset(res.df,qv<fdr))
                ggplot2::ggplot(subset(res.m, nbBinCons>2), ggplot2::aes(x=cn.coeff*2)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw() + ggplot2::xlim(0,5)
            })
            output$freq = shiny::renderPlot({
                fdr = as.numeric(input$fdr)
                nb.samp = length(unique(res.df$sample))
                freq.df = dplyr::summarize(dplyr::group_by(subset(res.df,qv<fdr), chr, start),
                    freq.n = dplyr::n(), freq.p = dplyr::n()/nb.samp)
                ggplot2::ggplot(freq.df, ggplot2::aes(x=freq.n)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw()
            })
        }))
}
