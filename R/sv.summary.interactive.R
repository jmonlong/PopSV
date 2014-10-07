##' Interactive summary of the individual bins with abnormal coverage as detected per 'call.abnomal.cov'. 
##' @title Interactive summary
##' @param res.df the data.frame with the results.
##' @return a shiny app in the web browser. 
##' @author Jean Monlong
##' @export
sv.summary.interactive <- function(res.df){
    shiny::runApp(list(


        ui = shiny::fluidPage(
            shiny::headerPanel("PopSV - Results summary"),
            shiny::sidebarPanel(
                shiny::textInput("fdr", "False Discovery rate", "0.05"),
                conditionalPanel(condition = "input.conditionPanels == 'Copy number estimates'",
                                 shiny::numericInput("nbc", "Minimum number of consecutive bins", 3, 0, Inf, 1),
                                 shiny::helpText("This take a bit of time to update, please be patient."))
                ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("Number of calls", shiny::plotOutput("nb.calls")),
                    shiny::tabPanel("Copy number estimates", shiny::plotOutput("cn")),
                    shiny::tabPanel("Frequency distribution", shiny::plotOutput("freq")), id="conditionPanels"
                    )
                )
            ),


        server = function(input, output) {
            plot.df <- reactive({
                subset(res.df, qv<as.numeric(input$fdr))
            })

            res.m <- reactive({
                dplyr::do(dplyr::group_by(plot.df(), sample),mergeConsBin.simple(.))
            })
            output$nb.calls = shiny::renderPlot({
                ggplot2::ggplot(plot.df(), ggplot2::aes(x=sample)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw()
            })
            output$cn = shiny::renderPlot({
                ggplot2::ggplot(subset(res.m(), nbBinCons>input$nbc), ggplot2::aes(x=cn.coeff*2)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw() + ggplot2::xlim(0,5)
            })
            output$freq = shiny::renderPlot({
                p.df = plot.df()
                nb.samp = length(unique(p.df$sample))
                freq.df = dplyr::summarize(dplyr::group_by(p.df, chr, start),
                    freq.n = length(start), freq.p = length(start)/nb.samp)
                ggplot2::ggplot(freq.df, ggplot2::aes(x=freq.n)) +
                    ggplot2::geom_bar() + ggplot2::theme_bw()
            })
        }))
}
