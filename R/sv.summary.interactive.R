##' Interactive summary of the individual bins with abnormal coverage as detected per 'call.abnomal.cov'. 
##' @title Interactive summary
##' @param res.df the data.frame with the results.
##' @param merge.cons.bin TRUE if the single-bin calls has been merged
##' (Default, 'merge' option in 'call.abnormal.cov'). FALSE for single-bin calls. 
##' @return a shiny app in the web browser. 
##' @author Jean Monlong
##' @export
sv.summary.interactive <- function(res.df, merge.cons.bin=TRUE){
    freq.gr <- function(cnv.o){
        gr =  with(cnv.o, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
        gr.d = GenomicRanges::disjoin(gr)
        ol = GenomicRanges::findOverlaps(gr.d, gr)
        freq.df = as.data.frame(table(IRanges::queryHits(ol)))
        colnames(freq.df) = c("freq.n", "nb")
        freq.df$prop = freq.df$n / length(unique(cnv.o$sample))
        freq.df
    }

    res.df$cnDev = abs(res.df$cn.coeff-1)
    
    if(merge.cons.bin){
###
### MERGED BIN CALLS
###
        
        shiny::runApp(list(
            
            ui = shiny::fluidPage(
                shiny::headerPanel("PopSV - Results summary"),
                shiny::sidebarPanel(
                    shiny::textInput("fdr", "False Discovery rate", "0.05"),
                    shiny::textInput("cnD", "Deviation from CN 2", "0"),
                    shiny::conditionalPanel(condition = "input.conditionPanels != 'Frequency distribution'",
                                            shiny::radioButtons("col", "Colour by ", c("event type","event size","sample"))),
                    shiny::conditionalPanel(condition = "input.conditionPanels == 'Copy number estimates'",
                                            shiny::numericInput("nbc", "Minimum number of consecutive bins", 3, 0, Inf, 1),
                                            shiny::numericInput("cnMin", "Minimum CN shown", 0, 0, Inf, 1),
                                            shiny::numericInput("cnMax", "Maximum CN shown", 5, 1, Inf, 1)),
                    shiny::conditionalPanel(condition = "input.conditionPanels == 'Frequency distribution'",
                                            shiny::numericInput("nbMin", "Minimum number of samples shown", 0, 0, Inf, 1))
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
                plot.df <- shiny::reactive({
                    subset(res.df, qv<as.numeric(input$fdr) & cnDev>=input$cnD)
                })
                
                output$nb.calls = shiny::renderPlot({
                    pdf = plot.df()
                    pdf$sample = factor(pdf$sample, levels=names(sort(table(pdf$sample))))
                    if(input$col=="event type"){
                        pdf$col = ifelse(pdf$cn.coeff>1, "duplication","deletion")
                        extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
                    } else if(input$col=="event size"){
                        pdf$col = cut(pdf$nbBinCons, breaks=c(0,1,2,3,5,10,Inf))
                        extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
                    } else if(input$col=="sample"){
                        pdf$col = as.character(pdf$sample)
                        extra.gg = ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),ceiling(length(unique(pdf$sample))/9)))
                    }
                    gp = ggplot2::ggplot(pdf, ggplot2::aes(x=sample, fill=col)) +
                        ggplot2::geom_bar() + ggplot2::theme_bw() +
                            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90),
                                           legend.position=c(0,1),legend.justification=c(0,1)) + extra.gg + ggplot2::ylab("number of calls")
                    if(input$col=="sample"){
                        gp = gp  + ggplot2::guides(fill=FALSE)
                    }
                    gp
                })
                output$cn = shiny::renderPlot({
                    pdf = subset(plot.df(), nbBinCons>=input$nbc)
                    if(input$col=="event type"){
                        pdf$col = ifelse(pdf$cn.coeff>1, "duplication","deletion")
                        extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
                    } else if(input$col=="event size"){
                        pdf$col = cut(pdf$nbBinCons, breaks=c(0,1,2,3,5,10,Inf))
                        extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
                    } else if(input$col=="sample"){
                        pdf$col = as.character(pdf$sample)
                        extra.gg = ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),ceiling(length(unique(pdf$sample))/9)))
                    }
                    gp = ggplot2::ggplot(pdf, ggplot2::aes(x=cn.coeff*2, fill=col)) +
                        ggplot2::geom_bar() + ggplot2::theme_bw() + ggplot2::xlim(input$cnMin,input$cnMax)+
                            ggplot2::theme(legend.position=c(1,1),legend.justification=c(1,1)) + extra.gg + ggplot2::xlab("Copy Number estimates") + ggplot2::ylab("number of calls")
                    if(input$col=="sample"){
                        gp = gp  + ggplot2::guides(fill=FALSE)
                    }
                    gp
                })
                output$freq = shiny::renderPlot({
                    freq.df = freq.gr(plot.df())
                    ggplot2::ggplot(subset(freq.df, nb>=input$nbMin), ggplot2::aes(x=nb)) +
                        ggplot2::geom_bar() + ggplot2::theme_bw() +
                            ggplot2::ylab("number of unique genomic region") +
                                ggplot2::xlab("number of samples")
                })
            }))
        
        
    } else {
###
### SINGLE BIN CALLS
###
        
        shiny::runApp(list(
            
            
            ui = shiny::fluidPage(
                shiny::headerPanel("PopSV - Results summary"),
                shiny::sidebarPanel(
                    shiny::textInput("fdr", "False Discovery rate", "0.05"),
                    shiny::conditionalPanel(condition = "input.conditionPanels == 'Copy number estimates'",
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
                plot.df <- shiny::reactive({
                    subset(res.df, qv<as.numeric(input$fdr))
                })

                res.m <- shiny::reactive({
                    dplyr::do(dplyr::group_by(plot.df(), sample),mergeConsBin.simple(.))
                })
                output$nb.calls = shiny::renderPlot({
                    pdf = plot.df()
                    pdf$sample = factor(pdf$sample, levels=names(sort(table(pdf$sample))))
                    ggplot2::ggplot(pdf, ggplot2::aes(x=sample)) +
                        ggplot2::geom_bar() + ggplot2::theme_bw() +
                            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90))
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
}
