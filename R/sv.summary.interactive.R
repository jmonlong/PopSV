##' Interactive summary of the individual bins with abnormal coverage as detected per 'call.abnomal.cov'. 
##' @title Interactive summary
##' @param res.df the data.frame with the results.
##' @param merge.cons.bin TRUE if the single-bin calls has been merged
##' (Default, 'merge' option in 'call.abnormal.cov'). FALSE for single-bin calls. 
##' @param height the height of the plot in the webpage. Default is "500px".
##' @return a shiny app in the web browser. 
##' @author Jean Monlong
##' @export
sv.summary.interactive <- function(res.df, merge.cons.bin=TRUE, height="500px"){
  sample = gen.kb = col = cn = chr = nb = fc = type = start = end = prop = . = freq.n = NULL ## Uglily appease R checks

  nb.samp = length(unique(res.df$sample))
  freq.chr.gr <- function(cnv.o){
    gr =  with(cnv.o, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    gr.d = GenomicRanges::disjoin(gr)
    ol = GenomicRanges::findOverlaps(gr.d, gr)
    thits = base::table(IRanges::queryHits(ol))
    freq.df = data.frame(win.id = as.numeric(names(thits)), nb=as.numeric(thits))
    win.df = GenomicRanges::as.data.frame(gr.d[as.numeric(freq.df$win.id)])[,1:3]
    freq.df$chr = as.character(win.df$seqnames)
    freq.df$start = as.numeric(win.df$start)
    freq.df$end = as.numeric(win.df$end)
    freq.df$prop = freq.df$nb / nb.samp
    freq.df
  }

  res.df$gen.kb = with(res.df, (end-start)/1e3)
  
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
          shiny::textInput("sing.kb", "Maximum amount of single bins (Kb)", "Inf"),
          shiny::conditionalPanel(condition = "input.conditionPanels == 'Copy number estimates' | input.conditionPanels == 'Number of calls'",
                                  shiny::radioButtons("col", "Colour by ", c("event type","event size","sample"))),
          shiny::conditionalPanel(condition = "input.conditionPanels == 'Copy number estimates'",
                                  shiny::numericInput("nbc", "Minimum number of consecutive bins", 3, 0, Inf, 1),
                                  shiny::numericInput("cnMin", "Minimum CN shown", 0, 0, Inf, 1),
                                  shiny::numericInput("cnMax", "Maximum CN shown", 5, 1, Inf, 1)),
          shiny::conditionalPanel(condition = "input.conditionPanels == 'Frequency across the genome'",
                                  shiny::selectInput("chr","Chromosome",c("all",1:22)),shiny::radioButtons("fchr.rep","Representation",c("Stacked","Sample"))),
          shiny::conditionalPanel(condition = "input.conditionPanels == 'Frequency across the genome' | input.conditionPanels == 'Frequency distribution'",
                                  shiny::radioButtons("freq.rep","Frequency:", c("Number of samples"="nb","Proportion of samples"="prop"))),
          shiny::conditionalPanel(condition = "input.conditionPanels == 'Frequency distribution'",
                                  shiny::numericInput("nbMin", "Minimum frequency (number of samples) shown", 0, 0, Inf, 1),
                                  shiny::hr(),
                                  shiny::helpText("Colours represent chromosomes."),
                                  shiny::hr(),
                                  shiny::helpText("Frequency computation might take a few seconds."))
          ),
        shiny::mainPanel(
          shiny::tabsetPanel(
            shiny::tabPanel("Number of calls", shiny::plotOutput("nb.calls", height=height)),
            shiny::tabPanel("Copy number estimates", shiny::plotOutput("cn", height=height)),
            shiny::tabPanel("Frequency distribution", shiny::plotOutput("freq", height=height)),
            shiny::tabPanel("Frequency across the genome", shiny::plotOutput("freq.chr", height=height)), id="conditionPanels"
            )
          )
        ),
      
      
      server = function(input, output) {
        plot.df <- shiny::reactive({
          pdf = res.df[which(res.df$qv<as.numeric(input$fdr) & res.df$cn2.dev>=input$cnD),]
          sing.kb = aggregate(gen.kb~sample, data=pdf[pdf$nb.bin.cons==1,], sum)
          pdf = pdf[pdf$sample %in% sing.kb$sample[sing.kb$gen.kb<as.numeric(input$sing.kb)], ]
          samp.o = aggregate(gen.kb~sample, data=pdf, sum)
          pdf$sample = factor(pdf$sample, levels=samp.o$sample[order(samp.o$gen.kb)])
          pdf
        })
        freq.df <- shiny::reactive({
          pdf = plot.df()
          rbind(data.frame(type="deletion",freq.chr.gr(pdf[which(pdf$fc<1),])),
                data.frame(type="duplication",freq.chr.gr(pdf[which(pdf$fc>1),])))
        })
        
        output$nb.calls = shiny::renderPlot({
          pdf = plot.df()
          if(input$col=="event type"){
            pdf$col = ifelse(pdf$fc>1, "duplication","deletion")
            extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
          } else if(input$col=="event size"){
            pdf$col = cut(pdf$nb.bin.cons, breaks=c(0,1,2,3,5,10,Inf))
            extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
          } else if(input$col=="sample"){
            pdf$col = pdf$sample
            extra.gg = ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),ceiling(length(unique(pdf$sample))/9)))
          }
          pdf.s = dplyr::summarize(dplyr::group_by(pdf, sample, col), gen.kb=sum((end-start)/1e3))
          gp = ggplot2::ggplot(pdf.s, ggplot2::aes(x=sample, y=gen.kb, fill=col)) +
            ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() +
              ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90),
                             legend.position=c(0,1),legend.justification=c(0,1)) + extra.gg + ggplot2::ylab("abnormal genome (Kb)")
          if(input$col=="sample"){
            gp = gp  + ggplot2::guides(fill=FALSE)
          }
          gp
        })
        output$cn = shiny::renderPlot({
          pdf = plot.df()
          pdf = pdf[which(pdf$nb.bin.cons>=input$nbc),]
          if(input$col=="event type"){
            pdf$col = ifelse(pdf$fc>1, "duplication","deletion")
            extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
          } else if(input$col=="event size"){
            pdf$col = cut(pdf$nb.bin.cons, breaks=c(0,1,2,3,5,10,Inf))
            extra.gg = ggplot2::scale_fill_brewer(name=input$col,palette="Set1")
          } else if(input$col=="sample"){
            pdf$col = pdf$sample
            extra.gg = ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),ceiling(length(unique(pdf$sample))/9)))
          }
          pdf$cn = pdf$fc*2
          if(any(pdf$cn>input$cnMax)) pdf$cn[pdf$cn>input$cnMax] = input$cnMax
          gp = ggplot2::ggplot(pdf, ggplot2::aes(x=cn, fill=col)) +
            ggplot2::geom_bar() + ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks=seq(input$cnMin,input$cnMax,1), labels=c(seq(input$cnMin,input$cnMax-1,1), paste0(">",input$cnMax)), limits=c(input$cnMin,input$cnMax+.2)) + 
              ggplot2::theme(legend.position=c(1,1),legend.justification=c(1,1)) + extra.gg + ggplot2::xlab("Copy Number estimates") + ggplot2::ylab("number of calls")
          if(input$col=="sample"){
            gp = gp  + ggplot2::guides(fill=FALSE)
          }
          gp
        })
        output$freq = shiny::renderPlot({
          f.df = freq.df()
          f.df = dplyr::summarize(dplyr::group_by(f.df, chr, start, end), nb=sum(nb), prop=sum(prop), gen.kb=head((end-start)/1e3, 1))
          if(input$freq.rep=="nb"){
            ggp = ggplot2::ggplot(dplyr::arrange(f.df[which(f.df$nb>=input$nbMin),], chr), ggplot2::aes(x=nb,  y=gen.kb, fill=chr)) + ggplot2::xlab("number of samples") 
          } else {
            ggp = ggplot2::ggplot(dplyr::arrange(f.df[which(f.df$nb>=input$nbMin),], chr), ggplot2::aes(x=prop, y=gen.kb, fill=chr)) + ggplot2::xlab("proportion of samples") 
          }
          ggp + ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() +
            ggplot2::ylab("abnormal genome (Kb)") +
              ggplot2::guides(fill=FALSE) +
                ggplot2::scale_fill_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),3)) 
        })
        output$freq.chr = shiny::renderPlot({
          if(input$chr=="all"){
            chr.df = freq.df()
            chr.df$chr = factor(chr.df$chr, levels=1:22)
            pdf = plot.df()
            pdf$chr = factor(pdf$chr, levels=1:22)
            facet.o = ggplot2::facet_wrap(~chr,scales="free")
          } else {
            chr.df = freq.df()
            chr.df = chr.df[which(chr.df$chr==input$chr),]
            pdf = plot.df()
            pdf = pdf[which(pdf$chr==input$chr),]
            facet.o = NULL
          }
          pdf$type = ifelse(pdf$fc>1, "duplication","deletion")
          pdf$sample = factor(pdf$sample)
          if(input$fchr.rep=="Stacked"){
            if(input$freq.rep=="nb"){
              ggp = ggplot2::ggplot(chr.df, ggplot2::aes(xmin=start/1e6, xmax=end/1e6, ymin=0, ymax=nb, fill=type)) + ggplot2::ylab("number of samples")  +
                ggplot2::geom_segment(ggplot2::aes(x=(end+start)/2e6,xend=(end+start)/2e6, y=0, yend=nb, colour=type), alpha=.2, data=chr.df[which(chr.df$end-chr.df$start<max(chr.df$end/1e3)),])
            } else {
              ggp = ggplot2::ggplot(chr.df, ggplot2::aes(xmin=start/1e6, xmax=end/1e6, ymin=0, ymax=prop, fill=type)) + ggplot2::ylab("proportion of samples") +
                ggplot2::geom_segment(ggplot2::aes(x=(end+start)/2e6,xend=(end+start)/2e6, y=0, yend=prop, colour=type), alpha=.2, data=chr.df[which(chr.df$end-chr.df$start<max(chr.df$end/1e3)),])
            }
            return(ggp +
                   ggplot2::scale_fill_brewer(palette="Set1") + 
                   ggplot2::scale_x_continuous(breaks=seq(0,max(chr.df$end)/1e6,20)) + 
                   ggplot2::geom_rect(alpha=.7) + ggplot2::theme_bw() +
                   ggplot2::theme(legend.position="top") + 
                   ggplot2::xlab("position (Mb)") + facet.o)
          } else {
            widths = pdf$end-pdf$start
            wmin = max(chr.df$end/1e3)
            if(any(widths<wmin))widths[widths<wmin] = wmin
            pdf$end = pdf$start + widths
            return(ggplot2::ggplot(pdf, ggplot2::aes(xmin=start/1e6, xmax=end/1e6, ymin=as.numeric(sample)-.5, ymax=as.numeric(sample)+.5, fill=type)) +
                   ggplot2::scale_fill_brewer(palette="Set1") +
                   ggplot2::scale_x_continuous(breaks=seq(0,max(chr.df$end)/1e6,20)) + 
                   ggplot2::geom_rect() + ggplot2::theme_bw() +
                   ggplot2::ylab("sample") + ggplot2::theme(axis.text.y=ggplot2::element_blank(),legend.position="top") + 
                   ggplot2::xlab("position (Mb)") + facet.o)
          }
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
          res.df[which(res.df$qv<as.numeric(input$fdr)),]
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
          pdf = res.m()
          ggplot2::ggplot(pdf[which(pdf$nb.bin.cons>input$nbc),], ggplot2::aes(x=fc*2)) +
            ggplot2::geom_bar() + ggplot2::theme_bw() + ggplot2::xlim(0,5)
        })
        output$freq = shiny::renderPlot({
          p.df = plot.df()
          freq.df = dplyr::summarize(dplyr::group_by(p.df, chr, start),
            freq.n = length(start), freq.p = length(start)/nb.samp)
          ggplot2::ggplot(freq.df, ggplot2::aes(x=freq.n)) +
            ggplot2::geom_bar() + ggplot2::theme_bw()
        })
      }))
  }
}
