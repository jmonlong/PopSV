##' Interactive tool to fine tune breakpoint location from the BAM files. An interactive Shiny app will open in your web browser. There, the bp coverage of the `test.sample` will be shown along with the coverage of the reference samples (`ref.samples`). The user can then slide the breakpoint where the tested samples starts to deviate from the reference. If the deviation is not clear a few parameters (e.g. mapping quality, resolution) can be changed live. Once satisfied the user pushes a button and the function returns the new breakpoint locations. In practice this function can be looped over to quickly fine-tune the breakpoints of several CNVs.
##'
##' Access to the BAM files is necessary, but also to a web-browser (for the application). Hence I would recommend to mount the cluster/server with the BAM on your personnal laptop to use this. Then you might need to tweak a bit the paths in the 'bam' column of 'files.df' to include your mounting location.
##' @title Interactive breakpoint finder
##' @param chr the chromosome
##' @param start the start position
##' @param end the end position
##' @param test.sample the sample with abnormal coverage
##' @param files.df a data.frame with the path to the bam files
##' @param ref.samples other samples to visualize
##' @param proper should proper mapping be counted.
##' @param nb.cores number of cores
##' @param bp.res default resolution (bp)
##' @param flanks default flank size (bp)
##' @param related.samples samples related to the `test.sample` (e.g. parents) to be drawn in a special color.
##' @return a list with
##' \item{bk.df}{a data.frame with the breakpoint coordinates}
##' \item{graph}{the ggplot object with the final graph}
##' @author Jean Monlong
##' @export
breakpoint.finder.interactive <- function(chr,start,end, test.sample, files.df, ref.samples, proper = TRUE, nb.cores=1, bp.res=1, flanks=2000, related.samples=NULL) {

  if(!all(c(test.sample, ref.samples, related.samples) %in% files.df$sample)){
    stop("One of the 'test.sample' or 'ref.samples' is not in 'files.df'.")
  }
  
  gr.bam <- function(bam.file, gr, proper = TRUE, map.quality = 30) {
    bai.file = sub("bam$", "bai", bam.file, perl = TRUE)
    if (!file.exists(bai.file)) {
      bai.file = paste0(bam.file, ".bai")
      if (!file.exists(bai.file)) {
        stop("Index file is missing (neither '.bai' nor '.bam.bai').")
      }
    }
    param = Rsamtools::ScanBamParam(which = gr, what = c("rname", "pos", "qwidth", "mapq"), flag = Rsamtools::scanBamFlag(isProperPair = proper, isDuplicate = FALSE, isNotPassingQualityControls = FALSE, isUnmappedQuery = FALSE))
    bam = Rsamtools::scanBam(bam.file, index = bai.file, param = param)
    bam = bam[which(unlist(lapply(bam, function(ee)length(ee[[1]])>0)))]
    bam.df = do.call(rbind, lapply(names(bam),function(bin) data.frame(bin=bin, as.data.frame(bam[[bin]]))))
    bam.df = bam.df[which(bam.df$mapq >= map.quality),]
    if(is.null(bam.df)) return(GenomicRanges::GRanges())
    return(with(bam.df, GenomicRanges::GRanges(rname, IRanges::IRanges(pos, width = qwidth), bin = bin)))
  }
  cov.reads <- function(reads.gr, win.gr){
    GenomicRanges::countOverlaps(win.gr, reads.gr)
  }

  shiny::runApp(list(

    ui = shiny::fluidPage(
      shiny::headerPanel("PopSV - Breakpoint finder"),
      shiny::sidebarPanel(
        shiny::numericInput("flanks", "Flanks", flanks, step=500, min=0),
        shiny::numericInput("bp.res", "Resolution (bp)", bp.res, step=1, min=0),
        shiny::numericInput("map.quality", "Minimum mapping quality", 30, step=10, min=0),
        shiny::hr(),
        shiny::textOutput("stats"),
        shiny::hr(),
        shiny::textInput("comment", "Comment:"),
        shiny::actionButton("exp","Done"), shiny::textOutput("export")),
      shiny::mainPanel(shiny::plotOutput("cov"),
                       shiny::wellPanel(shiny::sliderInput("start", "Start", start-2*flanks,end+2*flanks, value=start),
                                        shiny::sliderInput("end", "End", start-2*flanks,end+2*flanks, value=end)))),

    server = function(input, output) {

      reads.comp <- shiny::reactive({
        bkpts = c(start, end)
        st.fl = start-input$flanks
        end.fl = end+input$flanks
        mapq = input$map.quality
        gr.f = GenomicRanges::GRanges(chr, IRanges::IRanges(st.fl, end.fl))
        list(gr.f=gr.f, reads.l=parallel::mclapply(c(ref.samples, related.samples,test.sample), function(samp.i){
                                       gr.bam(files.df$bam[files.df$sample==samp.i], gr.f, proper=proper, map.quality=mapq)
                                     }, mc.cores=nb.cores))
      })

      gp.flanks <- shiny::reactive({
        reads.out = reads.comp()
        gr.f = reads.out$gr.f
        reads.l = reads.out$reads.l
        gr.bk = seq(GenomicRanges::start(gr.f), GenomicRanges::end(gr.f), input$bp.res)
        gr.frag = GenomicRanges::GRanges(chr, IRanges::IRanges(gr.bk[-length(gr.bk)], width = input$bp.res))
        cov.l = parallel::mclapply(1:length(reads.l), function(ii) {
          cov.reads(reads.l[[ii]], gr.frag)
        }, mc.cores=nb.cores)
        cov = matrix(0, length(gr.frag), length(cov.l))
        for (ii in 1:length(cov.l)) cov[, ii] = cov.l[[ii]]
        norm.f = apply(cov[c(1:floor(input$flanks/input$bp.res/2),(nrow(cov)-floor(input$flanks/input$bp.res/2)+1):nrow(cov)),],2,mean, na.rm=TRUE)
        cov = cov %*% diag(mean(norm.f)/norm.f)
        colnames(cov) = c(ref.samples, related.samples, test.sample)
        rownames(cov) = GenomicRanges::start(gr.frag)
        cov.df = reshape::melt.array(cov)
        colnames(cov.df) = c("position","sample","cov")
        cov.df$status = ifelse(cov.df$sample==test.sample, "tested","reference")
        cov.df$sample = factor(cov.df$sample, levels=c(ref.samples, related.samples,test.sample))
        if(!is.null(related.samples)) cov.df$status[which(cov.df$sample%in% related.samples)] = "related"
        cov.df$status = factor(cov.df$status, levels=c("reference","related","tested"))
        status = position = NULL ## Uglily appease R checks
        gp = ggplot2::ggplot(cov.df, ggplot2::aes(x=position,y=cov, group=sample)) + ggplot2::geom_line(ggplot2::aes(alpha=status, size=status,col=status)) + ggplot2::theme_bw() + ggplot2::scale_size_manual(values=1:nlevels(cov.df$status)) + ggplot2::scale_alpha_manual(values=seq(.5,1,length.out=nlevels(cov.df$status))) + ggplot2::theme(legend.position="top")
      })

      output$cov = shiny::renderPlot({
        pdf = gp.flanks()
        pdf + ggplot2::geom_vline(xintercept=c(input$start, input$end),linetype=2)
      })

      output$stats = shiny::renderText({
        paste0("Size: ", input$end-input$start+1," bp")
      })

      output$export = shiny::renderText({
        if(input$exp){
          pdf = gp.flanks()
          pdf = pdf + ggplot2::geom_vline(xintercept=c(input$start, input$end),linetype=2)
          shiny::stopApp(list(bk.df=data.frame(chr=chr,start=input$start, end=input$end, comment=input$comment),
                              graph=pdf))
        }
        ""
      })

    }))
}
