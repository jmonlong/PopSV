---
layout: page
title: Visualization
permalink: /3-Visualization.md/
---

Soon.

## Data quality before analysis

## Summary of the calls

## Example of a region

## Interactive breakpoint fine-tuning

PopSV calls are defined at the bin resolution, hence hundreds to thousands of bp. In order to fine-tune the breakpoints location we implemented an interactive app. `breakpoint.finder.interactive` takes as input a genomic region, some samples and the path to the bam files. It then opens an application in the web-browser where the user can slide the breakpoints where desired. The interface looks like this:

![Interactive breakpoint fine-tuning](public/bkptInteractive.jpg)


```r
breakpoint.finder.interactive(chr.i, start.i, end.i, test.sample=test.samp, files.df=files.df, ref.samples=ref.samp)
```

Once the user have decided on breakpoint location, he pushes the *Done* button and the function returns the final breakpoint location as well as the final graph. In practice you might want to save this output, and maybe even loop across several regions:

```r
## Interactive breakpoint fine-tuning
bk.l = NULL
while(length(bk.l)<nrow(val.df)){
  ii = length(bk.l) + 1
  bk.l = c(bk.l, list(breakpoint.finder.interactive(val.df$chr[ii], val.df$start[ii], val.df$end[ii], val.sample[ii], files.df, ref.samples=controls.samp)))
}
## Merge results
bk.df = do.call(rbind, lapply(bk.l, function(ll)ll$bk.df))
## Final graph for each region
pdf("bkptResults.pdf",14,5)
lapply(bk.l, function(ll)ll$graph)
dev.off()
```

You can also specify related samples (e.g. parents) that will be displayed with a specific color. This might help understanding the mode of inheritance or veracity of the variant. Then it would look like this:

![Interactive breakpoint fine-tuning trios](public/bkptInteractive-trios.jpg)
