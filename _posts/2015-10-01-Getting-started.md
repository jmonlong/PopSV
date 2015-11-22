---
layout: post
title: Getting started
---

### Installation
To install the latest development version:

```r
devtools::install_github("jmonlong/PopSV")
```

This command requires [*devtools* package](https://github.com/hadley/devtools) which can be easily installed with :

```r
install.packages("devtools")
``` 

Some [Bioconductor](http://bioconductor.org/) packages are also necessary and not installed automatically. Running the following commands should be sufficient :

```r
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19", "Rsamtools")
```

**R 3.1 or higher** is required.

### Usage

PopSV package is used as any R package, by simply loading the library and using the provided funtions.

```r
library(PopSV)
...
```
