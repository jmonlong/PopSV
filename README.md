PopSV
=====
PopSV is a structural variation (SV) detection method from high-throughput sequencing. 
Abnormal Read-Depth signal is detected by using a population of samples ase reference. Thanks to this population
view the whole genome can be robustly interrogated, including regions of low mappability. Moreover, any divergence from
the reference samples are detected, even if the signal is incomplete, e.g. tumoral aberrations or SV involving repeats.

To install the latest development version: `devtools::install_github("jmonlong/PopSV")`. This requires `devtools` package (more information [here](https://github.com/hadley/devtools) which can be easily installed with `install.packages("devtools")`. 
