# scDD: An R package to identify differentially distributed genes in scRNA-seq 

## Installation 
The Bioconductor landing page is 
[https://bioconductor.org/packages/scDD](https://bioconductor.org/packages/scDD).
To install, make sure you have the current version of Bioconductor, and 
use the following commands:

```R
install.packages("BiocManager")
BiocManager::install("scDD")
```


## Quick Start
For examples and tips on using the package, please see the vignette PDF 
[here](http://bioconductor.org/packages/devel/bioc/vignettes/scDD/inst/doc/scDD.pdf)
which you can alternatively bring up by typing 

```R
browseVignettes("scDD")
```

after installing and loading the package into R.  

## Getting Help
Please send bug reports and feature requests by sending a message to Bioconductor
Support at [https://support.bioconductor.org](https://support.bioconductor.org)
or by opening a new issue on 
[this page](https://github.com/kdkorthauer/scDD/issues)

## Citation
Please cite the following publication if you use scDD in your work:
> Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
Kendziorski C. A statistical approach for identifying differential 
distributions in single-cell RNA-seq experiments. Genome Biology. 
2016 Oct 25;17(1):222. 
[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y)

You can obtain this citation within R by typing

```R
citation("scDD")
```
