# scDD
Version 0.99.0 of R Package 'scDD' (for submission to Bioconductor).  To access 
previous versions of the package (including 1.1.0 as used in the [Genome 
Biology publication]
(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y), 
see the 
[Releases page](https://github.com/kdkorthauer/scDD/releases)).

To Install:

```
install.packages("devtools")

devtools::install_github("kdkorthauer/scDD")
```

For examples and tips on using the package, please see the vignette PDF, 
which you can bring up by typing 

```
vignette("scDD"")
```

after installing and loading the package into R.  You can also preview the
text and code in the vignette (without the results it produces) in the raw
.Rnw file without first installing scDD by clicking [here](https://github.com/kdkorthauer/scDD/blob/master/vignettes/scDD.Rnw).

Please cite the following publication if you use scDD in your work:
> Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
Kendziorski C. A statistical approach for identifying differential 
distributions in single-cell RNA-seq experiments. Genome Biology. 
2016 Oct 25;17(1):222. 
[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y]
(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y)

Note: The package was built using the development version of R (3.4.0,
2016-10-26 r71594)
