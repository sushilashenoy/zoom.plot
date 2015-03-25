# zoom.plot R package

##Installation
Do this:
```{r Install icreport, eval = FALSE}
devtools::install_github("sushilashenoy/zoom.plot")
```
If that didn't work it may be because you don't have the devtools package installed, in which case, do this:
```{r Install icreport, eval = FALSE}
install.packages("devtools") 
devtools::install_github("sushilashenoy/zoom.plot")
```
If that didn't work I'm sorry.


##About
This R package is a bunch of functions for making plots that I find useful. Originally I made it for "zoom-in" plots which in this case are association p-values for SNPs in a genetic region (a few Mb) with the genes plotted below that and LD (linkage disequilibrium) plottted below that. It can still do that plus a few random other things (like cosine interpolation) that are very poorly documented (but hopefully getting better).
