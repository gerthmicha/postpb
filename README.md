## About 
*postpb* can be used to evaluate the results of Bayesian phylogentics software packages (such as [*Phylobayes*](http://www.atgc-montpellier.fr/phylobayes/) and [*MrBayes*](http://nbisweden.github.io/MrBayes/)) in an interactive and graphical way. Although both programmes come with built-in methods for assessing convergence and mixing of parameters, additional analyses are often helpful when interpreting the outcomes of Bayesian phylogenetics software. *postpb* bundles a number of common post-analyses and exploratory plotting approaches in one shiny app.


## Notes
*postpb* is optimized and tested for usage with *Phylobayes* and *MrBayes* output files, but should work for any other Bayesian phylogenetic software that produces tabulated parameter files and tree files in Newick or Nexus format. It currently supports analysis of up to 4 chains, which should have been run under identical models and the same dataset. 


## Installation
To run *postpb* locally on your machine, open `R` and follow these instructions: 

```R
# Install packages that aren't installed already
packages <- c( "colourpicker",
                "data.table",
                "distory",
                "doParallel",
                "dplyr",
                "DT",
                "foreach",
                "gridExtra",
                "markdown",
                "phytools",
                "shinycssloaders",
                "shinyjs",
                "shinythemes",
                "shinyWidgets")
                
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
                    
# If you want to use the rwty functionality, also install this package
install.packages("rwty", dependencies = TRUE)

# Run Shiny app
shiny::runGitHub("postpb","gerthmicha")
```

*postpb* is also hosted on [shinyapps.io](https://gerthmicha.shinyapps.io/postpb/) and can be accessed without installation. This is only recommended for testing though, as it will be slow for bigger files. There is also usage limit of 25h per month and if this is used up the app will be offline.

A command line version of the app with reduced functionality can be found on [github](https://github.com/gerthmicha/pbplot/README_rscript.md).
<br>

## Quick start 
There are two modules in *postpb*, accessible through the tabs in the navbar on top of this page: **Traces** and **Trees**. Global options can be found on the sidebar on the left, and options specific for each analysis type will be in the main window. Usage is fairly self-explanatory. For a quick start, load one of the examples (see below), and click through the interface to explore the different plots and statistics.   

## Tutorial

A detailed turorial can be found in the wiki. 

## Examples
*postpb* comes with two sets of example data that can be loaded by clicking the buttons under the file selection dialogs in the trace and tree tabs. 

**Example 1** is taken from [Kocot et al. (2017)](https://doi.org/10.1093/sysbio/syw079), and it contains the first 5,000 iterations of their *Phylobayes* analysis of the "LB_106" matrix under the CAT+GTR model. The example contains traces and trees for 4 chains. Please refer to the paper for details on how the dataset was compiled. All data associated with this paper is accessible under https://doi.org/10.5061/dryad.30k4v.2.

**Example 2** is a *MrBayes* analysis of the "primates" alignment that is distributed with *MrBayes* (https://github.com/NBISweden/MrBayes/tree/develop/examples) and is originally from [Hayasaka et al. (1988)](https://doi.org/10.1093/oxfordjournals.molbev.a040524). The analysis was run under a HKY+G model using 2x4 chains for 1,000,000 generations with a sampling frequency of 50. 

## Citation 
If you use *postpb* in your research, I would appreciate if you could cite it as 

* Michael Gerth (2022) Postpb. Available under https://github.com/gerthmicha/postpb

A preprint describing *postpb* is in preparation. Please also cite the appropriate R packages listed below: 

**Tree plots, RF and consensus calculations:** 

  * Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289–290. https://cran.r-project.org/package=ape
  
  * Revell LJ (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol Evol 3: 217-223. https://cran.r-project.org/package=phytools
  
**MCMC diagnostics:**

  * Plummer M, Best N, Cowles K, Vines K (2006) CODA: convergence diagnosis and output analysis for MCMC. R News 6: 7–11. https://cran.r-project.org/package=coda
  
**RWTY:**
  * Warren DL, Anthony JG, Lanfear RL (2017) "RWTY (R We There Yet): an R package for examining convergence of Bayesian phylogenetic analyses." Mol Biol Evol 34: 1016-1020. https://cran.r-project.org/package=rytn
  
 
## Questions, comments, bugs

Please raise an issue at [github](https://github.com/gerthmicha/postpb/issues).

<br>

