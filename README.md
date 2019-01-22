### About ###
This app can be used to evaluate the results of Bayesian phylogentics software packages (such as [*Phylobayes*](http://www.atgc-montpellier.fr/phylobayes/) and [*MrBayes*](http://nbisweden.github.io/MrBayes/)) in an interactive and graphical way. Although both programmes come with built-in methods for assessing convergence and mixing of parameters, additional analyses are often helpful when interpreting the outcomes of Bayesian phylogenetics software. *postpb* bundles a number of common post-analyses and exploratory plotting approaches in one shiny app.


### Notes ###
This app is optimized and tested for usage with *Phylobayes* and *MrBayes* output files, but should work for any other Bayesian phylogenetic software that produces tabulated parameter files and tree files in Newick or Nexus format. It currently supports analysis of up to 4 chains, which should have been run under identical models and the same dataset. 


### Installation ###
To run *postpb* locally on your machine, open `R` and follow these instructions: 

```R
# Install required libraries (skip if already installed)
install.packages(c( "distory",
                    "dplyr",
                    "DT",
                    "gridExtra",
                    "markdown",
                    "phytools",
                    "shinycssloaders",
                    "shinyjs",
                    "shinythemes",
                    "shinyWidgets"), 
                    dependencies = TRUE)

# Run Shiny app
library(shiny)
shiny::runGitHub("postpb","gerthmicha")
```

The app is also hosted on [shinyapps.io](https://gerthmicha.shinyapps.io/postpb/) and can be accessed without installation. This is only recommended for testing though, as it will be slow for bigger files. There is also usage limit of 25h per month and if this is used up the app will be offline.

A commandline version of the app with reduced functionality can be found on [github](https://github.com/gerthmicha/pbplot/README_rscript.md).
<br>


### Usage ###
There are two modules in this app, accessible through the tabs in the navbar on top of this page: **Traces** and **Trees**. Global options can be found on the sidebar on the left, and options specific for each analysis type will be in the main window. 


#### Traces ####
Here, parameter files from Phylobayes or MrBayes can be summarized (file name ending in `.trace` or `.p`). Upload 1–4 of trace files from mcmc chains that were run under identical options. Click through the tabs to access the different funtions. Use the sliders in the sidebar to adjust the options. Currently implemented functions:

**Global options for all tabs (in sidebar)**

* Upload and select which trace files to analyse
* Thinning interval [default: consider every 10th generation]
* burnin adjustment [default: 20% of shortest chain]
* Plot adjustments (size of plot area, size of lines & points)
* Download of all plots in single PDF

**XY tab**

* XY plots of all parameters over mcmc generation

**Violin tab**

* Violin plots of all parameters

**Density tab** 

* Density plots of all parameters

**Summary statistics tab**

* Various statistics for assessing convergence, stationarity, and mixing of chains
* Export in several formats


#### Trees ####
Here, tree samples obtained with Phylobayes or MrBayes can be analysed (file name ending in `.treelist` or `.t`). Typically, one would upload 2–4 tree sample files from mcmc chains run under identical options, although loading a single tree sample file is also possible. Click through the tabs to access the different funtions. Check the sidebar for options. Currently implemented functions:

**Global options for all tabs (in sidebar)**

* Upload and select which tree files to analyse
* Thinning interval [default: consider every 10th generation]
* Burnin adjustment [default: 20% of shortest chain]
* Select outgroup(s) for rooting, reroot with outgroup, midpoint, or unroot
* Highlighting taxa (three different colors)
* Formating of tree labels
* Plot adjustments (size of plot area, size of lines & points)
* Download of currently displayed tree plot

**Trees tab**

* Plot trees for chains side by side, browse through mcmc generations (can be animated) 

**Consensus tab**

* Consensus tree including branch lengths and posterior probabilities (pp)
* Collapse branches below interactively chosen pp
* Export in newick file format

**Difference**
* Highlight conflicting nodes between consensus trees of all chains (requires at least 2 tree files)
 
**Pairwise Robinson-Foulds**
* Plot Robinson-Foulds distances of trees within and between chains over generation

**Bipartition support** 
* Posterior probabilities for any taxonomic group in posterior sample of trees. 

  
### Example ###
tba
<br>
<br>

### Citation ###
If you use postpb in your research, I would if you could cite it as 

* Michael Gerth (2019) Postpb. Available under https://github.com/gerthmicha/postpb

A preprint describing its functionality is in preparation. Please also cite the appropriate R packages listed below: 

**Tree plots, RF and consensus calculations:** 

  * Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289–290. https://cran.r-project.org/package=ape
  
  * Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223. https://cran.r-project.org/package=phytools
  
  * Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating Phylogenetic and Heirarchical Clustering Trees. arXiv:1006.1015v1. https://cran.r-project.org/package=distory

**MCMC diagnostics:**

  * Plummer, M., Best, N., Cowles, K. and Vines, K. (2006) CODA: convergence diagnosis and output analysis for MCMC. R News, 6, 7–11. https://cran.r-project.org/package=coda
  

<br>
<br>
### Questions, comments, bugs ###
Please raise an issue at [github](https://github.com/gerthmicha/postpb/issues)
<br>
