### About ###
This app can be used to evaluate the results of *Phylobayes* runs in an interactive and graphical way. [*Phylobayes*](http://www.atgc-montpellier.fr/phylobayes/) is a Bayesian phylogenetics package intended mainly for protein sequence analyses. Although it comes with built-in methods for assessing convergence and mixing of parameters (`tracecomp` & `bpcomp`), additional analyses are often helpful to evaluate the outcomes of a *Phylobayes* run. *postpb* bundles a number of common post-analyses and exploratory plotting approaches for *Phylobayes* in one shiny app.
<br>

### Installation ###
To run *postpb* locally on your machine, open `R` and follow these instructions: 

```R
# Install required libraries (skip if already installed)
install.packages(c( "ape",
                    "coda",
                    "distory",
                    "dplyr",
                    "DT",
                    "ggplot2",
                    "gridExtra",
                    "markdown",
                    "phytools",
                    "shiny",
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

A standalone version of the app with reduced functionality can be found on [github](https://github.com/gerthmicha/pbplot/README_rscript.md).
<br>

### Usage ###
There are two modules in this app, accessible through the tabs in the navbar on top of this page: **Traces** and **Trees**. Global options can be found on the sidebar on the left, and options specific for each analysis type will be in the main window. 

#### Traces ####
Here, Phylobayes trace files can be summarized (file name ending in `.trace`). Upload 1 or 2 trace files from Phylobayes chains that were run under identical options. Click through the tabs to access the different funtions. Use the sliders in the sidebar to adjust the options. Currently implemented functions:
  
  - Plots of trace parameters (XY, violin, density) 
  - interactive burnin adjustment
  - adjustment of plot appearance and size
  - summary statistics for assessment of convergence, stationarity, and mixing
  - pdf export

#### Trees ####
Here, tree samples obtained with Phylobayes can be analysed (file name ending in `.treelist`). Typically, one would upload 2 tree sample files from Phylobayes chains run under identical options, although loading a single tree sample file is also possible. Click through the tabs to access the different funtions. Check the sidebar for options. Currently implemented functions:

  - plotting of each tree sample for both chains side by side
  - consensus tree including branch lengths and posterior probabilities
  - consensus export in newick file format
  - highlighting conflicting nodes between the consensus trees of 2 chains
  - Robinson-Foulds distances of trees within and between chains over run iterations
  - posterior probabilities for any taxonomic group
  - basic tree manipulations: re-rooting, taxon coloring
  - adjustment of plot appearance and size
  - pdf export
  
### Example ###
tba
<br>
<br>

### Notes ###
Most of this app is optimized for usage with 2 *Phylobayes* chains run under identical models with the same dataset, as this is the most typical usage. Support for analysing more chains is planned.
<br>
<br>
For questions & comments please raise an issue at [github](https://github.com/gerthmicha/postpb/issues)
