### Motivation ###
This app can be used to evaluate the results of Phylobayes runs in an interactive and graphical way. [Phylobayes](http://www.atgc-montpellier.fr/phylobayes/) is a Bayesian phylogenetics package intended mainly for protein sequence analyses. Although it comes with built-in methods for assessing convergence and mixing of parameters (`tracecomp` & `bpcomp`), plotting of parameters and trees is often very helpful to evaluate the outcomes of a Phylobayes run. This is what this app is meant to do.
<br>

### Installation / access ###
The app is hosted on shinyapps.io: https://gerthmicha.shinyapps.io/pbplot/ and can be accessed without installation. There is a usage limit of 25h per month however and if this is used up the app will be offline. To run it locally on your machine, open `R` and follow these instructions: 
```R
# Install required libraries
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
shiny::runGitHub("pbplot","gerthmicha")
```
A standalone version of the app with reduced functionality can be found on [github](https://github.com/gerthmicha/pbplot/README_rscript.md).
<br>

### Usage ###
There are two modules in this app, accessible through the tabs in the navbar on top of this page: **Traces** and **Trees**

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
  - highlighting of conflicting nodes between the 2 consensus trees
  - Robinson-Foulds distances of trees within and between chains
  - basic tree manipulations: re-rooting, taxon coloring
  - adjustment of plot appearance and size
  - pdf export

  
### Example ###
tba
<br>
<br>
For questions & comments please raise an issue at [github](https://github.com/gerthmicha/pbplot/issues)
