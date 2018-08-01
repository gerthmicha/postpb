### Motivation ###
This app can be used to evaluate the results of Phylobayes runs in an interactive and graphical way. [Phylobayes](http://www.atgc-montpellier.fr/phylobayes/) is a Bayesian phylogenetics package intended mainly for protein sequence analyses. Although it comes with built-in methods for assessing convergence and mixing of parameters (`tracecomp` & `bpcomp`), plotting of parameters and trees is often very helpful to evaluate the outcomes of a Phylobayes run. This is what this app is meant to do.
<br>

### Installation / access ###
The app is hosted on shinyapps.io: https://gerthmicha.shinyapps.io/pbplot/ and can be accessed without installation. There is a usage limit of 25h per month however and if this is used up the app will be offline. To run it locally on your machine, open `R` and follow these instructions: 
```R
# Install required libraries
install.packages(c("ape", "distory", "ggplot2", "gridExtra", "LaplacesDemon", "markdown","phytools", "shiny", "shinythemes"), dependencies = TRUE)

# Run Shiny app
library(shiny)
shiny::runGitHub("pbplot","gerthmicha")
```
A standalone version of the app with reduced functionality can be found on github: https://github.com/gerthmicha/pbplot.
<br>

### Usage ###
There are two modules in this app, accessible through the tabs in the navbar on top of this page: **Traces** and **Trees**

#### Traces ####
Here, Phylobayes trace files can be summarized (file name ending in `.trace`). Upload 1–4 trace files from Phylobayes chains that were run under identical options. Click through the tabs to access the different funtions. Use the sliders in the sidebar to adjust the options. Currently implemented functions:
  
  - plotting of trace parameters as dot and/or line plot 
  - plotting of the distribution of trace parameters
  - interactive adjustment of the number of burnin generations
  - interactive adjustment of plot appearance and size
  - export of plots as pdf
  - summary statistics of runs (similar to `tracecomp` fundtion of Phylobayes)

#### Trees ####
Here, tree samples obtained with Phylobayes can be analysed (file name ending in `.treelist`). Typically, one would upload 2 tree sample files from Phylobayes chains run under identical options, although loading a single tree sample file is also possible. Click through the tabs to access the different funtions. Check the sidebar for options. Currently implemented functions:

  - plotting of each tree sample for both chains side by side
  - consensus trees for each of the chain, including branch lengths and posterior probabilities
  - highlighting of conflicting nodes between the 2 consensus trees
  - interactive adjustment of the number of burnin generations
  - interactive rooting for all trees
  - interactive coloring of taxa 
  
### Example ###
tba
<br>
<br>
For questions & comments please raise an issue at github: https://github.com/gerthmicha/pbplot