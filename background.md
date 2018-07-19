<br><br>
This app can be used to plot trace files from Phylobayes runs which is helpful in assessing convergence and mixing of parameters. Input must be at least one, and at most four Phylobayes trace files from chains run with identical parameters.

The app can be accessed through

* https://gerthmicha.shinyapps.io/pbplot/ (Shinyapps) or 
* https://github.com/gerthmicha/pbplot (R script to be run in the command line)

and can be run as Shiny app locally: 

```R
# Install required libraries
install.packages(c("ggplot2", "gridExtra", "LaplacesDemon", "markdown", "shiny"), dependencies = TRUE)

# Run Shiny app
library(shiny)
shiny::runGitHub("pbplot","gerthmicha")
```
<br><br>
Questions & comments welcome: gerth@liv.ac.uk 