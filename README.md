#### pbplot

This app can be used to plot trace files from Phylobayes runs which is helpful in assessing convergence and mixing of parameters. Input must be at least one, and at most four Phylobayes trace files from chains run with identical parameters.

The app is hosted on shinyapps (limited to 25h of usage/month) https://gerthmicha.shinyapps.io/pbplot/ or can be run locally: 

```R
# Install required libraries
install.packages(c("ggplot2", "gridExtra", "LaplacesDemon", "markdown", "shiny", "shinythemes"), dependencies = TRUE)

# Run Shiny app
library(shiny)
shiny::runGitHub("pbplot","gerthmicha")
```
<br>


#### pbpplot.R

This is the standalone version to be run as executable R script

```
Usage: /usr/bin/pbplot.R [options]

Options:
	-b NUMBER, --burnin=NUMBER
		number of burnin iterations to discard [required]

	-f, --file
		Shall plot be saved to pdf file instead of disp? [default: FALSE]

	-h, --help
		Show this help message and exit

NOTE:	This script plots the parameters estimated during a phylobayes run.
	It should be run in a directory containing 1–4 trace files created by phylobayes.
	Requires R packages ggplot2, optparse, and gridExtra.
```
