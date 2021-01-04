# all required libraries can be installed with 
# install.packages(c( "colourpicker",
#                     "distory",
#                     "data.table",
#                     "doParallel",
#                     "dplyr",
#                     "DT",
#                     "foreach",
#                     "gridExtra",
#                     "markdown",
#                     "phytools",
#                     "shinycssloaders",
#                     "shinyjs",
#                     "shinythemes",
#                     "shinyWidgets"), 
#                  dependencies = TRUE)


# load libraries
library(shiny, warn.conflicts = FALSE)
library(colourpicker, warn.conflicts = FALSE)
library(coda, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(distory, warn.conflicts = FALSE)
library(doParallel, warn.conflicts = FALSE)
library(DT, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(foreach, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)
library(markdown, warn.conflicts = FALSE)
library(phytools, warn.conflicts = FALSE)
library(shinycssloaders, warn.conflicts = FALSE)
library(shinyjs, warn.conflicts = FALSE, quietly = TRUE)
library(shinythemes, warn.conflicts = FALSE)
library(shinyWidgets, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)