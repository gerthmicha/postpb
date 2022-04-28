# load libraries
packages <- c("shiny", "colourpicker", "coda", "data.table", "distory", 
              "doParallel", "DT", "dplyr", "ggplot2", "foreach", "gridExtra", 
              "markdown", "phytools", "shinycssloaders", "shinyjs", "shinythemes",
              "shinyWidgets", "stringr", "tidyr")
invisible(lapply(packages, library, character.only = TRUE))