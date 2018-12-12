#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


# load all libraries necessary for 
# ui
library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
# & server
library(coda)
library(distory)
library(dplyr)
library(DT)
library(ggplot2)
library(gridExtra)
library(markdown)
library(phytools)
library(tidyr)


# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),  # defines theme of app
          # Layout of app: 1 navbar with 3 tab panels. The first 2 tab panels have 1 sidebar + 1 main display (including multiple tabs each)
          # Final tab panel is markdown help file 
          navbarPage("Phylobayes run stats",
                     # first tab panel for trace analyses
                     tabPanel("Parameters",
                              fluidRow(
                                column(4, wellPanel(  # the following defines the elements of the sidebar
                                  fileInput("tracefile", "Select Phylobayes trace files",
                                            multiple = TRUE,
                                            accept = c(".trace")),
                                  hr(),
                                  uiOutput("burnin"),
                                  uiOutput("whichchain"),
                                  numericInput("prop", "Plot every Nth iteration",
                                               value = 10,
                                               width = "30%"),
                                  sliderInput("cex", "Scaling factor for points and lines", 
                                              min = 0.1, 
                                              max = 10, 
                                              value = 1, 
                                              step = 0.1),
                                  sliderInput("height", "Height of plot in pixels", 
                                              min = 100, 
                                              max = 5000, 
                                              value = 1000, 
                                              step = 100),
                                  sliderInput("width", "Width of plot in pixels", 
                                              min = 100, 
                                              max = 5000, 
                                              value = 1000, 
                                              step = 100),
                                  hr(),
                                  conditionalPanel(condition = 'output.whichchain', 
                                                   downloadButton("downloadPDF", "Download pdf of all trace plots"))
                                )),
                                column(8,   # the following defines the elements of the main display (4 tabs in total)
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Trace", 
                                                            prettyCheckboxGroup(inputId = "traceplotstyle", 
                                                                                label = "", 
                                                                                choices = c("lines", "points"), 
                                                                                selected = "lines", 
                                                                                inline = TRUE, 
                                                                                status = "primary",
                                                                                icon = icon("check")),
                                                            hr(),
                                                            uiOutput("tracePlot.ui")),  
                                                   tabPanel("Violin", 
                                                            prettyCheckboxGroup(inputId = "violinplotstyle", 
                                                                                label = "", 
                                                                                choices = c("boxplot", "points"), 
                                                                                selected = c("boxplot"), 
                                                                                inline = TRUE, 
                                                                                status = "primary",
                                                                                icon = icon("check")),
                                                            hr(),
                                                            uiOutput("violinPlot.ui")),
                                                   tabPanel("Density", 
                                                            uiOutput("densePlot.ui")),
                                                   tabPanel("Summary statistics", 
                                                            withSpinner(DT::dataTableOutput("table", width = "80%"), 
                                                                        color = "#2C4152", size=0.5), 
                                                            conditionalPanel(condition = 'output.table',
                                                                             br(), 
                                                                             useShinyjs(), 
                                                                             actionButton("explanation", "Toggle explanations", style='padding:5px 10px; font-size:90%; background-color:white; color:black'),
                                                                             hidden(div(id="stats", includeMarkdown("stats.md"))))))
                                )
                              )
                     ),
                     # second tab panel for tree analyses, set up as "trace" tab panel 
                     tabPanel("Trees", 
                              fluidRow(
                                column(4, 
                                       wellPanel(
                                         fileInput("treefile", "Select Phylobayes tree files",
                                                   multiple = TRUE,
                                                   accept = ".treelist"),
                                         hr(),
                                         numericInput("treethin", "Consider every Nth tree",
                                                      value = 10,
                                                      width = "30%"),
                                         uiOutput("treegens"),
                                         uiOutput("conburnin"),
                                         uiOutput("outgroups"),
                                         conditionalPanel(condition='output.outgroups',
                                                          actionGroupButtons(c("reroot", "midpoint", "unroot"),
                                                                       c("Reroot with outgroup", "Midpoint rooting", "Unroot"), 
                                                                       status = "primary",
                                                                       size = "sm",
                                                                       fullwidth = TRUE)),
                                         hr(),
                                         uiOutput("highlight"),
                                         uiOutput("highlight2"),
                                         conditionalPanel(condition='output.outgroups',
                                                          helpText(HTML("<b>NOTE: Taxa selected twice will be highlighted in <font color=\"#4daf4a\">green</font></b>"))),
                                         hr(),
                                         sliderInput("treecex", 
                                                     "Scaling factor for labels and lines", 
                                                     min = 0.1, 
                                                     max = 10, 
                                                     value = 1, 
                                                     step =  0.1),
                                         sliderInput("treeheight", 
                                                     "Height of plot in pixels", 
                                                     min = 100, 
                                                     max = 5000, 
                                                     value = 800, 
                                                     step = 100),
                                         sliderInput("treewidth", 
                                                     "Width of plot in pixels", 
                                                     min = 100, 
                                                     max = 5000, 
                                                     value = 1000, 
                                                     step = 100),
                                       conditionalPanel(condition = 'output.conburnin', 
                                                        downloadButton("downloadtreePDF", "Download pdf of current tree plot")))
                                ),
                                column(8, 
                                       tabsetPanel(type = "tabs", 
                                                   tabPanel("Trees", 
                                                            uiOutput("treePlot.ui")),
                                                   tabPanel("Consensus", 
                                                            uiOutput("consensusPlot.ui"),
                                                            conditionalPanel(condition = 'output.consensusPlot', 
                                                                             downloadButton("newick", "Export tree in newick format", style='padding:5px 10px; font-size:90%; background-color:white; color:black'))),
                                                   tabPanel("Difference", 
                                                            uiOutput("differencePlot.ui"),
                                                            br(),
                                                            helpText(HTML("Conflicting nodes between consensus trees will be highlighted in <font color=\"#FF1493\"><b>pink</b></font>."))),
                                                   tabPanel("Pairwise Robinson-Foulds", 
                                                            uiOutput("rfPlot.ui"),
                                                            br(),
                                                            conditionalPanel(condition = 'output.rfPlot', 
                                                                             downloadButton("downloadrfplot", "Download plot as pdf", style='padding:5px 10px; font-size:90%; background-color:white; color:black')),
                                                            br(),
                                                            helpText("Calculation of pairwise RF distances can be computationally intensive. Reduce the number of considered trees to speed this up.")))
                                )
                              )          
                     ),
                     tabPanel("About",
                              includeMarkdown("README.md"))
          ))


# Define server
server <- function(input, output, session) {
  
  ###################
  # GENERAL OPTIONS #
  ###################
  # increase maximum upload size of allocated memory 
  options(shiny.maxRequestSize = 500*1024^2)
  
  # shut down R session when browser window is closed
  session$onSessionEnded(function() {
    stopApp()
  })
  
  
                                                ############################
                                                #>>>>>>>>>>>>><<<<<<<<<<<<<#
#################################################>>>>>>> TRACE TAB <<<<<<<<#################################################
                                                #>>>>>>>>>>>>><<<<<<<<<<<<<#
                                                ############################
  
  #########################
  # TRACE FILE PROPERTIES #
  #########################
  
  # get number of generations from tracefiles
  ngen  <- reactive({
    file1 <- input$tracefile$datapath[1]
    trace1 <- read.table(file1, 
                         sep = "\t", 
                         header = TRUE, 
                         comment.char = "") # this should allow to read in trace files from non-mpi phylobayes runs 
    ngen <- nrow(trace1)
    
    # in case 2 trace files were provided, use lowest number of iterations from these two traces
    if(length(input$tracefile$datapath) > 1){
      file2 <- input$tracefile$datapath[2]
      trace2 <- read.table(file2, 
                           sep = "\t", 
                           header = TRUE, 
                           comment.char = "")
    ngen <- min(c(nrow(trace1), nrow(trace2)))
    }
    ngen
  })
  
  # Get names of chains from file names provided
  chainnames <- reactive({
    # get names of chains from trace file names
    strsplit(input$tracefile$name[1:length(input$tracefile$name)], ".trace", fixed = TRUE)
  })

  # Reading the trace files
  tracedata <- reactive({
    req(input$tracefile)
    # read in trace files (at least 1, try up to 4) 
    file1 <- input$tracefile$datapath[1]
    trace1 <- read.table(file1, 
                         sep = "\t", 
                         header = TRUE, 
                         comment.char = "")
    try(file2 <- input$tracefile$datapath[2], silent = TRUE)
    try(trace2 <- read.table(file2, 
                             sep = "\t", 
                             header = TRUE, 
                             comment.char = ""), 
        silent = TRUE)
    
    # rename first variable to "iter" for consistency (sequential phylobayes version calls this "#ngen")
    colnames(trace1)[1] <- "iter"
    try(colnames(trace2)[1] <- "iter", silent = TRUE)
    
    # remove burnin 
    burnin <- input$burnin
    trace1 <- trace1[burnin+1:nrow(trace1), ]
    
    # subsample according to "Plot every Nth generation" slider in UI
    trace1 <- trace1[seq(from = 0, to = nrow(trace1), by= input$prop),]
    
    # repeat same 2 steps for trace 2
    try(trace2 <- trace2[burnin+1:nrow(trace2), ], silent = TRUE)
    try(trace2 <- trace2[seq(from = 0, to = nrow(trace2), by= input$prop),], silent = TRUE)

    
    # Create novel variable "trace", specifying which trace file the data originates from
    trace1$trace <- chainnames()[[1]][1]
    
    # Do for all trace files
    try(trace2$trace <- chainnames()[[2]][1], silent=T) 
    
    # remove irrelavant (in terms of analyses here) second and third variable ("time" & "topo") from trace file 
    plotDF <- trace1[,c(1,4:ncol(trace1))]
    
    # same for other trace files
    try(plotDF <- rbind(trace1[,c(1,4:ncol(trace1))],trace2[,c(1,4:ncol(trace2))]), silent=T)

    # finally 
    plotDF <- gather(plotDF, variable, value, -trace, -iter, na.rm = TRUE, factor_key = TRUE)
    
    plotDF
  })
  
  # filter tracedata to only plot selected trace files in checkbox (default= select all)
  traceDF <- reactive({
    req(input$tracefile)
    
    if(length(input$tracefile$datapath)==1){
      traceDF <- tracedata()
    }
    
    if(length(input$tracefile$datapath)>1){
      traceDF <- filter(tracedata(), trace %in% input$whichchain)
    }
    
    # This adds prettier error message in case no trace file is selected
    validate(
      need(nrow(traceDF) > 0, "Please select at least one trace file!")
    )
    
    traceDF
  })
  #########################  
  
  #############################
  # UI ELEMENTS FOR TRACE TAB #
  #############################
  
  # display burnin slider using the number of generations read from trace file
  output$burnin <- renderUI({
    req(input$tracefile)
    sliderInput("burnin", "Burnin [# of iterations]:", 
                min = 0, 
                max = ngen(),  
                value = ngen()/5, # default = 20% of iterations
                step = 100)
  })
  
  # display checkbox to select which trace file to plot
  output$whichchain <- renderUI({
    req(input$tracefile$datapath[2])
    names <- lapply(chainnames(), `[[`, 1)
    prettyCheckboxGroup(inputId = "whichchain", 
                        label = "Select trace file[s] to plot", 
                        choices = names, 
                        selected = names, 
                        inline = FALSE, 
                        status = "primary",
                        icon = icon("check"))
  })
  
  # Button that toggles explanations for the statistics
  observeEvent(input$explanation, {
    toggle('stats')
  })   
  
  # Create pdf download handle for plots
  output$downloadPDF <- downloadHandler(filename = "traceplots.pdf",
                                        content = function(file) {
                                          pdf(file,  height=input$height/72, width=input$width/72)
                                          grid.arrange(tP(), ncol=1)
                                          grid.arrange(vP(), ncol=1)
                                          grid.arrange(dP(), ncol=1)
                                          dev.off()
                                        })
  #############################
  
  ################
  # CREATE PLOTS #
  ################
  
  # create color vector for consistent color schemes across all plots
  tracecolors <- reactive({
    traces <- unique(tracedata()$trace)
    colorvector <- c("#377eb8", "#ff7f00")
    colorvector <- colorvector[1:length(traces)]
    names(colorvector) <- traces
    colorvector
  })
  
  # XY plot for traces - tab 1
  tP <- reactive({
    tP1 <- ggplot(traceDF(), aes(y = value, x = iter, fill = trace))+
      facet_wrap(~variable, scales = "free", ncol = 2)+
      theme_light()+
      theme(axis.title = element_blank(),
            legend.title = element_blank(),
            legend.spacing.x = unit(0.2, "cm"),
            axis.text = element_text(size = 11*scalefactor()*0.8),
            legend.text = element_text(size = 12*scalefactor()*1.1),
            legend.key = element_rect(size = 12*scalefactor()*0.8),
            strip.text = element_text(size = 12*scalefactor()))+
      scale_color_manual(values = tracecolors())
    
    # This adds points to XY plos only if this option was chosen in checkbox
    if("points" %in% input$traceplotstyle){
      tP1<- tP1 + geom_point(data = traceDF(), aes(y = value, x = iter, color = trace), size = input$cex)
    }
    
    # This adds lines to XY plos only if this option was chosen in checkbox
    if("lines" %in% input$traceplotstyle){
      tP1<- tP1 + geom_line(data = traceDF(), aes(y = value, x = iter, color = trace), size = input$cex/2)
    }
    
    tP1
  })
  
  # Violin plot for traces - tab 2
  vP <- reactive({
    vP1 <- ggplot(traceDF(), aes(y = value, x = trace, fill = trace))+
      facet_wrap(~variable, scales = "free", ncol = 2)+
      theme_light()+
      theme(axis.title = element_blank(),
            axis.text.x = element_blank(),
            legend.title = element_blank(),
            legend.spacing.x = unit(0.2, "cm"),
            axis.text = element_text(size = 11*scalefactor()*0.8),
            legend.text = element_text(size = 12*scalefactor()*1.1),
            legend.key = element_rect(size = 12*scalefactor()*0.8),
            strip.text = element_text(size = 12*scalefactor()))+
      scale_fill_manual(values = tracecolors())+
      guides(scale_color_manual())
    
    # Users can add Boxplots and/or datapoints to violin plot
    # Datapoints must always be first layer, so each combination of Violin/Boxplot/points is iterated below
    
    # No points, no boxplots
    if(!"boxplot" %in% input$violinplotstyle & !"points" %in% input$violinplotstyle){
      vP1 <- vP1 + 
        geom_violin(trim = TRUE, alpha = 0.5, color = NA) 
    }
    
    # With points, no boxplots
    if(!"boxplot" %in% input$violinplotstyle & "points" %in% input$violinplotstyle){
      vP1 <- vP1 + 
        geom_jitter(height = 0, width = 0.1, alpha = 0.2, color="gray", show.legend=FALSE)+
        geom_violin(trim = TRUE, alpha = 0.5, color = NA) 
    }
    
    # No points, with boxplots
    if("boxplot" %in% input$violinplotstyle & !"points" %in% input$violinplotstyle){
      vP1 <- vP1 + 
        geom_violin(trim = TRUE, alpha = 0.5, color = NA) + 
        geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = input$cex/2)
    }
    
    # With points and boxplots
    if("points" %in% input$violinplotstyle & "boxplot" %in% input$violinplotstyle){
      vP1 <- vP1 + 
        geom_jitter(height = 0, width = 0.1, alpha = 0.2, color = "gray", show.legend = FALSE)+
        geom_violin(trim = TRUE, alpha = 0.5, color = NA)+
        geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = input$cex/2)
    }
    
    vP1
  })

  # density plot for traces – tab 3
  dP <- reactive({
    ggplot(traceDF())+
      geom_density(aes(x = value, fill = trace), alpha = 0.5, size = 0)+
      facet_wrap(~variable, scales = "free", ncol = 2)+
      theme_light()+
      theme(axis.title = element_blank(),
            axis.text.y = element_blank(),
            legend.title = element_blank(),
            legend.spacing.x = unit(0.2, "cm"),
            axis.text = element_text(size=11*scalefactor()*0.8),
            legend.text = element_text(size=12*scalefactor()*1.1),
            legend.key = element_rect(size=12*scalefactor()*0.8),
            strip.text = element_text(size=12*scalefactor()))+
      scale_fill_manual(values = tracecolors())
  })
  ###############
  
  #################
  # DISPLAY PLOTS #
  #################
  
  # get height & width of plot area from slider inputs
  plotheight <- reactive({
    input$height
    })
  
  plotwidth <- reactive({
    input$width
    })
  
  # calculate scalefactor from height & width given – is used for e.g., line width, cex, etc in plots
  scalefactor <- reactive({(
    input$height+input$width)/2000
    })
  
  # render all 3 plots
  output$tracePlot <- renderPlot({
    tP()
    })
  
  output$tracePlot.ui <- renderUI({
    withSpinner(plotOutput("tracePlot", 
                           height = plotheight(), 
                           width = plotwidth()), 
                color = "#2C4152", size = 0.5)
  })
  
  output$violinPlot <- renderPlot({
    vP()
  })
  
  output$violinPlot.ui <- renderUI({
    withSpinner(plotOutput("violinPlot", 
                           height = plotheight(), 
                           width = plotwidth()), 
                color = "#2C4152", size = 0.5)
  })
  
  output$densePlot <- renderPlot({
    dP()
  })
  
  output$densePlot.ui <- renderUI({
    withSpinner(plotOutput("densePlot", 
                           height = plotheight(), 
                           width = plotwidth()), 
                color = "#2C4152", size = 0.5)
  })
  #################
  
  ######################
  # SUMMARY STATISTICS #  Tab 4
  ######################
  
  # Use DT package for flexible table formatting
  output$table <- DT::renderDataTable(selection = 'none', # do not select any column/row
                                      extensions = 'Buttons', # add download buttons
                                      options = list(dom = "Bt", # show Buttons on top, then table
                                                     ordering = FALSE, # don't add buttons to sort table
                                                     autoWidth = TRUE, 
                                                     buttons = c('copy', 'csv', 'pdf', 'print'), # options for downloading buttons
                                                     rowCallback = DT::JS( # the following defines the looks of the table: row names bold, and values outside of recommended range in red (see stats.md)
                                                       'function(row, data) {
                                                                           $("td:eq(0)", row).css("fontWeight", "bold");
                                                                           if (parseFloat(data[3]) < 100)
                                                                           $("td:eq(3)", row).css("color", "red");
                                                                           if (parseFloat(data[4]) > 2)
                                                                           $("td:eq(4)", row).css("color", "red");
                                                                           if (parseFloat(data[5]) > 2)
                                                                           $("td:eq(5)", row).css("color", "red");
                                                                           if (parseFloat(data[6]) > 0.3)
                                                                           $("td:eq(6)", row).css("color", "red");
                                                                           if (parseFloat(data[7]) > 1.2)
                                                                           $("td:eq(7)", row).css("color", "red");
                                                                           if (parseFloat(data[8]) > 1.2)
                                                                           $("td:eq(8)", row).css("color", "red");
                                                                           }'
                                                     )),
                                         {
                         req(input$tracefile)
                         
                         # spread the trace file                   
                         traceDF <- spread(tracedata(), variable, value)
                         
                         # calculate means & sds for all numeric columns
                         Mean <- summarise_if(traceDF,is.numeric, mean)
                         SD <- summarise_if(traceDF,is.numeric, sd)
                         
                         # same for ess
                         ess1 <- traceDF %>% 
                           filter(trace==chainnames()[[1]][1]) %>% 
                           summarize_if(is.numeric, effectiveSize)
                         
                         # ess for chain 2 (if present, merge 2 estimates)
                         try(ess2 <- traceDF %>% 
                               filter(trace==chainnames()[[2]][1]) %>% 
                               summarize_if(is.numeric, effectiveSize), 
                             silent = TRUE)
                         ess <- ess1
                         try(ess <- ess1 + ess2, silent = TRUE)
                         
                         # calculate geweke for chain 1 only
                         traceDF1 <- traceDF %>% filter(trace==chainnames()[[1]][1])
                         traceDF1 <- subset(traceDF1, select = -trace)
                         geweke1 <- geweke.diag(traceDF1)[[1]]
                         geweke1 <- abs(geweke1)
                         
                         # merge means, sd, ess, and geweke into results dataframe
                         results <- data.frame(t(rbind(Mean, SD, ess, geweke1)))
                         names(results) <- c("Mean", "SD", "ESS", paste("Geweke", chainnames()[[1]][1], sep=" "))
                         
                         # remove 'iter' & 'time' variables
                         results <- results[2:nrow(results),]
                         
                         if(length(input$tracefile$datapath)==2){
                           # calculate geweke for chain 2 only
                           traceDF2 <- traceDF %>% filter(trace==chainnames()[[2]][1])
                           traceDF2 <- subset(traceDF2, select = -trace)
                           geweke2 <- geweke.diag(traceDF2)[[1]]
                           geweke2 <- abs(geweke2)
                           geweke2 <- data.frame(geweke2[2:length(geweke2)])
                           names(geweke2) <- paste("Geweke", chainnames()[[2]][1], sep=" ")
                          
                           # add geweke results to dataframe
                           results <- cbind(results, geweke2)
                           
                           # calculate discrepancy according to phylobayes manual
                           # first, get means and sd for each variable and each chain sepately
                           means1 <- traceDF %>% filter(trace==chainnames()[[1]][1]) %>% summarize_if(is.numeric, mean)
                           means2 <- traceDF %>% filter(trace==chainnames()[[2]][1]) %>% summarize_if(is.numeric, mean)
                           sd1 <- traceDF %>% filter(trace==chainnames()[[1]][1]) %>% summarize_if(is.numeric, sd, na.rm = TRUE)
                           sd2 <- traceDF %>% filter(trace==chainnames()[[2]][1]) %>% summarize_if(is.numeric, sd, na.rm = TRUE)
                           
                           # do the calculation and add to results dataframe
                           Discrepancy <- 2*abs(means1 - means2) / (sd1 + sd2)
                           Discrepancy <- data.frame(t(Discrepancy[2:length(Discrepancy)]))
                           names(Discrepancy) <- "Discrepancy"
                           results <- cbind(results, Discrepancy)
                           
                           # finally, do Gelman & Rubin convergence estimate
                           # for this, chains must be of identical length, so first determine shorter chain
                           min_chainlength <- min(c(nrow(traceDF1), nrow(traceDF2)))
                           
                           # then limit the data to the minimal generation and combine as mcmc list
                           tracelist <- list(mcmc(traceDF1[1:min_chainlength,2:ncol(traceDF1)]), mcmc(traceDF2[1:min_chainlength,2:ncol(traceDF2)]))
                           
                           # Calculate Gelman & Rubin, extract point estimates and ci, rename, and combine with results dataframe 
                           gel.res<- gelman.diag(tracelist, autoburnin = F)
                           gel.point <- data.frame(gel.res$psrf[,1])
                           gel.ci <- data.frame(gel.res$psrf[,2])
                           names(gel.point) <- "Gelman & Rubin point estimate"
                           names(gel.ci) <- "Gelman & Rubin 95% CI"
                           gel.point <- tibble::rownames_to_column(gel.point)
                           gel.ci <- tibble::rownames_to_column(gel.ci)
                           results <- tibble::rownames_to_column(results)
                           results <- full_join(results, gel.point, by = "rowname")
                           results <- full_join(results, gel.ci, by = "rowname")
                           row.names(results) <- results$rowname
                           results <- results[,2:ncol(results)]
                         }
                         
                         # for the numeric values in dataframe, round using 2 decimals
                         is.num <- sapply(results, is.numeric)
                         results[is.num] <- lapply(results[is.num], round, 2)
                         
                         # call dataframe
                         results
                         })
  ######################
  

                                                  ###########################
                                                  #>>>>>>>>>>>>><<<<<<<<<<<<#
  #################################################>>>>>>> TREE TAB <<<<<<<<#################################################
                                                  #>>>>>>>>>>>>><<<<<<<<<<<<#
                                                  ###########################  
  
  
  ########################
  # TREE FILE PROPERTIES #
  ########################
  
  # read in tree file 1
  trees1 <- reactive({
    req(input$treefile)
    treefile1 <- input$treefile$datapath[1]
    trees1 <- read.tree(treefile1)
    # thin out trees by number given in tree thinning slider
    trees1.1 <- trees1[seq(from = 1, to = length(trees1), by= input$treethin)]
    trees1.1
  })
  
  # read also tree file2 if present  
  trees2 <- reactive({
    req(input$treefile)
    if(length(input$treefile$datapath)>1){
      treefile2 <- input$treefile$datapath[2]
      trees2 <- read.tree(treefile2)
      trees2.1 <- trees2[seq(from = 1, to = length(trees2), by= input$treethin)]
      trees2.1
    }
  })

  # determine number of tree generations (smaller number from both files) 
  ngentree <- reactive({
    treegen <- length(trees1())
    if(length(input$treefile$datapath)>1){
      treegen <- min(c(length(trees1()), length(trees2())))
    }
    treegen
    })
  
  # extract the tree tip labels (= taxon names) 
  tipnames <- reactive({
    req(input$treefile)
    tree1 <- trees1()[[1]]
    tiplabs <- list()
    labels <-sort(tree1$tip.label)
    for(i in 1:length(labels)){
      tiplabs[[i]] <- labels[i]
    }
    names(tiplabs) <- tiplabs
    tiplabs
  })
  
  # get the names of the chains from the tree file names
  treenames <- reactive({
    req(input$treefile)
    strsplit(input$treefile$name[1:length(input$treefile$name)], ".treelist", fixed=T)
  })
  
  ########################  
  
  ############################
  # UI ELEMENTS FOR TREE TAB #
  ############################
  
  # Slider that determines the tree generation currently displyed
  output$treegens <- renderUI({
    req(input$treefile)
    sliderInput("generation", "Tree generation", 
                min = 1, 
                max = ngentree(), 
                value = 1, 
                step = 1, 
                animate = animationOptions(interval = 1000,  # this adds a button that animates an iteration through all tree generations
                                           loop = TRUE, 
                                           playButton = tags$h4(tags$b("> PLAY <")), pauseButton = tags$h4(tags$b("> PAUSE <"))))
                                           # custom play buttons with larger font  
  })
  
  # Slider that determines the number of burnin generations
  output$conburnin <- renderUI({
    req(input$treefile)
    sliderInput("conburnin", "Burnin [# of iterations]:", 
                min = 0, 
                max = ngentree(),  
                step = 100, 
                value = ngentree()/5) # default burnin = 20% of all trees
  })
  
  # Picker input to chose outgroup
  output$outgroups <- renderUI({
    req(input$treefile)
    pickerInput("outgroup", 
                label = "Select outgroup(s) for rooting",
                choices=c(tipnames()), 
                multiple=TRUE, 
                options = list(`selected-text-format` = "count > 2",        # show names of up to 2 outgroup taxa
                               `count-selected-text` = "{0} outgroup taxa", # show number of outgroups if more than 2 taxa are selected
                               `live-search` = TRUE),                       # this enables a search box in the selector
                inline = FALSE)
  })
  
  # Similar picker input to chose which taxa to highlight
  output$highlight <- renderUI({
    req(input$treefile)
    pickerInput("highlight", 
                label = HTML("Select taxa to highlight <font color=\"#377eb8\">[blue]</font>"),
                choices=tipnames(), 
                multiple=TRUE, 
                options = list(`selected-text-format` = "count > 2",
                               `count-selected-text` = "{0} taxa",
                               `actions-box` = TRUE,     # enables 'select-all' and 'select none' buttons
                               `live-search` = TRUE),  
                inline = FALSE)
  })
  
  # And another picker for another highlight
  output$highlight2 <- renderUI({
    req(input$treefile)
    pickerInput("highlight2", 
                label = HTML("Select taxa to highlight <font color=\"#ff7f00\">[orange]</font>"),
                choices=tipnames(), 
                multiple=TRUE, 
                options = list(`selected-text-format` = "count > 2",
                               `count-selected-text` = "{0} taxa",
                               `actions-box` = TRUE,
                               `live-search` = TRUE), 
                inline = FALSE)
  })
  
  # PDF download handle for tree plots
  output$downloadtreePDF <- downloadHandler(filename = "treeplots.pdf",
                                            contentType = "application/pdf",
                                            content = function(file1) {
                                              file.copy("treeplot.pdf", file1)
                                        })
  ############################  
  
  ################################
  # FURTHER TREE DISPLAY OPTIONS #
  ################################
  
  # Reset outgroup only when button is presed, default is no outgroup
  og <- reactiveValues(outgroup = "<None>")
  
  # observe midpoint button
  observeEvent(input$midpoint, {
    og$outgroup <- "<Midpoint>"
  })
  
  # observe unroot button
  observeEvent(input$unroot, {
    og$outgroup <- "<None>"
  })
  
  # observe reroot button
  observeEvent(input$outgroup, {
    og$outgroup <- input$outgroup
  })
  
  # this gets height and width for plots from the slider inputs and creates a scale factor to be used in plots
  treeheight <- reactive({
    input$treeheight
  })
  
  treewidth <- reactive({
    input$treewidth
  })
  
  treescalefactor <- reactive({(
    input$treeheight+input$treewidth)/1800
  })
  
  ################################  
  
  ##############
  # TREE PLOTS #
  ##############
  
  ###### TAB 1 
  # plot 1 is for single tree plots per generation
  output$treeplot <- renderPlot({
    
    # require these for interactive rooting
    req(input$treefile)
    input$reroot
    input$midpoint
    input$unroot
    
    # get outgroup from ui, but isolate so that rerooting is only done when button is pressed
    isolate(outgroup <- og$outgroup)
    
    # change plot layout to 2 columns if 2 treefiles are present
    if(length(input$treefile$name)==2){
      par(mfrow=c(1,2))
    }
    
    # unroot tree if no root was chosen
    if(outgroup=="<None>"){
      currtree1 <- unroot(trees1()[[input$generation]])
    }
    
    # midpoint root if chosen
    if(outgroup=="<Midpoint>"){
      currtree1 <- midpoint.root(trees1()[[input$generation]])
    }
    
    # root with outgroup if chosen
    if(outgroup!="<Midpoint>" & outgroup!="<None>"){
      currtree1<- root(trees1()[[input$generation]], 
           outgroup = outgroup, 
           resolve.root = TRUE)
    }
    
    # Create tip label color vector from highlight slider options (need to do this after rooting, as this impacts tip label order)
    concolors1.1 <- currtree1$tip.label %in% input$highlight
    concolors1.2 <- currtree1$tip.label %in% input$highlight2
    concolvec1 <- rep("black", length=length(currtree1$tip.label))
    concolvec1[which(concolors1.1==TRUE)] <- "#377eb8"
    concolvec1[which(concolors1.2==TRUE)] <- "#ff7f00"
    concolvec1[which(concolors1.1==TRUE & concolors1.2==TRUE)] <- "#4daf4a"
    
    # Plot
    plot(ladderize(currtree1), 
         main = paste(treenames()[[1]], "iteration", input$generation), 
         cex = input$treecex, 
         cex.main = input$treecex*1.1,
         edge.width = input$treecex, 
         tip.color = concolvec1)
    add.scale.bar(lwd = input$treecex)
    
    # now repeat all for plot from chain2!
    if(length(input$treefile$name)==2){
      
      # no root
      if(outgroup=="<None>"){
        currtree2 <- unroot(trees2()[[input$generation]])
      }
      
      # midpoint root
      if(outgroup=="<Midpoint>"){
        currtree2 <- midpoint.root(trees2()[[input$generation]])
      }
      
      # outgroup root
      if(outgroup!="<Midpoint>" & outgroup!="<None>"){
        currtree2 <- root(trees2()[[input$generation]], 
                         outgroup = outgroup, 
                         resolve.root = TRUE)
      }
      
      # color vector
      concolors2.1 <- currtree2$tip.label %in% input$highlight
      concolors2.2 <- currtree2$tip.label %in% input$highlight2
      concolvec2 <- rep("black", length=length(currtree2$tip.label))
      concolvec2[which(concolors2.1==TRUE)] <- "#377eb8"
      concolvec2[which(concolors2.2==TRUE)] <- "#ff7f00"
      concolvec2[which(concolors1.1==TRUE & concolors1.2==TRUE)] <- "#4daf4a"
     
      # plot
      plot(ladderize(currtree2), 
           cex = input$treecex, 
           main = paste(treenames()[[2]], "iteration", input$generation),
           cex.main = input$treecex*1.1,
           edge.width = input$treecex, 
           tip.color = concolvec2)
      add.scale.bar(lwd = input$treecex)
    }  
    
    # Copy plot to device
    dev.copy2pdf(file = "treeplot.pdf", 
                 height=treeheight()/72, 
                 width=treewidth()/72)
    
  })
  
  
  # render the plot with spinner & using the height and widths from ui
  output$treePlot.ui <- renderUI({
    req(input$treefile)
    withSpinner(plotOutput("treeplot", 
                           height = treeheight(), 
                           width = treewidth()), 
                color = "#2C4152", 
                size = 0.5)
  })
  
  ###### TAB 2 
  # Plot 2 is a consensus plot 
  
  # first calculate the tree
  contree <- reactive({
    # require these for interactive rooting
    req(input$treefile)
    input$reroot
    input$midpoint
    input$unroot
    isolate(outgroup <- og$outgroup)
    
    # here, remove burnin from trees before plotting 
    trees1 <- trees1()[(input$conburnin+1):length(trees1())]
    
    # if 2 tree files are present, merge
    if(length(input$treefile$name)==2){
      trees2 <- trees2()[(input$conburnin+1):length(trees2())]
      trees1 <- c(trees1, trees2)
    }
    
    # get consensus branch lengths 
    contree1 <- consensus.edges(trees1, method="least.squares")
    
    # and count how often the nodes are present in all trees (=pp) and writes this as nodelabels to the tree
    sv <- prop.clades(contree1, trees1)
    contree1$node.label <- sv/length(trees1)
    contree1$node.label <- formatC(contree1$node.label, digits = 2, format = "f") # 2 decimals for pp values
    
    # do the rooting
    # unroot
    if(outgroup=="<None>"){
      contree1 <- unroot(contree1)
    }
    
    # midpoint root
    if(outgroup=="<Midpoint>"){
      contree1 <- midpoint.root(contree1)
    }
    
    # root with outgroup
    if(outgroup!="<Midpoint>" & outgroup!="<None>"){
      contree1<- root(contree1, 
                       outgroup = outgroup, 
                       resolve.root = TRUE)
    }
    contree1
  })
  
  # now render consensus plot
  output$consensusPlot <- renderPlot({
    # colorvector again
    concolors1.1 <- contree()$tip.label %in% input$highlight 
    concolors1.2 <- contree()$tip.label %in% input$highlight2
    concolvec1 <- rep("black", length=length(contree()$tip.label))
    concolvec1[which(concolors1.1==TRUE)] <- "#377eb8"
    concolvec1[which(concolors1.2==TRUE)] <- "#ff7f00"
    concolvec1[which(concolors1.1==TRUE & concolors1.2==TRUE)] <- "#4daf4a"
    
    # plot
    plot(ladderize(contree()), 
         main = paste("Consensus tree, burnin =",input$conburnin),  # this adds the burnin to tree title
         cex.main = input$treecex*1.1,
         cex = input$treecex, 
         edge.width = input$treecex, 
         tip.color = concolvec1)
    nodelabels(text = contree()$node.label,  # some formatting for the pp values
               frame = "none", 
               adj = 0, 
               cex = input$treecex*0.8)
    add.scale.bar(lwd = input$treecex)  
    
    # Copy plot to device
    dev.copy2pdf(file = "treeplot.pdf", 
                 height=treeheight()/72, 
                 width=treewidth()/72)
    
  })
  
  # render plot 2 with spinner
  output$consensusPlot.ui <- renderUI({
    withSpinner(plotOutput("consensusPlot", 
                           height=treeheight(), 
                           width=treewidth()*0.75), 
                color="#2C4152", 
                size=0.5)
  }) 
  
  # also export the tree as newick if wanted
  output$newick <- downloadHandler(
    filename = "consensus.tre",
    content = function(file4) {
      write.tree(contree(), file = file4, append = FALSE,
                 digits = 10, tree.names = FALSE)
    }
  )


  
  ###### TAB 3 
  # Plot 3 shows 2 consensus plots with comparison 
  
  output$differencePlot <- renderPlot({
    
    # required for interactive rooting
    input$reroot
    input$midpoint
    input$unroot
    isolate(outgroup <- og$outgroup)
    
    # check if 2 tree files are present and if rooting other than midpoint is selected. Throw error if not.
    validate(
      need(length(input$treefile$name)==2, "Upload tree samples from 2 Phylobayes runs to display differences between consensus topologies.")
    )
    validate(
      need(outgroup!="<Midpoint>", "Midpoint rooting not possible for consensus cladograms.\nPlease choose different root.")
    )

    # create majority rule consensus for tree from both chains
    trees1 <- trees1()[(input$conburnin+1):length(trees1())]
    contree1 <- consensus(trees1, p = 0.5)
    trees2 <- trees2()[(input$conburnin+1):length(trees2())]
    contree2 <- consensus(trees2, p = 0.5)
    
    # reroot outgroup
    if(outgroup!="<Midpoint>" & outgroup!="<None>"){
      contree1 <- root(contree1, outgroup=outgroup)
      contree2 <- root(contree2, outgroup=outgroup)
    }
    
    # unroot tree  
    if(outgroup=="<None>"){
      contree1 <- unroot(contree1)
      contree2 <- unroot(contree2)
    }
    
    # colorvector
    colors1.1 <- contree1$tip.label %in% input$highlight 
    colors1.2 <- contree1$tip.label %in% input$highlight2
    colvec1 <- rep("black", length=length(contree1$tip.label))
    colvec1[which(colors1.1==TRUE)] <- "#377eb8"
    colvec1[which(colors1.2==TRUE)] <- "#ff7f00"
    colvec1[which(colors1.1==TRUE & colors1.2==TRUE)] <- "#4daf4a"
    colors2.1 <- contree2$tip.label %in% input$highlight 
    colors2.2 <- contree2$tip.label %in% input$highlight2
    colvec2 <- rep("black", length=length(contree1$tip.label))
    colvec2[which(colors2.1==TRUE)] <- "#377eb8"
    colvec2[which(colors2.2==TRUE)] <- "#ff7f00"
    colvec2[which(colors2.1==TRUE & colors2.2==TRUE)] <- "#4daf4a"
    
    
    # function that plots 2 trees next to each other and hifhlights differences 
    # modified from http://blog.phytools.org/
    phylo.diff.new <-function (x, y, main1, main2, coltip1, coltip2, ...){
      uniqT1 <- distinct.edges(x, y)
      uniqT2 <- distinct.edges(y, x)
      treeA.cs <- rep("black", dim(x$edge)[1])
      treeA.cs[uniqT1] = "#FF1493"
      treeA.lw <- rep(input$treecex, dim(x$edge)[1])
      treeA.lw[uniqT1] = input$treecex*2
      treeB.cs <- rep("black", dim(y$edge)[1])
      treeB.cs[uniqT2] = "#FF1493"
      treeB.lw <- rep(input$treecex, dim(x$edge)[1])
      treeB.lw[uniqT1] = input$treecex*2
      par(mfrow = c(1, 2))
      plot(x, edge.color = treeA.cs, main=main1, tip.color=coltip1, edge.width= treeA.lw, ...)
      plot(y, edge.color = treeB.cs, main=main2, tip.color=coltip2, edge.width= treeB.lw, direction = "leftwards", ...)
      invisible()
    }
    
    # plot 2 consensus tree next to each other, highlighting the nodes in each that are uniq
    phylo.diff.new(contree1, contree2, 
                   cex = input$treecex, 
                   cex.main = input$treecex*1.1,
                   main1 = treenames()[[1]], 
                   main2 = treenames()[[2]], 
                   coltip1 = colvec1, 
                   coltip2 = colvec2)
    
    # Copy plot to device
    dev.copy2pdf(file = "treeplot.pdf", 
                 height=treeheight()/72, 
                 width=treewidth()/72)
    
  })
  
  # finally, render this plot with spinner
  output$differencePlot.ui <- renderUI({
    req(input$treefile[[2]])
    withSpinner(plotOutput("differencePlot", height=treeheight(), width=treewidth()), color="#2C4152", size=0.5)
  })   
  

  
  ###### TAB 4 
  # Plot 4 shows differences between trees across iterations and between chains 
  
  rP <- reactive({
    
    req(input$treefile)

    # read in tree file and create vector for results
    trees1 <- trees1()[(input$conburnin+1):length(trees1())]
    rf <- vector(mode = "numeric")
    
    # now loop through each tree and calculate R-F distance between tree N and tree (N-1), store results in vector
    for(i in 2:length(trees1)){
      rf[i-1] <- dist.topo(trees1[i-1], trees1[i])  
    }
    
    # label with generation number (including burnin), merge into dataframe, and create new variable with chain name
    x <- (2+input$conburnin):(length(trees1)+input$conburnin)
    RFdf<- as.data.frame(cbind(x, rf))
    RFdf$chain <- paste("Iteration n vs. Iteration (n-1),", treenames()[[1]])
    
    # create xy plot of RF distance over generation
    RFggplot <- ggplot(data=RFdf, aes(x=x, y=rf, color=chain)) +
      geom_line(size = input$treecex, alpha=0.7)+
      theme_light() +
      xlab("Tree generation") +
      ylab("Robinson-Foulds distance")+
      theme(axis.title = element_text(size = 12*treescalefactor()),
            legend.title = element_blank(),
            legend.spacing.x = unit(0.2, "cm"),
            axis.text = element_text(size = 11*treescalefactor()),
            legend.position = "bottom",
            legend.direction = "vertical",
            legend.text = element_text(size = 12*treescalefactor()))+
      scale_color_manual(values = c("#377eb8", "#ff7f00", "#4daf4a"))
    
    # do all of these steps again for second tree file if present
    if(length(input$treefile$name)==2){
      
      # read in tree file and create vector for results
      trees2 <- trees2()[(input$conburnin+1):length(trees2())]
      rf2 <- vector(mode = "numeric")
      
      # loop through each tree and calculate R-F distance between tree N and tree (N-1), store results in vector
      for(i in 2:length(trees2)){
        rf2[i-1] <- dist.topo(trees2[i-1], trees2[i])  
      }
      
      # label with generation number (including burnin), merge into dataframe, and create new variable with chain name
      x <- (2+input$conburnin):(length(trees2)+input$conburnin)
      RFdf2<- as.data.frame(cbind(x, rf2))
      RFdf2$chain <-  paste("Iteration n vs. Iteration (n-1),", treenames()[[2]])
      
      # determine min generation time from both chains
      maxgen <- min(c(length(trees1), length(trees2)))
      
      # results vector
      rf3 <- vector(mode = "numeric")
      
      # for each generation there was a tree generated for both chains, calculate the RF distance between the trees of the 2 chains
      for(i in 1:maxgen){
        rf3[i] <- dist.topo(trees1[i], trees2[i])
      }
      
      # label with generation number (including burnin), merge into dataframe, and create new variable with chain name
      x <- (1+input$conburnin):(length(rf3)+input$conburnin)
      RFdf3<- as.data.frame(cbind(x, rf3))
      RFdf3$chain <- paste(treenames()[[1]], treenames()[[2]], sep = " vs. ")
      
      # add 2 additional lines to already existing plot
      RFggplot <- RFggplot + 
        geom_line(data=RFdf2, aes(x=x, y=rf2), size=input$treecex, alpha=0.7)+ 
        geom_line(data=RFdf3, aes(x=x, y=rf3), size=input$treecex, alpha=0.7)
    }
    
    # call plot
    RFggplot 
  
  })
  
  # render plot
  output$rfPlot <- renderPlot({
    rP()
  })
  
  output$rfPlot.ui <- renderUI({
    req(input$treefile[[1]])
    withSpinner(plotOutput("rfPlot", height=treeheight()/2, width=treewidth()), color="#2C4152", size=0.5)
  }) 
  
  # create download button for plot
  output$downloadrfplot <- downloadHandler(filename = "RFplot.pdf",
                                        content = function(file3) {
                                          pdf(file3,  height=input$height/144, width=input$width/72)
                                          grid.arrange(rP(), ncol=1)
                                          dev.off()
                                        })
  
} 

# Run the application 
shinyApp(ui = ui, server = server)


#   LICENSE INFORMATION
#
#   Copyright 2018 Michael Gerth
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy of 
#   this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights to 
#   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
#   of the Software, and to permit persons to whom the Software is furnished to do 
#   so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all 
#   copies or substantial portions of the Software.
# 
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
#   SOFTWARE.
