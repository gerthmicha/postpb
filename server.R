server <- function(input, output, session) {

  # GENERAL OPTIONS ######
  source("app.R", local = TRUE)

  # increase maximum upload size
  options(shiny.maxRequestSize = 500 * 1024^2)

  # shut down R session when browser window is closed
  session$onSessionEnded(function() {
    stopApp()
  })

  # Prefer regular to scientific number connotation, except for very large numbers
  options(scipen = 5)

  # >>> TRACE TAB <<<######

  #  TRACE FILE PROPERTIES ------------------------------

  # get input files
  # create reactive value for input files
  tracefile <- reactiveValues()

  # if input files are uploaded, use file paths and names from these
  observeEvent(input$tracefile, {
    tracefile$datapath <- input$tracefile$datapath
    tracefile$name <- input$tracefile$name
  }) 

  # if 'example 1' button is pressed, load example 1 from example folder
  observeEvent(input$exampletrace1, {
    tracefile$datapath <- list.files("example/", "\\.trace\\>", full.names = TRUE)
    tracefile$name <- list.files("example/", "\\.trace\\>", full.names = FALSE)
  })

  # if 'example 2' button is pressed, load example 2 from example folder
  observeEvent(input$exampletrace2, {
    tracefile$datapath <- list.files("example/", "\\.p\\>", full.names = TRUE)
    tracefile$name <- list.files("example/", "\\.p\\>", full.names = FALSE)
  })

  # if 'info' button is pressed, show popup that explains examples
  observeEvent(input$exampletraceinf, {
    showModal(modalDialog(
      title = "postpb comes with two sets of example data:",
      includeMarkdown("example/examples.md"),
      size = "m"
    ))
  })

  # get number of generations from tracefiles
  ngen <- reactive({
    req(tracefile$datapath)
    get.ngen(tracefile)
  })

  # Get names of chains from file names provided
  chainnames <- reactive({
    # get names of chains from trace file names
    strsplit(tracefile$name[1:length(tracefile$name)], "\\.trace\\>|\\.p\\>")
  })

  # Reading the trace files
  tracedata_raw <- reactive({
    req(tracefile$datapath)
    read.trace(tracefile, chainnames())
  })

  # Thinning the tracefiles
  tracedata_thin <- reactive({
    # Applies thinning
    tracelist <- lapply(tracedata_raw(), thin.trace, tracethin())

    # rename first column in all dfs to "iter"
    lapply(tracelist, function(dflist) {
      names(dflist)[1] <- "iter"
      dflist
    })
  })

  # Removing burnin from the tracefiles, and merge all into df
  tracedata <- reactive({
    tracelist.as.df(lapply(tracedata_thin(), burn.trace, input$burnin))
  })

  # filter tracedata to only plot selected trace files in checkbox (default= select all)
  traceDF <- reactive({
    req(tracefile$datapath)
    choose.trace(tracedata(), input$whichchain)
  })


  #| UI ELEMENTS FOR TRACE TAB SIDEBAR ------------------------------

  # display burnin slider using the number of generations read from trace file
  output$burnin <- renderUI({
    req(tracefile)
    sliderInput("burnin", "Burnin [# of iterations]:",
      min = 0,
      max = trunc(ngen() / 10),
      value = trunc(ngen() / 10 * 0.2), # default = 20% of iterations
      # default for step 10 for 100–999 gens, 100 for 1000-9999 etc
      step = 10^(ceiling(log10(ngen() / 10 * 0.02)))
    )
  })

  # Set the default value of trace file thinning to 10
  tracethin <- reactiveVal(10)

  # Update the tree thinning value for calculating the consensus
  # only after burnin has changed
  observeEvent(input$burnin, {
    tracethin(input$tracethin)
  })

  # Update the burnin after tree thinning button is pressed
  observeEvent(input$trecalc, {
    thinval <- input$tracethin
    # Update burnin slider
    updateSliderInput(session, "burnin",
      label = "Burnin [# of iterations]:",
      min = 0,
      max = trunc(ngen() / thinval),
      step = 10^(ceiling(log10(ngen() / thinval * 0.02))),
      value = ngen() / thinval * 0.2
    )
  })

  # display checkbox to select which trace file to plot
  output$whichchain <- renderUI({
    req(tracefile$datapath[2])
    names <- lapply(chainnames(), `[[`, 1)
    prettyCheckboxGroup(
      inputId = "whichchain",
      label = "Select trace file[s] to plot",
      choices = names,
      selected = names,
      inline = FALSE,
      status = "primary",
      icon = icon("check")
    )
  })

  # Button that toggles graphics parameters for the statistics
  observeEvent(input$traceplotsopts, ignoreInit = TRUE, {
    toggle("tplotopts")
  })

  # Button that toggles explanations for the statistics
  observeEvent(input$explanation, {
    toggle("stats")
  })


  #| CREATE PLOTS ------------------------------

  #| General display options ------------------------------

  # get height & width of plot area from slider inputs
  plotheight <- reactive({
    input$height
  })
  plotwidth <- reactive({
    input$width
  })
  plotcex <- reactive({
    input$cex
  })

  # calculate scale factor from height & width given – is used for e.g., line width, cex, etc in plots
  scalefactor <- reactive({
    (input$height + input$width) / 2000
  })

  # Set theme for all trace plots
  tracetheme <- reactive({
    req(scalefactor())
    theme_light() +
      theme(
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.2, "cm"),
        axis.text = element_text(size = 11 * scalefactor() * 0.8),
        legend.text = element_text(size = 12 * scalefactor() * 1.1),
        legend.key = element_rect(size = 12 * scalefactor() * 0.8),
        strip.text = element_text(size = 12 * scalefactor())
      )
  })

  # create color vector for consistent color schemes across all plots
  tracecolors <- reactive({
    set.trace.colors(tracedata())
  })

  #| #  Tab 1 (Trace) -----
  tP <- reactive({
    xyplot(traceDF(), tracecolors(), tracetheme(), input$traceplotstyle, input$facetcol, input$cex)
  })
  output$tracePlot <- renderPlot({
    tP()
  })
  output$tracePlot.ui <- renderUI({
    shinycssloaders::withSpinner(plotOutput("tracePlot",
                                            height = input$height,
                                            width = input$width),
                                 color = "#2C4152", size = 0.5)
    })

  #| # Tab 2 (Violin) -----
  vP <- reactive({
    violinplot(traceDF(), tracecolors(), tracetheme(), input$facetcol, input$cex, input$violinplotstyle)
  })
  output$violinPlot <- renderPlot({
    vP()
  })
  output$violinPlot.ui <- renderUI({
    shinycssloaders::withSpinner(plotOutput("violinPlot",
                                            height = input$height,
                                            width = input$width),
                                 color = "#2C4152", size = 0.5)
  })

  #| # Tab 3 (Density) -----
  # density plot for traces
  dP <- reactive({
    densityplot(traceDF(), tracecolors(), tracetheme(), input$facetcol)
  })
  output$densePlot <- renderPlot({
    dP()
  })
  output$densePlot.ui <- renderUI({
    shinycssloaders::withSpinner(plotOutput("densePlot",
                                            height = input$height,
                                            width = input$width),
                                 color = "#2C4152", size = 0.5)
  })

  # Create pdf download handle for plots
  output$downloadPDF <- downloadHandler(
    filename = "traceplots.pdf",
    content = function(file) {
      pdf(file, height = input$height / 72, width = input$width / 72)
      grid.arrange(tP(), ncol = 1)
      grid.arrange(vP(), ncol = 1)
      grid.arrange(dP(), ncol = 1)
      dev.off()
    }
  )

  #| SUMMARY STATISTICS (Tab 4) ------------------------------

  output$table <- DT::renderDataTable({
    req(tracefile)

    # calculate summary stats
    sum.stats <- calc.sum.stats(traceDF())

    # Use DT package for flexible table formatting
    style.table(sum.stats, traceDF())
  })



  # >>> TREE TAB <<<#####


  #| TREE FILE PROPERTIES #------------------------------

  treefile <- reactiveValues()
  example <- reactiveValues()

  # if input files are uploaded, use file paths and names from these
  observeEvent(input$treefile, {
    treefile$datapath <- input$treefile$datapath
    treefile$name <- input$treefile$name
    example$click <- 0
  })

  # if 'example' button is pressed, load example from example folder
  observeEvent(input$exampletree1, {
    treefile$datapath <- list.files("example/", "\\.treelist\\>", full.names = TRUE)
    treefile$name <- list.files("example/", "\\.treelist\\>", full.names = FALSE)
    example$click <- 1
    showNotification("Loading example 1. This may take some time.",
      id = "ex2",
      duration = NULL
    )
  })

  observeEvent(input$exampletree2, {
    treefile$datapath <- list.files("example/", "\\.t\\>", full.names = TRUE)
    treefile$name <- list.files("example/", "\\.t\\>", full.names = FALSE)
    example$click <- 2
    showNotification("Loading example 2. This may take some time.",
      id = "ex2",
      duration = NULL
    )
  })

  # if 'info' button is pressed, show popup that explains examples
  observeEvent(input$exampletreeinf, {
    showModal(modalDialog(
      title = "postpb comes with two sets of example data:",
      includeMarkdown("example/examples.md"),
      size = "m"
    ))
  })


  # read in tree files
  completetrees <- reactive({
    req(treefile$datapath)

    # only update tree format when new files are uploaded
    input$treefile
    isolate(read.treefiles(treefile))
  })

  alltrees <- reactive({
    req(length(completetrees()) >= 1)
    if (length(completetrees()) == 1) {
      alltrees <- completetrees()
    }

    if (length(completetrees()) > 1 & !is.null(input$whichtree)) {
      alltrees <- completetrees()[which(treenames() %in% input$whichtree)]
    }

    if (is.null(input$whichtree) & length(completetrees()) > 1) {
      alltrees <- completetrees()
    }

    alltrees
  })

  # Tree thinning
  thintrees <- reactive({
    req(alltrees())
    thin.trees(alltrees(), treethin())
  })

  # After tree files have been read in, reset tree format radio buttons
  observeEvent(alltrees(), {
    updatePrettyRadioButtons(session, "treefiletype",
      choices = c("Newick (e.g., Phylobayes)", "Nexus (e.g., MrBayes)"),
      selected = "character(0)",
      inline = TRUE,
      prettyOptions = list(
        shape = "round",
        status = "primary"
      )
    )
    # Also, remove Notifications
    removeNotification("ex1")
    removeNotification("ex2")
  })

  # determine number of tree generations (smallest number from all files)
  ngentree <- reactive({
    alltrees()
    treegen <- min(unlist(lapply(alltrees(), length)))
    treegen
  })

  # extract the tree tip labels (= taxon names)
  tipnames <- reactive({
    req(treefile$datapath)
    tree1 <- completetrees()[[1]][[1]]
    labels <- sort(tree1$tip.label)

    # pickerInput requires list
    as.list(labels)
  })

  # get the names of the chains from the tree file names
  treenames <- reactive({
    req(treefile$datapath)
    strsplit(treefile$name[1:length(treefile$name)], "\\.treelist\\>|\\.t\\>")
  })


  #| UI ELEMENTS FOR TREE TAB SIDEBAR #------------------------------

  # Some tree calculations benefit from parallelization.
  # Get the number of cores interactively, and set the default to ( total cores -1 )

  totalcores <- reactive({
    parallel::detectCores()[1]
  })

  observe({
    totalcores()
    doParallel::registerDoParallel(cores = input$ncores)
  })

  output$ncores <- renderUI({
    numericInput("ncores",
      label = "Number of cores",
      min = 1,
      max = totalcores(),
      value = totalcores() - 1
    )
  })


  # display check box to select which tree file to plot
  output$whichtree <- renderUI({
    req(treefile$datapath[2])
    treenames <- lapply(treenames(), `[[`, 1)
    shinyWidgets::prettyCheckboxGroup(
      inputId = "whichtree",
      label = "Select tree file[s] to analyse",
      choices = treenames,
      selected = treenames,
      inline = FALSE,
      status = "primary",
      icon = icon("check")
    )
  })

  # Slider that determines the number of burnin generations
  output$conburnin <- renderUI({
    req(treefile$datapath)
    sliderInput("conburnin", "Burnin [# of iterations]:",
      min = 0,
      max = trunc(ngentree() / 10),
      step = 10^(ceiling(log10(ngentree() / 10 * 0.02))),
      value = ngentree() / 10 * 0.2
    ) # default burnin = 20% of all trees
  })

  # Set the default value of tree thinning to 10
  treethin <- reactiveVal(10)

  # Update the tree thinning value for calculating the consensus
  # only if burnin has changed
  observeEvent(input$conburnin, {
    treethin(input$treethin)
  })

  # Update the burnin after tree thinning button is pressed
  observeEvent(input$recalc, {
    thinval <- input$treethin
    treegens <- ngentree()
    # Burnin slider
    updateSliderInput(session, "conburnin",
      label = "Burnin [# of iterations]:",
      min = 0,
      max = trunc(treegens / thinval),
      step = 10^(ceiling(log10(ngentree() / 10 * 0.02))),
      value = treegens / thinval * 0.2
    )
  })


  # Picker input to chose outgroup
  output$outgroups <- renderUI({
    req(treefile$datapath)
    pickerInput("outgroup",
      label = "Select outgroup(s) for rooting",
      choices = c(tipnames()),
      multiple = TRUE,
      options = list(
        `selected-text-format` = "count > 2", # show names of up to 2 outgroup taxa
        `count-selected-text` = "{0} outgroup taxa", # show number of outgroups if more than 2 taxa are selected
        `live-search` = TRUE
      ), # this enables a search box in the selector
      inline = FALSE
    )
  })

  # After outgroup rooting has been performed, reset chosen outgroups
  observeEvent(roottree(), {
    updatePickerInput(session, "outgroup",
      label = "Select outgroup(s) for rooting",
      choices = c(tipnames())
    )
  })


  # Similar picker input to chose which taxa to highlight
  output$highlight <- renderUI({
    req(treefile$datapath)
    pickerInput("highlight",
      label = HTML("Select taxa to highlight"),
      choices = tipnames(),
      multiple = TRUE,
      options = list(
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} taxa",
        `actions-box` = TRUE, # enables 'select-all' and 'select none' buttons
        `live-search` = TRUE
      ),
      inline = FALSE
    )
  })

  # Button that toggles graphics parameters for the tree plots
  observeEvent(input$treeplotsopts, ignoreInit = TRUE, {
    toggle("trplotopts")
  })

  #| FURTHER TREE DISPLAY OPTIONS #------------------------------

  # Reset outgroup only when button is pressed, default is no outgroup
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

  # observe as new tree files are being read in
  # Also reset outgroups when new tree files are loaded
  observeEvent(completetrees(), {
    og$outgroup <- "<None>"
  })

  # get further tree options from UI ("Tree plots")
  treeopts <- reactive({
    set.tree.opts(input$treeopts)
  })

  # get further tree options from UI ("Tree labels")
  treefont <- reactive({
    set.tree.font(input$treefont)
  })

  # this gets height and width for plots from the slider inputs and creates a scale factor to be used in plots
  treeheight <- reactive({
    input$treeheight
  })

  treewidth <- reactive({
    input$treewidth
  })

  treescalefactor <- reactive({
    (input$treeheight + input$treewidth) / 1800
  })


  #| TREE PLOTS ------------------------------

  #| Interactive features -----
  # Enable interactive selection of taxon selection for rooting
  
  # create a dataframe of tip labels
  tipDF <- reactive({
    req(input$plot_brush)
    get.tip.df(roottree()) 
  })

  # get the names of selected tips from the interactive 'brush' click+drag
  # important here are only the ymin & ymax values: ape plots each tree from 0 (bottom) to Ntaxa, and difference between taxa is always 1
  conroot <- reactive({
    req(input$plot_brush)
    # get names
    tipDF <- tipDF() %>%
      filter(tiporder >= input$plot_brush$ymin) %>%
      filter(tiporder <= input$plot_brush$ymax)

    as.character(tipDF$tips)
  })

  # Once selected, use popup to decide what to do with selected taxa (highlight or rooting)
  observeEvent(input$plot_brush, {
    showModal(modalDialog(
      title = NULL,
      HTML(paste("You have selected ", length(conroot()), " taxa\n", "<br><br>")),
      align = "center",
      actionButton("rootbutton", HTML("Re-root with<br>selection"),
        width = "150px",
        style = "border:2px solid; border-radius: 4px; margin:5px; background-color:white; color:black; font-weight:bold"
      ),
      actionButton("colbutton",
        HTML("Higlight<br>selection"),
        width = "150px",
        style = "border:2px solid; border-radius: 4px; margin:5px; background-color:white; color:black; font-weight:bold"
      ),
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel")
      )
    ))
  })

  # If any of the buttons is pressed, close popup
  observeEvent(input$rootbutton, {
    removeModal()
  })
  observeEvent(input$colbutton, {
    removeModal()
  })

  # With this selection made, update selection in outgroup pickerInput...
  observeEvent(input$rootbutton, {
    updatePickerInput(session, "outgroup",
      label = "Select outgroup(s) for rooting",
      choices = c(tipnames()),
      selected = conroot()
    )
    session$resetBrush("plot_brush")
    og$outgroup <- conroot()
  })

  # ... or highlight pickerInput
  observeEvent(input$colbutton, {
    updatePickerInput(session, "highlight",
      choices = c(tipnames()),
      selected = unique(conroot())
    )
    session$resetBrush("plot_brush")
  })

  # create reactive values for highlighting
  highcol <- reactiveValues()

  # Whenever new tree files are loaded, reset color vector
  observeEvent(completetrees(), ignoreInit = FALSE, {
    # When rooted tree is first plotted, create a dataframe with tip names and colors (here: all black)
    # Only do this once
    observeEvent(roottree(), once = T, {
      highcol$df <- as_tibble(cbind(roottree()$tip.label, rep("black", length(roottree()$tip.label))))
    })
  })

  # If taxa are selected, update the colors for the selected taxa in the dataframe
  observeEvent(input$highlight, {
    highcol$df <- highcol$df %>%
      mutate(V2 = ifelse(V1 %in% input$highlight, input$high1, V2))
  })

  # Same when colors are changed
  observeEvent(input$high1, {
    req(input$highlight)
    highcol$df <- highcol$df %>%
      mutate(V2 = ifelse(V1 %in% input$highlight, input$high1, V2))
  })
  
  
  #| # Tab 1 (Consensus) -----
  # Plot 1 is a consensus plot

  # merge trees from all chains into single object
  treesall <- reactive({
    req(length(alltrees()) >= 1)
    if (length(completetrees()) > 1) {
      req(length(input$whichtree) > 0)
    }
    combine.trees(thintrees(), input$conburnin)
  })
  
  # calculate consensus tree
  contree.p <- reactive({
    calc.cons(treesall())
  })
  
  # collapse nodes lower than threshold chosen by user
  contree <- reactive({
    collapse.nodes(contree.p(), input$postprop)
  })

  # root the tree
  roottree <- reactive({
    
    # require these for interactive rooting
    input$reroot
    input$midpoint
    input$unroot
    
    # only update when button is pressed
    isolate(outgroup <- og$outgroup)
    
    # root the tree
    root.tree(contree(), outgroup)
  })

  # now render consensus plot
  consensusplot <- reactive({
    req(highcol$df)
    req(roottree())
    render.contree(roottree(), highcol$df, thintrees(), input$treecex, treeopts(), treefont(), input$annot)
    recordPlot()
  })
  output$consensusPlot <- renderPlot({
    consensusplot()
  })
  
  # render consensus plot with spinner
  output$consensusPlot.ui <- renderUI({
    withSpinner(plotOutput("consensusPlot",
      height = treeheight(),
      # brush here enables interactive selection of taxa
      brush = brushOpts(
        id = "plot_brush", # this enables selecting
        resetOnNew = TRUE,
        delay = 2000,
        direction = "y", # selection using y axis only
        fill = "#ccc",
        delayType = "debounce"
      ),
      width = treewidth() * 0.75
    ),
    color = "#2C4152",
    size = 0.5
    )
  })

  # PDF download handle for consensus plot
  output$consensusPDF <- downloadHandler(
    filename = "consensus.pdf",
    contentType = "application/pdf",
    content = function(file1) {
      pdf(file1, height = input$treeheight / 72, width = input$treewidth / 72)
      replayPlot(consensusplot())
      dev.off()
    }
  )
  # also export the tree as newick if wanted
  output$newick <- downloadHandler(
    filename = "consensus.nwk",
    content = function(file4) {
      write.tree(roottree(),
        file = file4, append = FALSE,
        digits = 10, tree.names = FALSE
      )
    }
  )

  #| # Tab 2 (Trees) -----

  # Slider that determines the tree generation currently displayed
  output$treegens <- renderUI({
    req(treefile$datapath)
    sliderInput("generation", "Tree generation",
      min = 1,
      max = trunc(ngentree() / treethin()),
      value = 1,
      step = 1,
      ticks = FALSE,
      animate = animationOptions(
        interval = 1000, # this adds a button that animates an iteration through all tree generations
        loop = TRUE,
        # this make the buttons pretty
        playButton = tags$button("PLAY", style = "border:2px solid; border-radius: 4px; margin:5px; padding:2px 10px; font-size:90%; background-color:white; color:black; font-weight:bold"),
        pauseButton = tags$button("PAUSE", style = "border:2px solid; border-radius: 4px; margin:5px; padding:2px 10px; font-size:90%; background-color:white; color:black; font-weight:bold")
      )
    )
  })

  # This plot shows a single tree per iteration and chain
  treeplot <- reactive({
    
    # require these for interactive rooting
    input$reroot
    input$midpoint
    input$unroot
    
    # Make sure trees for plotting are selected
    if (length(completetrees()) > 1) {
      req(length(input$whichtree) > 0)
    }
    
    # get outgroup from ui, but isolate so that rerooting is only done when button is pressed
    isolate(outgroup <- og$outgroup)
    
    # plot tree
    render.singletrees(thintrees(), outgroup, input$generation, highcol$df, input$treecex, treefont(), treeopts(), input$whichtree)
    recordPlot()
  })
  output$treeplot <- renderPlot({
    treeplot()
  })

  # render the plot with spinner & using the height and widths from ui
  output$treePlot.ui <- renderUI({
    req(treefile$datapath)
    withSpinner(plotOutput("treeplot",
      height = treeheight(),
      width = treewidth()
    ),
    color = "#2C4152",
    size = 0.5
    )
  })
  
  # PDF download handle for single trees plot
  output$singletreePDF <- downloadHandler(
    filename = "single_tree.pdf",
    contentType = "application/pdf",
    content = function(file2) {
      pdf(file2, height = input$treeheight / 72, width = input$treewidth / 72)
      replayPlot(treeplot())      
      dev.off()
    }
  )

  #| # Tab 3 (Difference) -----
  # Plot 3 shows 2 consensus plots with comparison

  differenceplot <- reactive({
    
    # require these for interactive rooting
    input$reroot
    input$midpoint
    input$unroot
    
    # check if at least 2 tree files are present and if rooting other than midpoint is selected. Throw error if not.
    isolate(outgroup <- og$outgroup)
    render.treediff(thintrees(), alltrees(), outgroup, input$conburnin, highcol$df, input$treecex, treefont(), input$whichtree)
    recordPlot()
  })
  output$differencePlot <- renderPlot({
    differenceplot()
    })
  
  # PDF download handle for difference plot
  output$differencePDF <- downloadHandler(
    filename = "difference.pdf",
    contentType = "application/pdf",
    content = function(file4) {
      pdf(file4, height = input$treeheight / 72, width = input$treewidth / 72)
      replayPlot(differenceplot())      
      dev.off()
    }
  )
  
  # finally, render this plot with spinner
  output$differencePlot.ui <- renderUI({
    shinycssloaders::withSpinner(plotOutput("differencePlot", height = treeheight(), width = treewidth()), color = "#2C4152", size = 0.5)
  })


  #| # Tab 4 (Pairwise Robinson-Foulds) -----
  # Plot 4 shows differences between trees across iterations and between chains

  # Calculate RF distances in parallel if multiple cores are available

  rP <- reactive({
    req(treefile$datapath)
    if (length(completetrees()) > 1) {
      req(length(input$whichtree) > 0)
    }
    render.rfplot(thintrees(), input$conburnin, input$whichtree, input$treecex, treescalefactor())
  })
  output$rfPlot <- renderPlot({
    rP()
  })

  output$rfPlot.ui <- renderUI({
    req(treefile$datapath)
    shinycssloaders::withSpinner(plotOutput("rfPlot", height = treeheight(), width = treewidth()),
      color = "#2C4152", size = 0.5
    )
  })

  # create download button for plot
  output$downloadrfplot <- downloadHandler(
    filename = "RFplot.pdf",
    content = function(file5) {
      pdf(file5, height = input$treeheight / 72, width = input$treewidth / 72)
      gridExtra::grid.arrange(rP(), ncol = 1)
      dev.off()
    }
  )


  #| # Tab 5 (Bipartition support) -----

  # Tab 5 counts the number of trees in which a particular group was monophyletic. It also prints a simplified tree of this group.

  # Picker input to chose which taxa to check (similar to highlight color choser)
  output$bpselect <- renderUI({
    req(treefile$datapath)
    pickerInput("bpselect",
      label = ("or select taxa"),
      choices = tipnames(),
      multiple = TRUE,
      options = list(
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} taxa",
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `width` = TRUE
      ),
      inline = FALSE,
      width = "400px"
    )
  })

  # Selection of taxa
  selectedtax <- eventReactive(input$check, {
    req(treefile$datapath)

    # Use taxa from inputpicker if nothing is entered into the text field
    if (length(input$bpselect) > 0 & input$bptext == "") {
      selectedtax <- input$bpselect
    }

    # Use taxa from text field and overwrite any potantially selected taxa in the picker
    if (input$bptext != "") {
      selectedtax <- as.character(stringr::str_split(input$bptext, "\n", simplify = TRUE))
      selectedtax <- unique(selectedtax[selectedtax != ""]) # remove empty lines and multiply selected taxa
    }

    selectedtax
  })


  # Count support for selected bipartitions
  bpsupport <- reactive({
    req(thintrees())
    calc.bpsupport(selectedtax(), tipnames(), thintrees(), input$whichtree, input$conburnin)
  })

  output$bipart <- renderText({
    req(treefile$datapath)
    req(selectedtax())
    validate(
      need(length(selectedtax()) > 1,
        message = "Please chose at least two taxa! Note that entering taxon names into the field will overwrite any selection from the drop-down menu."
      ),
      need(selectedtax() %in% unlist(tipnames()),
        message = "Taxon name(s) not in tree file. Please check."
      )
    )

    # Print result as HTML formatted text
    paste0(
      "The ", "<b>", as.character(length(selectedtax())), "</b>",
      " selected taxa are monophyletic in ",
      "<b>", as.character(bpsupport()["absolute"]), "</b>",
      " out of ",
      "<b>", as.character(bpsupport()["total"]), "</b>",
      " trees (",
      "<b>", as.character(bpsupport()["relative"]), "</b>",
      "%)."
    )
  })

  # also plot the tested groups as monophyletic. This is just for visual confirmation that the right taxa were selected.
  output$bipartPlot <- renderPlot({
    req(treefile$datapath, selectedtax(), bpsupport())
    render.bipartplot(tipnames(), selectedtax())
  })

  # render
  output$bipartPlot.ui <- renderUI({
    req(treefile$datapath, selectedtax())
    validate(
      need(length(selectedtax()) > 1, message = FALSE)
    )
    plotOutput("bipartPlot", height = treeheight(), width = treewidth())
  })


  #| # Tab 6 (RWTY) -----
  # Some useful plots from the RWTY package (https://cran.r-project.org/package=rwty)

  # convert tree files to rwty format
  rwtytrees <- reactive({
    req(treefile$datapath)
    prepare.rwty.trees(thintrees(), input$treethin)
  })
 
  # generate rwty plots
  rwtyP <- reactive({
    req(treefile$datapath)
    rwty.wrapper(input$ncores, treescalefactor(), rwtytrees())
  })
  
  # and plot 
  output$rwtyPlot <- renderPlot({
    rwtyP()
  })
  output$rwtyPlot.ui <- renderUI({
    req(treefile$datapath)
    plotOutput("rwtyPlot", height = treeheight(), width = treewidth())
  })

  output$downloadrwtyplot <- downloadHandler(
    filename = "RWTY.pdf",
    content = function(file6) {
      pdf(file6, height = input$height / 72, width = input$width / 72)
      grid.arrange(rwtyP(), ncol = 1)
      dev.off()
    }
  )
}
