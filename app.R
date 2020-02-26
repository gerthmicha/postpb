# get number of generations from trace files
# read all tracefiles
get.ngen <- function(tracefile){
  tracefilelist <- lapply(tracefile$datapath, data.table::fread,
                            sep = "\t",
                            header = TRUE
  )
  # determine shortest trace file
  ngen <- min(unlist(lapply(tracefilelist, nrow)))
  return(ngen)
}

# read and re-format the trace files
read.trace <- function(tracefile)({
  tracefilelist <- lapply(tracefile$datapath, data.table::fread,
                            sep = "\t",
                            header = TRUE,
  )
  # add chain name as parameter
  tracefilelist <- Map(cbind, tracefilelist, trace = chainnames())
  # apply thinning and burnin
  tracefilelist <- lapply(tracefilelist, prep.trace)
  # combine all into 1 dataframe
  traceDF <- do.call("rbind", tracefilelist)
  traceDF <- gather(traceDF, variable, value, -trace, -iter, na.rm = TRUE, factor_key = TRUE)
  return(traceDF)
})

# apply thinning and burnin to trace files 
prep.trace <- function(trace) {
  trace <- as_tibble(trace)
  # rename all first columns as 'iter'
  colnames(trace)[1] <- "iter"
  # apply thinning
  trace <- trace[seq(from = 0, to = nrow(trace), by = input$prop), ]
  # remove burnin
  trace <- trace[input$burnin + 1:nrow(trace), ]
  # remove meaningless variables (in terms of analysis here)
  trace <- select(trace, -matches("time|topo"))
  return(trace)
}

# chose only the traces that are currently selected
chose.trace <- function(traceDF){
  if (length(tracefile$datapath) == 1) {
    traceDF <- tracedata()
  }
  if (length(tracefile$datapath) > 1) {
    traceDF <- filter(tracedata(), trace %in% input$whichchain)
  }
  # This adds prettier error message in case no trace file is selected
  validate(
    need(nrow(traceDF) > 0, "Please select at least one trace file!")
  )
  return(traceDF)
}

# set colors for trace file plotting
set.trace.colors <- function(tracedata){
  traces <- unique(tracedata$trace)
  colorvector <- c("#377eb8", "#ff7f00", "#4daf4a", "#984ea3")
  colorvector <- colorvector[1:length(traces)]
  names(colorvector) <- traces
  return(colorvector)
}

# plotting functions
# Trace Plot
xyplot <- function(traceDF){
  tP1 <- ggplot(traceDF, aes(y = value, x = iter, fill = trace)) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    scale_color_manual(values = tracecolors())
  # This adds points to XY plos only if this option was chosen in checkbox
  if ("points" %in% input$traceplotstyle) {
    tP1 <- tP1 + geom_point(data = traceDF(), aes(y = value, x = iter, color = trace), size = input$cex)
  }
  # This adds lines to XY plos only if this option was chosen in checkbox
  if ("lines" %in% input$traceplotstyle) {
    tP1 <- tP1 + geom_line(data = traceDF(), aes(y = value, x = iter, color = trace), size = input$cex / 2)
  }
  return(tP1 + tracetheme())
}

# Violin plot
violinplot <- function(traceDF){
  vP1 <- ggplot(traceDF, aes(y = value, x = trace, fill = trace)) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    scale_fill_manual(values = tracecolors()) +
    guides(scale_color_manual())
  # Users can add Boxplots and/or datapoints to violin plot
  # Datapoints must always be first layer, so each combination of Violin/Boxplot/points is iterated below
  # No points, no boxplots
  if (!"boxplot" %in% input$violinplotstyle & !"points" %in% input$violinplotstyle) {
    vP1 <- vP1 +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA)
  }
  # With points, no boxplots
  if (!"boxplot" %in% input$violinplotstyle & "points" %in% input$violinplotstyle) {
    vP1 <- vP1 +
      geom_jitter(height = 0, width = 0.1, alpha = 0.2, color = "gray", show.legend = FALSE) +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA)
  }
  # No points, with boxplots
  if ("boxplot" %in% input$violinplotstyle & !"points" %in% input$violinplotstyle) {
    vP1 <- vP1 +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
      geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = input$cex / 2)
  }
  # With points and boxplots
  if ("points" %in% input$violinplotstyle & "boxplot" %in% input$violinplotstyle) {
    vP1 <- vP1 +
      geom_jitter(height = 0, width = 0.1, alpha = 0.2, color = "gray", show.legend = FALSE) +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
      geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = input$cex / 2)
  }
  return(vP1 + tracetheme() + theme(axis.text.x = element_blank()))
}

# density plot
densityplot <- function(traceDF){
  dP1 <- ggplot(traceDF) +
    geom_density(aes(x = value, fill = trace), alpha = 0.5, size = 0) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    scale_fill_manual(values = tracecolors()) +
    tracetheme() + theme(axis.text.y = element_blank())  
  return(dP1)
}

# render plots
render.traceplots <- function(traceplot){
  traceplot.ui <- renderUI({
    withSpinner(plotOutput(traceplot,
                                   height = plotheight(),
                                   width = plotwidth()
  ),
  color = "#2C4152", size = 0.5
  )
  })
  return(traceplot.ui)
}

# calculate summary statistics
sum.stats <- function(traceDF){
  
} 
  
  
# style data table
# call dataframe with DT::datatable to enable nice formatting
style.table <- function(table){
  styled_table <- datatable(table,
                      selection = "none",
                      extensions = "Buttons",
                      options = list(
                        searching = FALSE,
                        ordering = FALSE,
                        orientation = "landscape",
                        pageLength = nrow(table),
                        dom = "Bt",
                        buttons = c("copy", "csv", "print")
                      )
  ) %>%
    formatStyle("ESS", color = styleInterval(99.99, c("red", "black"))) %>%
    formatStyle(names(table)[grep("Geweke", names(table))], color = styleInterval(2, c("black", "red"))) %>%
    formatStyle(0, fontWeight = "bold")
  
  if (length(unique(traceDF()$trace)) > 1) {
    styled_table <- styled_table %>%
      formatStyle(c("GR point estimate", "GR 95% CI"), color = styleInterval(1.2, c("black", "red"))) %>%
      formatStyle("Discrepancy", color = styleInterval(0.3, c("black", "red")))
  }
  return(styled_table)
}


# read trees 
# only update tree format when new files are uploaded

read.treefiles <- function(treefile){
  
  if (example$click == 1) {
    treeformat <- "Newick (e.g., Phylobayes)"
  }
  
  if (example$click == 2) {
   treeformat <- "Nexus (e.g., MrBayes)"
  }
  
  if (example$click == 0) {
    treeformat <- input$treefiletype
  }
  
  
  treelist <- list()
  
  for (i in 1:length(treefile$datapath)) {
    treepath <- treefile$datapath[i]
    
    if (treeformat == "Newick (e.g., Phylobayes)") {
      firstline <- readLines(treepath, n = 1)
      validate(
        # add very simple check to make sure file IS NOT nexus format
        need(!is.element("#Nexus", firstline), "Error reading file(s). Please check format."),
        need(!is.element("#NEXUS", firstline), "Error reading file(s). Please check format."),
        need(!is.element("#nexus", firstline), "Error reading file(s). Please check format.")
      )
      treelist[[i]] <- ape::read.tree(treepath)
    }
    
    if (treeformat == "Nexus (e.g., MrBayes)") {
      
      # add very simple check to make sure file IS nexus format
      firstline <- readLines(treepath, n = 1)
      validate(
        need(any(firstline == c("#Nexus", "#NEXUS", "#nexus")), "Error reading file(s). Please check format.")
      )
      treelist[[i]] <- ape::read.nexus(treepath)
    }
    # thin treelist
    treelist[[i]] <- treelist[[i]][seq(from = 1, to = length(treelist[[i]]), by = input$treethin)]
  }
  return(treelist)
}

# combine all trees into single multiphylo object
combine.trees <- function(alltrees){
  treesall <- list()
  class(treesall) <- "multiPhylo"
  for (i in 1:length(alltrees)) {
    # get trees
    trees <- alltrees[[i]]
    trees <- trees[(input$conburnin + 1):length(trees)]
    if (length(treesall) >= 1) {
      treesall <- c(treesall, trees)
    }
    if (length(treesall) == 0) {
      treesall <- trees
    }
  }
  return(treesall)
}

# function that plots 2 trees next to each other and highlights differences
# modified from http://blog.phytools.org/
phylo.diff.new <- function(x, y, main1, main2, coltip1, coltip2, ...) {
  uniqT1 <- distinct.edges(x, y)
  uniqT2 <- distinct.edges(y, x)
  treeA.cs <- rep("black", dim(x$edge)[1])
  treeA.cs[uniqT1] <- "#FF1493"
  treeA.lw <- rep(input$treecex, dim(x$edge)[1])
  treeA.lw[uniqT1] <- input$treecex * 2
  treeB.cs <- rep("black", dim(y$edge)[1])
  treeB.cs[uniqT2] <- "#FF1493"
  treeB.lw <- rep(input$treecex, dim(x$edge)[1])
  treeB.lw[uniqT2] <- input$treecex * 2
  # par(mfrow = c(1, 2))
  plot(x, edge.color = treeA.cs, main = main1, tip.color = coltip1, edge.width = treeA.lw, align.tip.label = TRUE, ...)
  plot(y, edge.color = treeB.cs, main = main2, tip.color = coltip2, edge.width = treeB.lw, align.tip.label = TRUE, direction = "leftwards", ...)
  invisible()
}


