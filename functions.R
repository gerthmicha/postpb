# get number of generations from trace files
# read all tracefiles
get.ngen <- function(tracefile){
  tracefilelist <- mclapply(tracefile$datapath, read.table,
                            sep = "\t",
                            header = TRUE,
                            check.names = FALSE,
                            comment.char = "["
  )
  # determine shortest trace file
  ngen <- min(unlist(lapply(tracefilelist, nrow)))
  return(ngen)
}

# read and re-format the trace files
read.trace <- function(tracefile)({
  tracefilelist <- mclapply(tracefile$datapath, read.table,
                            sep = "\t",
                            header = TRUE,
                            check.names = FALSE,
                            comment.char = "["
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




# function that plots 2 trees next to each other and hifhlights differences
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