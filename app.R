# get number of generations from trace files
# read all tracefiles
get.ngen <- function(tracefile) {
  tracefilelist <- lapply(tracefile$datapath, data.table::fread,
    sep = "\t",
    header = TRUE
  )
  # determine shortest trace file
  ngen <- min(unlist(lapply(tracefilelist, nrow)))
  return(ngen)
}

# read and re-format the trace files
read.trace <- function(tracefile, chain.names) {
  ({
    tracefilelist <- lapply(tracefile$datapath, data.table::fread,
      sep = "\t",
      header = TRUE,
    )
    # add chain name as parameter
    tracefilelist <- Map(cbind, tracefilelist, trace = chain.names)
    lapply(tracefilelist, tibble::as_tibble)
  })
}

# trace thinning
thin.trace <- function(trace, trace.thin) {
  trace <- trace[seq(from = 0, to = nrow(trace), by = trace.thin), ]
  trace <- select(trace, -matches("time|topo"))
}

# remove burnin from trace
burn.trace <- function(trace, burnin) {
  trace <- trace[burnin + 1:nrow(trace), ]
  return(trace)
}

# merge all trace files into dataframe
tracelist.as.df <- function(tracelist) {
  traceDF <- do.call("rbind", tracelist)
  traceDF <- tidyr::gather(traceDF, variable, value, -trace, -iter, na.rm = TRUE, factor_key = TRUE)
  return(traceDF)
}

# chose only the traces that are currently selected
choose.trace <- function(traceDF, which.chain) {
  if (length(tracefile$datapath) == 1) {
    traceDF1 <- traceDF
  }
  if (length(tracefile$datapath) > 1) {
    traceDF1 <- filter(traceDF, trace %in% which.chain)
  }
  # This adds prettier error message in case no trace file is selected
  validate(
    need(nrow(traceDF1) > 0, "Please select at least one trace file!")
  )
  return(traceDF1)
}

# set colors for trace file plotting
set.trace.colors <- function(tracedata) {
  traces <- unique(tracedata$trace)
  colorvector <- c("#377eb8", "#ff7f00", "#4daf4a", "#984ea3")
  colorvector <- colorvector[1:length(traces)]
  names(colorvector) <- traces
  return(colorvector)
}

# plotting functions
# Trace Plot
xyplot <- function(traceDF, trace.colors, trace.theme, trace.style, facet.col, cex) {
  tP1 <- ggplot(traceDF, aes(y = value, x = iter, fill = trace)) +
    facet_wrap(~variable, scales = "free", ncol = facet.col) +
    scale_color_manual(values = trace.colors)
  # This adds points to XY plots only if this option was chosen in check box
  if ("points" %in% trace.style) {
    tP1 <- tP1 + geom_point(data = traceDF, aes(y = value, x = iter, color = trace), size = cex)
  }
  # This adds lines to XY plots only if this option was chosen in check box
  if ("lines" %in% trace.style) {
    tP1 <- tP1 + geom_line(data = traceDF, aes(y = value, x = iter, color = trace), size = cex / 2)
  }
  return(tP1 + trace.theme)
}

# Violin plot
violinplot <- function(traceDF, trace.colors, trace.theme, facet.col, cex, violinplot.style) {
  vP1 <- ggplot(traceDF, aes(y = value, x = trace, fill = trace)) +
    facet_wrap(~variable, scales = "free", ncol = facet.col) +
    scale_fill_manual(values = trace.colors) +
    guides(scale_color_manual())
  # Users can add box plots and/or data points to violin plot
  # Data points must always be first layer, so each combination of violin/box plot/points is iterated below
  # No points, no box plots
  if (!"boxplot" %in% violinplot.style & !"points" %in% violinplot.style) {
    vP1 <- vP1 +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA)
  }
  # With points, no box plots
  if (!"boxplot" %in% violinplot.style & "points" %in% violinplot.style) {
    vP1 <- vP1 +
      geom_jitter(height = 0, width = 0.1, alpha = 0.2, color = "gray", show.legend = FALSE) +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA)
  }
  # No points, with box plots
  if ("boxplot" %in% violinplot.style & !"points" %in% violinplot.style) {
    vP1 <- vP1 +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
      geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = cex / 2)
  }
  # With points and box plots
  if ("points" %in% violinplot.style & "boxplot" %in% violinplot.style) {
    vP1 <- vP1 +
      geom_jitter(height = 0, width = 0.1, alpha = 0.2, color = "gray", show.legend = FALSE) +
      geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
      geom_boxplot(fill = NA, width = 0.2, color = "darkgray", outlier.shape = NA, size = cex / 2)
  }
  return(vP1 + trace.theme + theme(axis.text.x = element_blank()))
}

# density plot
densityplot <- function(traceDF, trace.colors, trace.theme, facet.col) {
  dP1 <- ggplot(traceDF) +
    geom_density(aes(x = value, fill = trace), alpha = 0.5, size = 0) +
    facet_wrap(~variable, scales = "free", ncol = facet.col) +
    scale_fill_manual(values = trace.colors) +
    trace.theme +
    theme(axis.text.y = element_blank())
  return(dP1)
}

# render plots
render.traceplots <- function(traceplot) {
  traceplot.ui <- renderUI({
    shinycssloaders::withSpinner(
    plotOutput(traceplot,
      height = input$height,
      width = input$width
    ),
    color = "#2C4152", size = 0.5
    )
  })
  return(traceplot.ui)
}

# Big function to calculate all of the summary statistics
calc.sum.stats <- function(trace.table) {
  # spread the trace file
  traceDF <- spread(trace.table, variable, value)
  tracenames <- unique(traceDF$trace)
  
  # calculate means, sds, and ess for all numeric columns
  Mean <- summarise_if(traceDF, is.numeric, mean)
  SD <- summarise_if(traceDF, is.numeric, sd)
  
  # same for ess (here, calculate ESS for each trace separately and then add values)
  ess <- group_by(traceDF, trace) %>%
    summarize_if(is.numeric, effectiveSize) %>%
    select(-trace) %>%
    summarise_all(sum)
  
  # calculate 95% HPD intervals
  hpd <- HPDinterval(mcmc(select(traceDF, -trace)), prob = 0.95)
  hpd <- as.data.frame(hpd)
  hpd <- hpd %>%
    mutate(lower = round(lower, 2)) %>%
    mutate(upper = round(upper, 2))
  names(hpd) <- c("95% HPD lower", "95% HPD upper")
  
  # calculate geweke
  geweke_res <- list()
  chainlengths <- vector()
  for (i in 1:length(tracenames)) {
    traceDF1 <- traceDF %>%
      filter(trace == tracenames[i]) %>%
      select(-trace)
    chainlengths[i] <- nrow(traceDF1)
    geweke <- coda::geweke.diag(coda::as.mcmc(traceDF1))[[1]] # first item in geweke.diag result are the actual values
    geweke <- data.frame(abs(geweke))
    names(geweke) <- paste("Geweke", tracenames[i], sep = " ")
    geweke_res[[i]] <- geweke
  }
  
  geweke_all <- do.call("cbind", geweke_res)
  min_chainlength <- min(chainlengths)
  
  # merge means, sd, ess, and geweke into results data frame
  results <- data.frame(t(rbind(Mean, SD, ess)))
  names(results) <- c("Mean", "SD", "ESS")
  results <- cbind(results, hpd)
  results <- dplyr::select(results, -SD, -ESS, SD, ESS)
  results <- cbind(results, geweke_all)
  
  # remove 'iter' & 'time' variables
  results <- results[2:nrow(results), ]
  
  # calculate discrepancy according to phylobayes manual
  if (length(tracenames) > 1) {
    # first, get means and sd for each variable and each chain separately
    means <- group_by(traceDF, trace) %>%
      summarize_if(is.numeric, mean) %>%
      select(-trace)
    sds <- group_by(traceDF, trace) %>%
      summarize_if(is.numeric, sd) %>%
      select(-trace)
    
    # then calculate discrepancy for 2 chains
    if (length(tracenames) == 2) {
      Discrepancy <- 2 * abs(means[1, ] - means[2, ]) / (sds[1, ] + sds[2, ])
      Discrepancy <- data.frame(t(Discrepancy[2:length(Discrepancy)]))
    }
    
    # for 3 & 4 chains, calculate discrepancy between any 2 chains, use average of these values
    if (length(tracenames) == 3) {
      Discrepancy <- 2 * abs(means[1, ] - means[2, ]) / (sds[1, ] + sds[2, ])
      Discrepancy[2, ] <- 2 * abs(means[1, ] - means[2, ]) / (sds[1, ] + sds[2, ])
      Discrepancy[3, ] <- 2 * abs(means[2, ] - means[3, ]) / (sds[2, ] + sds[3, ])
      Discrepancy <- dplyr::summarize_all(Discrepancy, mean)
      Discrepancy <- data.frame(t(Discrepancy[2:length(Discrepancy)]))
    }
    if (length(tracenames) == 4) {
      Discrepancy <- 2 * abs(means[1, ] - means[2, ]) / (sds[1, ] + sds[2, ])
      Discrepancy[2, ] <- 2 * abs(means[1, ] - means[3, ]) / (sds[1, ] + sds[3, ])
      Discrepancy[3, ] <- 2 * abs(means[1, ] - means[4, ]) / (sds[1, ] + sds[4, ])
      Discrepancy[4, ] <- 2 * abs(means[2, ] - means[3, ]) / (sds[2, ] + sds[3, ])
      Discrepancy[5, ] <- 2 * abs(means[2, ] - means[4, ]) / (sds[2, ] + sds[4, ])
      Discrepancy[6, ] <- 2 * abs(means[3, ] - means[4, ]) / (sds[3, ] + sds[4, ])
      Discrepancy <- dplyr::summarize_all(Discrepancy, mean)
      Discrepancy <- data.frame(t(Discrepancy[2:length(Discrepancy)]))
    }
    
    # add row name and merge with results vector
    names(Discrepancy) <- "Discrepancy"
    results <- cbind(results, Discrepancy)
    
    # create list of mcmc objects, 1 for each trace file
    tracelist <- list()
    for (j in 1:length(tracenames)) {
      tracex <- traceDF %>%
        filter(trace == tracenames[j]) %>%
        select(-trace)
      tracelist[[j]] <- coda::mcmc(tracex[1:min_chainlength, 2:ncol(tracex)])
    }
    
    # Calculate Gelman & Rubin, extract point estimates and ci, rename, and combine with results data frame
    gel.res <- coda::gelman.diag(tracelist, autoburnin = FALSE, multivariate = FALSE)
    gel.point <- data.frame(gel.res$psrf[, 1])
    gel.ci <- data.frame(gel.res$psrf[, 2])
    names(gel.point) <- "GR point estimate"
    names(gel.ci) <- "GR 95% CI"
    gel.point <- tibble::rownames_to_column(gel.point)
    gel.ci <- tibble::rownames_to_column(gel.ci)
    results <- tibble::rownames_to_column(results)
    results <- dplyr::full_join(results, gel.point, by = "rowname")
    results <- dplyr::full_join(results, gel.ci, by = "rowname")
    row.names(results) <- results$rowname
    results <- results[, 2:ncol(results)]
  }
  
  # for the numeric values in data frame, round using 2 decimals
  is.num <- sapply(results, is.numeric)
  results[is.num] <- lapply(results[is.num], round, 2)
  return(results)
}

# style data table
# call data frame with DT::data table to enable nice formatting
style.table <- function(table, trace.table) {
  styled_table <- DT::datatable(table,
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
    DT::formatStyle("ESS", color = styleInterval(99.99, c("red", "black"))) %>%
    DT::formatStyle(names(table)[grep("Geweke", names(table))], color = styleInterval(2, c("black", "red"))) %>%
    DT::formatStyle(0, fontWeight = "bold")

  if (length(unique(trace.table$trace)) > 1) {
    styled_table <- styled_table %>%
      DT::formatStyle(c("GR point estimate", "GR 95% CI"), color = styleInterval(1.2, c("black", "red"))) %>%
      DT::formatStyle("Discrepancy", color = styleInterval(0.3, c("black", "red")))
  }
  return(styled_table)
}


# set tree options and font
set.tree.opts <- function(tree.opts){
  treeopts <- vector()
  if ("align" %in% tree.opts) {
    treeopts[1] <- TRUE
  }
  if ("ignore" %in% tree.opts) {
    treeopts[2] <- FALSE
  }
  if (!("ignore" %in% tree.opts)) {
    treeopts[2] <- TRUE
  }
  if (!("align" %in% tree.opts)) {
    treeopts[1] <- FALSE
  }
  return(treeopts)
}

set.tree.font <- function(tree.font){
  treefont <- vector(mode = "numeric")
  if (is.null(tree.font)) {
    treefont <- 1
  } else {
    if (length(tree.font) == 1) {
      if (tree.font == "bold") {
        treefont <- 2
      }
      if (tree.font == "italic") {
        treefont <- 3
      }
    }
    if (length(tree.font) == 2) {
      treefont <- 4
    }
  }
  return(treefont)
}

# read trees
# only update tree format when new files are uploaded

read.treefiles <- function(treefile) {
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
  }
  return(treelist)
}

# tree thinning
thin.trees <- function(treelist, treethin) {
  thinlist <- list()
  for (i in 1:length(treelist)) {
    thinlist[[i]] <- treelist[[i]][seq(from = 1, to = length(treelist[[i]]), by = treethin)]
  }
  return(thinlist)
}

# combine all trees into single multiphylo object
combine.trees <- function(alltrees) {
  req(treefile$datapath)
  req(length(alltrees()) >= 1)
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

# Collapse nodes below a posterior probability (inspired by http://evoslav.blogspot.com/2015/01/how-to-collapse-unsupported-branches-of.html)
collapse.nodes <- function(con.tree, pp){
  contree <- con.tree
  # get position of nodes
  collapse_nodes <- which(contree$node.label < pp) + length(contree$tip.label)
  
  # get index of edges from these nodes
  collapse_indexes <- which(contree$edge[, 2] %in% collapse_nodes)
  
  # assign 0 branch length
  contree$edge.length[collapse_indexes] <- 0
  
  # use di2multi to collpase 0 branch lengths
  # important: use tiny number for tol in order for short branches with high pp not to be collapsed
  contree <- di2multi(contree, tol = 1e-10000)
  
  # remove support values for collapsed nodes
  contree$node.label[contree$node.label < pp] <- ""
  
  return(contree)
}

# Root the tree
root.tree <- function(con.tree, og){

  # unroot
  if ("<None>" %in% og) {
    con.tree <- unroot(con.tree)
  }
  
  # midpoint root
  if ("<Midpoint>" %in% og) {
    con.tree <- phytools::midpoint.root(con.tree)
  }
  
  # root with outgroup
  if (!("<Midpoint>" %in% og) & !("<None>" %in% og)) {
    validate(
      need(is.monophyletic(con.tree, og), "The specified outgroup is not monophyletic!")
    )
    con.tree <- root(con.tree,
                     outgroup = og,
                     resolve.root = TRUE
    )
  }
  
  # remove 'root label'
  con.tree$node.label[con.tree$node.label == "Root"] <- ""
  
  return(con.tree)
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
  plot(x, edge.color = treeA.cs, main = main1, tip.color = coltip1, edge.width = treeA.lw, align.tip.label = TRUE, ...)
  plot(y, edge.color = treeB.cs, main = main2, tip.color = coltip2, edge.width = treeB.lw, align.tip.label = TRUE, direction = "leftwards", ...)
  invisible()
}
