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
    lapply(tracefilelist, as.data.frame)
  })
}

# trace thinning
thin.trace <- function(trace, trace.thin) {
  trace <- trace[seq(from = 0, to = nrow(trace), by = trace.thin), ]
  trace <- select(trace, -matches("time|topo"))
  return(trace)
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

# set tree display options
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

# set tree label font
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

# create a dataframe of tip labels
# first get the tip labels and order them according to their appearance in plot (1= bottom taxon, length(tiplabels)=top taxon)
# get this info from tree$edge[,2] all numbers < ntaxa(tree) correspond to tip labels, order is as plotted
get.tip.df <- function(tree){
  # which edges belong to tips?
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  
  # order according
  ordered_tips <- tree$edge[is_tip, 2]
  
  # now just reorder the tip labels,  and add consecutive numbering as 2nd row in that dataframe
  tips <- as.data.frame(tree$tip.label[ordered_tips])
  tiporder <- as.data.frame(1:length(tree$tip.label))
  tipDF <- as.data.frame(cbind(tips, tiporder))
  names(tipDF) <- c("tips", "tiporder")
  
  # call dataframe
  return(tipDF)
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
combine.trees <- function(alltrees, burnin) {
  #req(treefile$datapath)
  req(length(alltrees) >= 1)
  treesall <- list()
  class(treesall) <- "multiPhylo"
  for (i in 1:length(alltrees)) {
    # get trees
    trees <- alltrees[[i]]
    trees <- trees[(burnin + 1):length(trees)]
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

# calculate consensus tree
calc.cons <- function(treelist){
  # get consensus branch lengths
  contree <- phytools::consensus.edges(treelist,
                                       method = "least.squares"
  )
  
  # and count how often the nodes are present in all trees (=pp) and writes this as node labels to the tree
  sv <- prop.clades(contree, treelist)
  contree$node.label <- sv / length(treelist)
  contree$node.label <- formatC(contree$node.label, digits = 2, format = "f") # 2 decimals for pp values
  
  # adjust node labels: 1) remove "root" label
  contree$node.label[contree$node.label == "Root"] <- ""
  
  # convert to numeric
  contree$node.label <- as.numeric(contree$node.label)
  
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
    
    # determine the node number of the outgroup clade
    root.node <- getMRCA(con.tree, tip = og)
    
    # root in the middle of the edge
    edge.position <- 0.5 * con.tree$edge.length[which(con.tree$edge[, 2] == root.node)]
    
    con.tree <- reroot(con.tree, root.node, edge.position)
    
  }
  
  # remove 'root label'
  con.tree$node.label[con.tree$node.label == "Root"] <- ""
  
  con.tree <- ladderize(con.tree)
  
  return(con.tree)
}

# render plot
render.contree <- function(root.tree, high.col, thin.trees, tree.cex, tree.opts, tree.font, tree.annot){
  
  # colorvector for tips 
  col.df <- high.col %>%
    arrange(factor(V1, levels = root.tree$tip.label))
  
  concolvec <- as.vector(col.df$V2)
  
  treetot <- sum(unlist(lapply(thin.trees, length)))
  
  # plot
  plot(root.tree,
       main = paste0(
         "Consensus of ", treetot - (input$conburnin * length(thin.trees)), " trees", # number of trees
         " (", length(thin.trees), " chains)"
       ), # chains
       cex.main = tree.cex * 1.1,
       cex = tree.cex,
       align.tip.label = tree.opts[1],
       use.edge.length = tree.opts[2],
       edge.width = tree.cex,
       label.offset = 0.01,
       font = tree.font,
       tip.color = concolvec
  )
  nodelabels(
    text = root.tree$node.label, # some formatting for the pp values
    frame = "none",
    adj = 0,
    cex = tree.cex * 0.8
  )
  add.scale.bar(lwd = tree.cex)
  
  # add custom plot annotations
  title(
    sub = tree.annot,
    adj = 0,
    line = 1,
    font = 2,
    cex.sub = tree.cex * 0.85
  )
  
}

render.singletrees <- function(thin.trees, og, tree.generation, high.col, tree.cex, tree.font, tree.opts, which.tree){
  
  validate(
    need(length(which.tree) > 1,
         message = "Provide tree files from at least 2 chains to display differences between consensus trees."
    )
  )
  
  # Adjusting plot layout 
  # change plot layout to 2 columns if 2 treefiles are present
  if (length(thin.trees) == 2) {
    par(mfrow = c(1, 2))
  }
  
  # 3 columns if 3 treefiles are present
  if (length(thin.trees) == 3) {
    par(mfrow = c(1, 3))
  }
  
  # and 2x2 for 4 treefiles
  if (length(thin.trees) == 4) {
    par(mfrow = c(2, 2))
  }
  
  # plot all trees in loop
  for (i in 1:length(thin.trees)) {
    
    # get trees
    trees <- thin.trees[[i]]
    
    # unroot tree if no root was chosen
    if ("<None>" %in% og) {
      currtree <- unroot(trees[[tree.generation]])
    }
    
    # midpoint root if chosen
    if ("<Midpoint>" %in% og) {
      currtree <- phytools::midpoint.root(trees[[tree.generation]])
    }
    
    # root with outgroup if chosen
    if (!("<Midpoint>" %in% og) & !("<None>" %in% og)) {
      validate(
        need(is.monophyletic(trees[[tree.generation]], og), "The specified outgroup is not monophyletic!")
      )
      # currtree <- root(trees[[tree.generation]],
      #                  outgroup = og,
      #                  resolve.root = TRUE
      # )
      
      root.node <- getMRCA(trees[[tree.generation]], tip = og)
      edge.position <- 0.5 * trees[[tree.generation]]$edge.length[which(trees[[tree.generation]]$edge[, 2] == root.node)]
      currtree <- reroot(trees[[tree.generation]], root.node, edge.position)
    }
    
    # Create tip label color vector from highlight picker options (need to do this after rooting, as this impacts tip label order)
    col.df <- high.col %>%
      arrange(factor(V1, levels = currtree$tip.label))
    concolvec <- as.vector(col.df$V2)
    
    # Plot
    plot(ladderize(currtree),
         main = paste(which.tree[i], "iteration", tree.generation),
         cex = tree.cex,
         cex.main = tree.cex * 1.1,
         edge.width = tree.cex,
         align.tip.label = tree.opts[1],
         use.edge.length = tree.opts[2],
         label.offset = 0.01,
         font = tree.font,
         tip.color = concolvec
    )
    add.scale.bar(lwd = tree.cex)
  }
}

render.treediff <- function(thin.trees, all.trees, og, burnin, high.col, tree.cex, tree.font, which.tree){
  
  validate(
    need(!("<Midpoint>" %in% og),
    message = "Midpoint rooting not possible for consensus cladograms.\nPlease choose different root!"
  ))
  
  # loop through all trees and create tree list as well as color list
  contrees <- list()
  colvecs <- list()
  
  for (i in 1:length(thin.trees)) {
    tree <- thin.trees[[i]]
    tree <- tree[(burnin + 1):length(tree)]
    contree <- consensus(tree, p = 0.5)
    
    # reroot outgroup
    if (!("<Midpoint>" %in% og) & !("<None>" %in% og)) {
      contree <- root(contree, outgroup = og)
    }
    
    # unroot tree
    if ("<None>" %in% og) {
      contree <- unroot(contree)
    }
    
    contrees[[i]] <- contree
    
    # colorvector
    
    col.df <- high.col %>%
      arrange(factor(V1, levels = contree$tip.label))
    colvecs[[i]] <- as.vector(col.df$V2)
  }
  
  
  # change plot layout according to number of chains analysed
  if (length(all.trees) == 2) {
    par(mfrow = c(1, 2))
  }
  if (length(all.trees) == 3) {
    par(mfrow = c(3, 2))
  }
  if (length(all.trees) == 4) {
    par(mfrow = c(3, 4))
  }
  
  # this nested loop plots consensus trees side by side, for all combinations in 2, 3, or 4 trees
  for (i in 1:length(thin.trees)) {
    for (j in i:length(thin.trees)) {
      if (i != j) {
        phylo.diff.new(contrees[[i]], contrees[[j]],
                       cex = tree.cex,
                       cex.main = tree.cex * 1.1,
                       label.offset = 0.01,
                       font = tree.font,
                       main1 = which.tree[i],
                       main2 = which.tree[j],
                       coltip1 = colvecs[[i]],
                       coltip2 = colvecs[[j]]
        )
      }
    }
  }
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

# Plot to display RF distances within and between chains
render.rfplot <- function(thin.trees, burnin, which.tree, tree.cex, tree.scalefactor){
  rflist <- list()
  for (i in 1:length(thin.trees)) {
    tree <- thin.trees[[i]]
    tree <- tree[(burnin + 1):length(tree)]
    rf <- vector(mode = "numeric")
    
    rf <- foreach(j = 2:length(tree), .packages = "ape", .combine = c) %dopar% {
      rf[j - 1] <- dist.topo(tree[j - 1], tree[j])
    }
    
    x <- (2 + burnin):(length(tree) + burnin)
    
    RFdf <- as.data.frame(cbind(x, rf))
    RFdf$chain <- paste("Iteration n vs. Iteration (n-1),", which.tree[i])
    rflist[[i]] <- RFdf
  }
  
  RFdf <- do.call("rbind", rflist)
  RFdf$difference <- "Tree distance within chain"
  
  # results vector
  if (length(thin.trees) > 1) {
    rflist3 <- list()
    counter <- 0
    
    for (h in 1:length(thin.trees)) {
      tree1 <- thin.trees[[h]]
      tree1 <- tree1[(burnin + 1):length(tree1)]
      for (j in h:length(thin.trees)) {
        if (h != j) {
          counter <- counter + 1
          tree2 <- thin.trees[[j]]
          tree2 <- tree2[(burnin + 1):length(tree2)]
          minlength <- min(c(length(tree1), length(tree2)))
          
          rf <- foreach(k = 1:minlength, .packages = "ape", combine = c) %dopar% {
            dist.topo(tree1[[k]], tree2[[k]])
          }
          
          x <- (1 + burnin):(length(rf) + burnin)
          
          RFdf3 <- as.data.frame(cbind(x, rf))
          RFdf3$chain <- paste(which.tree[h], which.tree[j], sep = " vs. ")
          rflist3[[counter]] <- RFdf3
        }
      }
    }
    
    RFdf3 <- do.call("rbind", rflist3)
    RFdf3 <- as.data.frame(lapply(RFdf3, unlist))
    RFdf3$difference <- "Tree distance between chains"
    
    RFdf <- rbind(RFdf, RFdf3)
  }
  
  RFggplot <- ggplot(data = RFdf, aes(x = x, y = rf, color = chain)) +
    geom_line(size = tree.cex / 2) +
    theme_light() +
    facet_wrap(~difference, nrow = 2, scales = "free_y") +
    xlab("Tree generation") +
    ylab("Robinson-Foulds distance") +
    theme(
      axis.title = element_text(size = 12 * tree.scalefactor),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.2, "cm"),
      axis.text = element_text(size = 11 * tree.scalefactor),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.text = element_text(size = 12 * tree.scalefactor),
      strip.text = element_text(size = 14 * tree.scalefactor, face = "bold")
    ) +
    scale_color_manual(values = c(
      "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
    ))
  return(RFggplot)
}

# calculate bipartition support for any clade 
calc.bpsupport <- function(selected.tax, tip.names, thin.trees, which.tree, burnin){
  validate(
    need(length(selected.tax) > 1,
         message = "Please chose at least two taxa! Note that entering taxon names into the field will overwrite any selection from the drop-down menu."
    ),
    need(length(intersect(selected.tax, tip.names)) == length(selected.tax),
         message = "Not all taxon name(s) are present in tree. Please check."
    )
  )
  
  if (length(thin.trees) > 1) {
    req(length(which.tree) > 0)
  }
  
  # read in trees in loop
  treesall <- list()
  for (i in 1:length(thin.trees)) {
    # get trees
    trees <- thin.trees[[i]]
    trees <- trees[(burnin + 1):length(trees)]
    if (i == 1) {
      treesall <- trees
    } else {
      treesall <- c(treesall, trees)
    }
  }
  
  selectedtax <- selected.tax
  
  # check if taxa from the list are monophyletic (for whole list of trees)
  bpcount <- vector()
  bpcount <- foreach(
    i = 1:length(treesall),
    .packages = "ape",
    combine = c,
    .inorder = FALSE
  ) %dopar% {
    is.monophyletic(treesall[[i]], selected.tax)
  }
  
  # count number of "TRUE"
  monotrue <- sum(unlist(bpcount))
  monotrue <- as.numeric(monotrue)
  
  # summarize results in list
  bpsupport <- list(
    absolute = monotrue, # in how many trees is the group monophyletic
    relative = formatC(monotrue / length(treesall) * 100,
                       digits = 2, format = "f"
    ), # in percent
    total = length(bpcount) # how many trees were analysed
  )
  
  return(bpsupport)
}

# simple plot of relationships tested
render.bipartplot <- function(tip.names, selected.tax){
  # get all tipnams, the ones that were selected and extract the non selected ones
  alltips <- tip.names
  selecttips <- selected.tax
  inverttips <- alltips[!alltips %in% selecttips]
  
  # this adds parentheses and commas to all selected and non-selected tip names s
  str1 <- paste("(", selecttips, "),", sep = "", collapse = "")
  str1 <- substr(str1, 1, nchar(str1) - 1)
  str2 <- paste("(", inverttips, "),", sep = "", collapse = "")
  str2 <- substr(str2, 1, nchar(str2) - 1)
  
  # to create newick tree, add paretheses around each monophyletic group and final semicolon
  bipartition <- paste0("(", str2, "),", "(", str1, ");")
  
  # read tree string
  bipartition_plot <- read.newick(text = bipartition)
  
  # make all edge lengths equal (improves display)
  bipartition_plot$edge.length <- rep(1, length(bipartition_plot$edge / 2))
  
  # Highlight the edges connecting selected trees in pink
  edge.highlight <- which.edge(bipartition_plot, selecttips)
  edgecolor <- rep("black", Nedge(bipartition_plot))
  edgecolor[edge.highlight] <- "#FF1493"
  
  # Also, increase edge width
  edgewidth <- rep(1, Nedge(bipartition_plot))
  edgewidth[edge.highlight] <- 2
  
  # and plot
  p <- plot(bipartition_plot,
       type = "phylogram",
       cex = input$treecex,
       edge.color = edgecolor,
       edge.width = edgewidth * input$treecex
  )
  return(p)
}

# reformat trees for rwty
prepare.rwty.trees <- function(thin.trees, treethin){
  rwtytrees <- list()
  for (i in 1:length(thin.trees)) {
    trees <- list()
    trees[[1]] <- thin.trees[[i]]
    trees[[2]] <- NULL
    trees[[3]] <- treethin
    names(trees) <- c("trees", "ptable", "gens.per.tree")
    class(trees) <- "rwty.chain"
    rwtytrees[[i]] <- trees
  }
  return(rwtytrees)
}

# wrapper for all rwty plots
rwty.wrapper <- function(ncores, tree.scalefactor, rwty.trees){
  # install if not present
  validate(
    need("rwty" %in% rownames(installed.packages()),
         message = "\nPackage rwty not found! Please install."
    )
  )
  library(rwty)
  
  # set number of processors fo rwty calculations
  rwty.processors <- ncores
  
  # set theme for all rwty plots
  theme_rwty <-
    theme_light() +
    theme(
      axis.title = element_text(size = 12 * tree.scalefactor),
      legend.title = element_blank(),
      axis.text = element_text(size = 11 * tree.scalefactor),
      legend.text = element_text(size = 12 * tree.scalefactor),
      title = element_text(size = 14 * tree.scalefactor),
      strip.text = element_text(size = 12 * tree.scalefactor, face = "bold")
    )
  
  # Autocorrelation
  if (input$rwtytype == "Autocorrelation") {
    autocorplot <- makeplot.autocorr(rwty.trees)
    autocorplot$autocorr.plot + theme_rwty
  }
  
  # Split frequencies
  else if (input$rwtytype == "Split frequencies") {
    cumsplitfreq <- makeplot.splitfreqs.cumulative(rwty.trees)
    slidesplitfreq <- makeplot.splitfreqs.sliding(rwty.trees)
    grid.arrange(cumsplitfreq$splitfreqs.cumulative.plot + theme_rwty,
                 slidesplitfreq$splitfreqs.sliding.plot + theme_rwty,
                 ncol = 1
    )
  }
  
  # Topology traces
  else if (input$rwtytype == "Topology trace") {
    topologyplot <- makeplot.topology(rwty.trees)
    trace <- topologyplot$trace.plot + theme_rwty
    dense <- topologyplot$density.plot + theme_rwty
    grid.arrange(trace, dense, ncol = 1)
  }
  
  # Tree space
  else if (input$rwtytype == "Tree space") {
    treespaceplot <- makeplot.treespace(rwty.trees)
    heat <- treespaceplot$treespace.heatmap + theme_rwty
    point <- treespaceplot$treespace.points.plot + theme_rwty
    grid.arrange(heat, point, ncol = 1)
  }
}

