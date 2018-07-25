#!/usr/bin/env Rscript
# load required libraries
library(ggplot2)
library(gridExtra)
library(optparse)


# turn off warnings for quiet output
options(warn=-1)

# define options
option_list = list(
  make_option(c("-b", "--burnin"), type="integer", default=NULL, 
              help="number of burnin iterations to discard [required]", metavar="Number"),
  make_option(c("-f", "--file"), action="store_true", default=FALSE,
              help="Shall plot be saved to pdf file instead of disp? [default: %default]")
); 

# parse options into list, and the into df
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# stop script if no burnin value is provided, and print help screen
if (is.null(opt$burnin)){
  print_help(opt_parser)
  stop("\nNOTE:\tThis script plots the parameters estimated during a phylobayes run.\n\tIt should be run in a directory containing 1–4 trace files created by phylobayes.\n\tRequires R packages ggplot2, optparse, and gridExtra.\n\n", call.=FALSE)
}

# define current directory as working directory 
wd <- getwd()
setwd(wd)

# search for phylobayes files in working directory, write names into vector
files <- list.files(".", pattern=".trace")

# strip file ending of trace files for plot annotation
chainnames <- strsplit(files, ".trace", fixed=T)

# read in the files (at least one should be there, up to 4 in total)
trace1<-read.table(files[1], sep="\t", header=T)
try(trace2<-read.table(files[2], sep="\t", header=T), silent=T)
try(trace3<-read.table(files[3], sep="\t", header=T), silent=T)
try(trace4<-read.table(files[4], sep="\t", header=T), silent=T)

# substract the burnin from all trace files
trace1 <- trace1[opt$burnin+1:nrow(trace1), ]
try(trace2 <- trace2[opt$burnin+1:nrow(trace2), ], silent=T)
try(trace3 <- trace3[opt$burnin+1:nrow(trace3), ], silent=T)
try(trace4 <- trace4[opt$burnin+1:nrow(trace4), ], silent=T)
# only use every 10th iteration for plotting (faster calculation & smaller output files)
  trace1 <- trace1[seq(0, nrow(trace1), 10),]
  try(trace2 <- trace2[seq(0, nrow(trace2), 10),], silent=T)
  try(trace3 <- trace3[seq(0, nrow(trace3), 10),], silent=T)
  try(trace4 <- trace4[seq(0, nrow(trace4), 10),], silent=T)


# extract the names of the phylobayes parameters from the headers of a trace file 
names <- colnames(trace1)

# plot the parameters using ggplot (using 1 geom_plot for each trace file)
# For each phylobayes variable, one plot of the trace and one plot of the density distribution is generated
# separate lists are created for trace plots and density plots
if(length(files)==4){
  cat("Found 4 trace files!\n")
  sink("/dev/null")
  traceplots <- list()
  densplots <- list() 
  for (i in 4:length(names)){
    p1 <- ggplot()+
    geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
               color='#377eb8',
               size=0.5)+
    geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
               color='#ff7f00',
               size=0.5)+
    geom_point(data=trace3, aes_string(x=names[1], y=names[i]),
                color='#4daf4a',
                size=0.5)+
    geom_point(data=trace4, aes_string(x=names[1], y=names[i]),
                color='#984ea3',
                size=0.5)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3)+
    xlab("")+
    ylab(names[i])+   
    theme_light()
    traceplots[[i]]<- p1
    
    p2 <- ggplot()+
      geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.5, size=0 )+
      geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.5, size=0)+
      geom_density(data=trace3,aes_string(names[i]), fill="#4daf4a", alpha=0.5, size=0)+
      geom_density(data=trace3,aes_string(names[i]), fill="#984ea3", alpha=0.5, size=0)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3)+
      ylab("")+
      theme_light()
    densplots[[i]] <- p2
  }
  sink()
}

# same as above, but for three trace files
if(length(files)==3){
  cat("Found 3 trace files!\n")
  sink("/dev/null")
  traceplots <- list()
  densplots <- list() 
  for (i in 4:length(names)){
    p1 <- ggplot()+
      geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                 color='#377eb8',
                 size=0.5)+
      geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
                 color='#ff7f00',
                 size=0.5)+
      geom_point(data=trace3, aes_string(x=names[1], y=names[i]),
                 color='#4daf4a',
                 size=0.5)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
      xlab("")+
      ylab(names[i])+   
      theme_light()
    traceplots[[i]]<- p1
    
    p2 <- ggplot()+
      geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
      geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.4, size=0)+
      geom_density(data=trace3,aes_string(names[i]), fill="#4daf4a", alpha=0.4, size=0)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
      ylab("")+
      theme_light()
    densplots[[i]] <- p2
  }
  sink()
}

# same as above, but for two trace files
if(length(files)==2){
  cat("Found 2 trace files!\n")
  sink("/dev/null")
  traceplots <- list()
  densplots <- list() 
  
  for (i in 4:length(names)){
    p1 <- ggplot()+
      geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                 color='#377eb8',
                 size=0.5)+
      geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
                 color='#ff7f00',
                 size=0.5)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      xlab("")+
      ylab(names[i])+   
      theme_light()
    traceplots[[i]]<- p1

    p2 <- ggplot()+
      geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
      geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.4, size=0)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      ylab("")+
      theme_light()
    densplots[[i]] <- p2
  }
  sink()
}

# same as above, but for one trace files
if(length(files)==1){
  print("Found 1 trace file!", quote=FALSE)
  sink("/dev/null")
  traceplots <- list()
  densplots <- list() 
  
  for (i in 4:ncol(trace1)){
    p1 <- ggplot()+
      geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                 color='#377eb8',
                 size=0.5)+
      xlab("")+
      ylab(names[i])+   
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      theme_light()
    traceplots[[i]]<- p1
    
    p2 <- ggplot()+
      geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.5, size=0 )+
      ylab("")+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      theme_light()
    densplots[[i]] <- p2
  }
  sink()
}

newtraceplots <- list()
for(i in 4:ncol(trace1)){
  newtraceplots[[i-3]] <- traceplots[[i]]  
} 

newdensplots <- list()
for(i in 4:ncol(trace1)){
  newdensplots[[i-3]] <- densplots[[i]]
} 



# check if option --file was used, if yes, print traceplots to pdf file
# arrange 8 traceplots in 2 columns using grid.arrange, and print to a4 sized potrait pdf
if(opt$file==TRUE){
  cat("Writing to 'plots.pdf'\n")
  sink("/dev/null")
  pdf(file= "plots.pdf", paper="a4", width=8, height=11)
  grid.arrange(grobs=newtraceplots, ncol = 2)
  #pdf(file= "density-plots.pdf", paper="a4", width=8, height=11)
  grid.arrange(grobs=newdensplots, ncol = 2)
  dev.off()
  sink()
  }


# if option --file was not used, plot using X11
if(opt$file==FALSE){
  X11(width=8, height=10)
  grid.arrange(grobs=newtraceplots, ncol = 2)
  message("Press Return for next page")
  invisible(readLines("stdin", n=1))
  grid.arrange(grobs=newdensplots, ncol = 2)
  message("Press Return to close")
  invisible(readLines("stdin", n=1))
}

