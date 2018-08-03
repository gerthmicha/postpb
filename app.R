#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)

# Define UI for application
ui <- fluidPage(theme = shinytheme("flatly"),
          navbarPage("Phylobayes run stats",
                     tabPanel("Traces",
                              fluidRow(
                                column(4, wellPanel(
                                  fileInput("file1", "Select Phylobayes trace files",
                                            multiple = TRUE,
                                            accept = c(".trace")),
                                  hr(),
                                  uiOutput("burnin"),
                                  numericInput("prop",
                                               "Plot every Nth iteration",
                                               value = 10,
                                               width="30%"),
                                  checkboxInput("dots", "Plot points", TRUE),
                                  checkboxInput("lines", "Plot lines", FALSE),
                                  sliderInput("cex", "Scaling factor for points and lines", min=0.1, max=10, value=1, step=0.1),
                                  sliderInput("height", "Height of plot in pixels", min=100, max=5000, value=1000, step=100),
                                  sliderInput("width", "Width of plot in pixels", min=100, max=5000, value=1000, step=100),
                                  hr(),
                                  conditionalPanel(condition = 'output.tracePlot', downloadButton("downloadPDF", "Download pdf of trace plots"))
                                )
                                ),
                                # Tabs, output
                                column(8,
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Trace", uiOutput("tracePlot.ui")),
                                                   tabPanel("Distribution", uiOutput("densePlot.ui")),
                                                   tabPanel("Summary statistics", tableOutput("table")), conditionalPanel(condition='output.table',
                                                                                                                          includeMarkdown("stats.md")))
                                )
                              )
                     ),
                     tabPanel("Trees", 
                              fluidRow(
                                column(4, 
                                       wellPanel(
                                         fileInput("treefile", "Select Phylobayes tree files",
                                                   multiple = TRUE,
                                                   accept = c(".treelist")),
                                         verbatimTextOutput('tipcolors'),
                                         hr(),
                                         uiOutput("treegens"),
                                         uiOutput("conburnin"),
                                         uiOutput("outgroups"),
                                         uiOutput("highlight"),
                                         uiOutput("highlight2"),
                                         sliderInput("treeheight", "Height of plot in pixels", min=100, max=5000, value=800, step=100),
                                         sliderInput("treewidth", "Width of plot in pixels", min=100, max=5000, value=1000, step=100),
                                         sliderInput("treecex", "Scaling factor for labels and lines", min=0.1, max=10, value=1, step=0.1))
                                ),
                                
                                column(8, 
                                       tabsetPanel(type = "tabs", 
                                                   tabPanel("Trees", uiOutput("treePlot.ui")),
                                                   tabPanel("Consensus", uiOutput("consensusPlot.ui")),
                                                   tabPanel("Difference", uiOutput("differencePlot.ui"), helpText("Conflicting nodes between consensus trees will be highlighted in red")))
                                )
                              )          
                     ),
                     tabPanel("About",includeMarkdown("README.md"))
          ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 50*1024^2)
  library(ggplot2)
  library(gridExtra)
  library(ape)
  library(phytools)
  library(markdown)
  library(LaplacesDemon)
  library(distory)

  myPlot <- function(plotlist){
    grid.arrange(grobs=plotlist, ncol=2)
  }
  
ngen  <- reactive({
    trace1 <- read.table(input$file1$datapath[1], sep="\t", header=T, comment.char = "")
    nrow(trace1)
})

  output$burnin <- renderUI({
    req(input$file1)
    sliderInput("burnin", "Burnin [# of iterations]:", min = 0, max = ngen(),  value = ngen()/5, step=100)
  })
  
  scalefactor <- reactive({(input$height+input$width)/2000})
  
 tracelist <- reactive({
     req(input$file1)
     # Remove error message without 
     if(is.null(input$file1))     return(NULL) 

     # read in trace files (at least 1, up to 4 must be present) 
     trace1 <- read.table(input$file1$datapath[1], sep="\t", header=T, comment.char = "")
     try(trace2<-read.table(input$file1$datapath[2], sep="\t", header=T, comment.char = ""), silent=T)
     try(trace3<-read.table(input$file1$datapath[3], sep="\t", header=T, comment.char = ""), silent=T)
     try(trace4<-read.table(input$file1$datapath[4], sep="\t", header=T, comment.char = ""), silent=T)
     
     # get names of categories from trace files
     names <- colnames(trace1)
     
     # get names of chains from trace file names
     chainnames <- strsplit(input$file1$name[1:length(input$file1$name)], ".trace", fixed=T)
     
     # remove burnin (reactive), and reduce by 80%
     burnin <- input$burnin
     trace1 <- trace1[burnin+1:nrow(trace1), ]
     trace1 <- trace1[seq(0, nrow(trace1), input$prop),]
     
     # same for trace 2–4
     burnin <- input$burnin
     try(trace2 <- trace2[burnin+1:nrow(trace2), ], silent=T)
     try(trace2 <- trace2[seq(0, nrow(trace2), input$prop),], silent=T)
     burnin <- input$burnin
     try(trace3 <- trace3[burnin+1:nrow(trace3), ], silent=T)
     try(trace3 <- trace3[seq(0, nrow(trace3), input$prop),], silent=T)
     burnin <- input$burnin
     try(trace4 <- trace4[burnin+1:nrow(trace4), ], silent=T)
     try(trace4 <- trace4[seq(0, nrow(trace4), input$prop),], silent=T)
     
     traceplots <- list() 
     denseplots <- list()
     
     
     if(length(input$file1$name)==1){
       for (i in 4:length(names)){
         p1 <- ggplot()+
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
           xlab("")+
           ylab(names[i])+   
           theme_light()+
           theme(axis.text = element_text(size=11*scalefactor()*0.8),
                 axis.title = element_text(size=16*scalefactor()*0.7))
         if(input$dots==TRUE){
           p1<- p1 + geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                      color='#377eb8',
                      size=input$cex)
         }
         if(input$lines==TRUE){
           p1<- p1 + geom_line(data=trace1, aes_string(x=names[1], y=names[i]),
                                color='#377eb8',
                                size=input$cex)
         }
         traceplots[[i]]<- p1
       }
     } 
     
       if(length(input$file1$name)==2){
       for (i in 4:length(names)){
         p1 <- ggplot()+
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
           annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
           xlab("")+
           ylab(names[i])+   
           theme_light()+
           theme(axis.text = element_text(size=11*scalefactor()*0.8),
                 axis.title = element_text(size=16*scalefactor()*0.7))
         if(input$dots==TRUE){
           p1<-p1 + geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                                         color='#377eb8',
                                         size=input$cex)+
             geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
                        color='#ff7f00',
                        size=input$cex)}
         if(input$lines==TRUE){
           p1 <- p1 + geom_line(data=trace1, aes_string(x=names[1], y=names[i]),
                                                     color='#377eb8',
                                                     size=input$cex)+
            geom_line(data=trace2, aes_string(x=names[1], y=names[i]),
                      color='#ff7f00',
                      size=input$cex)}
         traceplots[[i]]<- p1
       }
       }
     
     if(length(input$file1$name)==3){
       for (i in 4:length(names)){
         p1 <- ggplot()+
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
           annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
           annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3*scalefactor())+
           xlab("")+
           ylab(names[i])+   
           theme_light()+
           theme(axis.text = element_text(size=11*scalefactor()*0.8),
                 axis.title = element_text(size=16*scalefactor()*0.7))
         if(input$dots==TRUE){
           p1<- p1 + geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                      color='#377eb8',
                      size=input$cex)+
             geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
                        color='#ff7f00',
                        size=input$cex)+
             geom_point(data=trace3, aes_string(x=names[1], y=names[i]),
                        color='#4daf4a',
                        size=input$cex)
         }
         if(input$lines==TRUE){
           p1<- p1 + geom_line(data=trace1, aes_string(x=names[1], y=names[i]),
                                color='#377eb8',
                                size=input$cex)+
             geom_line(data=trace2, aes_string(x=names[1], y=names[i]),
                        color='#ff7f00',
                        size=input$cex)+
             geom_line(data=trace3, aes_string(x=names[1], y=names[i]),
                        color='#4daf4a',
                        size=input$cex)
         }
         traceplots[[i]]<- p1
       }
     }
     
     if(length(input$file1$name)==4){
     for (i in 4:length(names)){
       p1 <- ggplot()+
         annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
         annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
         annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3*scalefactor())+
         annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3*scalefactor())+
         xlab("")+
         ylab(names[i])+   
         theme_light()+
         theme(axis.text = element_text(size=11*scalefactor()*0.8),
               axis.title = element_text(size=16*scalefactor()*0.7))
       if(input$dots==TRUE){
        p1 <- p1 + geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                    color='#377eb8',
                    size=input$cex)+
           geom_point(data=trace2, aes_string(x=names[1], y=names[i]),
                      color='#ff7f00',
                      size=input$cex)+
           geom_point(data=trace3, aes_string(x=names[1], y=names[i]),
                      color='#4daf4a',
                      size=input$cex)+
           geom_point(data=trace4, aes_string(x=names[1], y=names[i]),
                      color='#984ea3',
                      size=input$cex)
       }
       if(input$lines==TRUE){
         p1 <- p1 + geom_line(data=trace1, aes_string(x=names[1], y=names[i]),
                    color='#377eb8',
                    size=input$cex)+
           geom_line(data=trace2, aes_string(x=names[1], y=names[i]),
                      color='#ff7f00',
                      size=input$cex)+
           geom_line(data=trace3, aes_string(x=names[1], y=names[i]),
                      color='#4daf4a',
                      size=input$cex)+
           geom_line(data=trace4, aes_string(x=names[1], y=names[i]),
                      color='#984ea3',
                      size=input$cex)
       }
       traceplots[[i]]<- p1
     }
     }
     newtraceplots <- list()
     for(i in 4:ncol(trace1)){
       newtraceplots[[i-3]] <- traceplots[[i]]  
     } 
    newtraceplots
  })
   
denselist <- reactive({
  req(input$file1)
  
  # Remove error message without 
  if(is.null(input$file1))     return(NULL) 
  
  # read in trace files (at least 1, up to 4 must be present) 
  trace1 <- read.table(input$file1$datapath[1], sep="\t", header=T, comment.char = "")
  try(trace2<-read.table(input$file1$datapath[2], sep="\t", header=T, comment.char = ""), silent=T)
  try(trace3<-read.table(input$file1$datapath[3], sep="\t", header=T, comment.char = ""), silent=T)
  try(trace4<-read.table(input$file1$datapath[4], sep="\t", header=T, comment.char = ""), silent=T)
  
  # get names of categories from trace files
  names <- colnames(trace1)
  
  # get names of chains from trace file names
  chainnames <- strsplit(input$file1$name[1:length(input$file1$name)], ".trace", fixed=T)
  
  # remove burnin (reactive), and reduce by 80%
  burnin <- input$burnin
  trace1 <- trace1[burnin+1:nrow(trace1), ]
  trace1 <- trace1[seq(0, nrow(trace1), input$prop),]
  
  # same for trace 2–4
  burnin <- input$burnin
  try(trace2 <- trace2[burnin+1:nrow(trace2), ], silent=T)
  try(trace2 <- trace2[seq(0, nrow(trace2), input$prop),], silent=T)
  burnin <- input$burnin
  try(trace3 <- trace3[burnin+1:nrow(trace3), ], silent=T)
  try(trace3 <- trace3[seq(0, nrow(trace3), input$prop),], silent=T)
  burnin <- input$burnin
  try(trace4 <- trace4[burnin+1:nrow(trace4), ], silent=T)
  try(trace4 <- trace4[seq(0, nrow(trace4), input$prop),], silent=T)
  
  denseplots <- list()
  
  if(length(input$file1$name)==1){  
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
        ylab("")+
        theme_light()+
        theme(axis.text = element_text(size=11*scalefactor()*0.8),
              axis.title = element_text(size=16*scalefactor()*0.7))
      denseplots[[i]] <- p2
    }
  }
  if(length(input$file1$name)==2){  
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
        geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.4, size=0)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
        annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
        ylab("")+
        theme_light()+
        theme(axis.text = element_text(size=11*scalefactor()*0.8),
              axis.title = element_text(size=16*scalefactor()*0.7))
      denseplots[[i]] <- p2
    }
  }
  if(length(input$file1$name)==3){
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.5, size=0 )+
        geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.5, size=0)+
        geom_density(data=trace3,aes_string(names[i]), fill="#4daf4a", alpha=0.5, size=0)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
        annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
        annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3*scalefactor())+
        ylab("")+
        theme_light()+
        theme(axis.text = element_text(size=11*scalefactor()*0.8),
              axis.title = element_text(size=16*scalefactor()*0.7))
      denseplots[[i]] <- p2
    }
  }
  if(length(input$file1$name)==4){
    for (i in 4:length(names)){
    p2 <- ggplot()+
      geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.5, size=0 )+
      geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.5, size=0)+
      geom_density(data=trace3,aes_string(names[i]), fill="#4daf4a", alpha=0.5, size=0)+
      geom_density(data=trace3,aes_string(names[i]), fill="#984ea3", alpha=0.5, size=0)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3*scalefactor())+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3*scalefactor())+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3*scalefactor())+
      annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3*scalefactor())+
      ylab("")+
      theme_light()+
      theme(axis.text = element_text(size=11*scalefactor()*0.8),
            axis.title = element_text(size=16*scalefactor()*0.7))
    denseplots[[i]] <- p2
    }
  }

  
  newdenseplots <- list()
  for(i in 4:ncol(trace1)){
    newdenseplots[[i-3]] <- denseplots[[i]]  
  } 
  newdenseplots
  })

plotheight <- reactive({input$height})

plotwidth <- reactive({input$width})

output$tracePlot <- renderPlot({
  myPlot(plotlist=tracelist())
})

output$densePlot <- renderPlot({
  myPlot(plotlist=denselist())
})

output$tracePlot.ui <- renderUI({
  plotOutput("tracePlot", height=plotheight(), width=plotwidth())
})

output$densePlot.ui <- renderUI({
  plotOutput("densePlot", height=plotheight(), width=plotwidth())
})

output$table <- renderTable(rownames=TRUE, {
  
  # Remove error message without 
  if(is.null(input$file1))     return(NULL) 
  
  # read in trace files (at least 1, up to 4 must be present) 
  trace1 <- read.table(input$file1$datapath[1], sep="\t", header=T, comment.char = "")
  # read in trace files (at least 1, up to 4 must be present) 
  trace1 <- read.table(input$file1$datapath[1], sep="\t", header=T, comment.char = "")
  try(trace2<-read.table(input$file1$datapath[2], sep="\t", header=T, comment.char = ""), silent=T)
  try(trace3<-read.table(input$file1$datapath[3], sep="\t", header=T, comment.char = ""), silent=T)
  try(trace4<-read.table(input$file1$datapath[4], sep="\t", header=T, comment.char = ""), silent=T)
  
  # get names of categories from trace files
  names <- colnames(trace1)
  
  # remove burnin (reactive)
  burnin <- input$burnin
  burnin <- burnin+1
  trace1 <- trace1[burnin:nrow(trace1), ]
  
  burnin <- input$burnin
  burnin <- burnin+1
  try( trace2 <- trace2[burnin:nrow(trace2), ] )
  try( reldif <- 2*abs((colMeans(trace1)-colMeans(trace2))) / (sapply(trace1, sd, na.rm = TRUE) + sapply(trace2, sd, na.rm = TRUE)))
  try(trace1 <- rbind(trace1,trace2))
  
  burnin <- input$burnin
  burnin <- burnin+1
  try( trace3 <- trace3[burnin:nrow(trace3), ] )
  try(trace1 <- rbind(trace1,trace3))
  
  burnin <- input$burnin
  burnin <- burnin+1
  try( trace4 <- trace2[burnin:nrow(trace4), ] )
  try(trace4 <- rbind(trace1,trace4))
  
  results <- as.data.frame(ESS(trace1))
  names(results) <- "ESS" 
  results$Mean <- colMeans(trace1)
  results$SD <- sapply(trace1, sd, na.rm = TRUE)
  try( results$Discrepancy <- reldif )
  tail(results, n=-3)
})

output$downloadPDF <- downloadHandler(filename = "traceplots.pdf",
  content = function(file) {
    pdf(file,  height=input$height/72, width=input$width/72)
    grid.arrange(grobs=tracelist(), ncol = 2)
    grid.arrange(grobs=denselist(), ncol = 2)
    dev.off()
  }
  
)
trees1 <- reactive({
  req(input$treefile)
  trees1 <- read.tree(input$treefile$datapath[1])
  trees1
})
trees2 <- reactive({
  req(input$treefile)
  trees2 <- read.tree(input$treefile$datapath[2])
  trees2
})

tipnames <- reactive({
  req(input$treefile)
  tree1 <- trees1()[[1]]
  tiplabels <- list()
  labels <-sort(tree1$tip.label)
  for(i in 1:length(labels)){
    tiplabels[[i]] <- labels[i]
  }
  names(tiplabels) <- tiplabels
})

output$treegens <- renderUI({
  req(input$treefile)
  sliderInput("generation", "Tree generation", min = 1, max = length(trees1()), value = 1, step = 1)
})

output$conburnin <- renderUI({
  req(input$treefile)
sliderInput("conburnin", "Burnin for consensus calculation [# of iterations]:", min = 0, max = length(trees1()),  value = length(trees1())/5)
})

output$outgroups <- renderUI({
  req(input$treefile)
  selectInput('outgroup', 'Select outgroup(s) for rooting', choices=c("<None>", "<Midpoint>", tipnames()), selected="<None>", multiple=TRUE, selectize=F, size=6)
})

output$highlight <- renderUI({
  req(input$treefile)
  selectInput('highlight', 'Select taxa to highlight [red]', choices=c("<None>", tipnames()), selected="<None>", multiple=TRUE, selectize=F, size=6)
})

output$highlight2 <- renderUI({
  req(input$treefile)
  selectInput('highlight2', 'Select taxa to highlight [blue]', choices=c("<None>", tipnames()), selected="<None>", multiple=TRUE, selectize=F, size=6)
})

treenames <- reactive({
  req(input$treefile)
  strsplit(input$treefile$name[1:length(input$treefile$name)], ".treelist", fixed=T)
  })

tipcolors <- reactive({
  req(input$treefile)
  colors1.1 <- trees1()[[input$generation]]$tip.label %in% input$highlight 
  colors1.2 <- trees1()[[input$generation]]$tip.label %in% input$highlight2
  colvec1 <- rep("black", length=length(trees1()[[input$generation]]$tip.label))
  colvec1[which(colors1.1==TRUE)] <- "#d7191c"
  colvec1[which(colors1.2==TRUE)] <- "#2c7bb6"
  try(colors2.1 <- trees2()[[input$generation]]$tip.label %in% input$highlight) 
  try(colors2.2 <- trees2()[[input$generation]]$tip.label %in% input$highlight2) 
  try(colvec2 <- rep("black", length=length(trees2()[[input$generation]]$tip.label)))
  try(colvec2[which(colors2.1==TRUE)] <- "#d7191c")
  try(colvec2[which(colors2.2==TRUE)] <- "#2c7bb6")
  cbind(colvec1,try(colvec2))
})

output$treeplot1 <- renderPlot({
  req(input$treefile)
  if(length(input$treefile$name)==2){
    par(mfrow=c(1,2))
    }
  if(input$outgroup=="<None>"){
  plot(ladderize(trees1()[[input$generation]]), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color=tipcolors()[,1])
  }
  if(input$outgroup=="<Midpoint>"){
    currtree1 <- midpoint.root(trees1()[[input$generation]])
    concolors1.1 <- currtree1$tip.label %in% input$highlight
    concolors1.2 <- currtree1$tip.label %in% input$highlight2
    concolvec1 <- rep("black", length=length(currtree1$tip.label))
    concolvec1[which(concolors1.1==TRUE)] <- "#d7191c"
    concolvec1[which(concolors1.2==TRUE)] <- "#2c7bb6"
    plot(ladderize(currtree1), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color=concolvec1)
    }
  if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
    plot(ladderize(root(trees1()[[input$generation]], outgroup=input$outgroup, resolve.root=T)), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color=tipcolors()[,1])
  }
  if(input$outgroup=="<None>"){
    try(plot(ladderize(trees2()[[input$generation]]), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, tip.color=tipcolors()[,2]))
  }
  if(input$outgroup=="<Midpoint>"){
    try(currtree2 <- midpoint.root(trees2()[[input$generation]]))
    try(concolors2.1 <- currtree2$tip.label %in% input$highlight)
    try(concolors2.2 <- currtree2$tip.label %in% input$highlight2)
    try(concolvec2 <- rep("black", length=length(currtree2$tip.label)))
    try(concolvec2[which(concolors2.1==TRUE)] <- "#d7191c")
    try(concolvec2[which(concolors2.2==TRUE)] <- "#2c7bb6")
    try(plot(ladderize(currtree2), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, tip.color=concolvec2))
    }
  if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
    try(plot(ladderize(root(trees2()[[input$generation]], outgroup=input$outgroup, resolve.root=T)), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, tip.color=tipcolors()[,2]))
  }
  })

treeheight <- reactive({input$treeheight})
treewidth <- reactive({input$treewidth})

output$treePlot.ui <- renderUI({
  req(input$treefile)
  plotOutput("treeplot1", height=treeheight(), width=treewidth())
})

output$consensusPlot <- renderPlot({
  req(input$treefile)
  if(length(input$treefile$name)==1){
    trees1 <- trees1()[2000:length(trees1())]
    contree1 <- consensus.edges(trees1,method="least.squares")
    sv <- prop.clades(contree1, trees1)
    contree1$nodelabel <- sv/length(trees1)
    concolors1.1 <- contree1$tip.label %in% input$highlight 
    concolors1.2 <- contree1$tip.label %in% input$highlight2
    concolvec1 <- rep("black", length=length(contree1$tip.label))
    concolvec1[which(concolors1.1==TRUE)] <- "#d7191c"
    concolvec1[which(concolors1.2==TRUE)] <- "#2c7bb6"
    if(input$outgroup=="<None>"){
      plot(ladderize(contree1), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1))
    }
    if(input$outgroup=="<Midpoint>"){
      contree1 <- midpoint.root(contree1)
      concolors1.1 <- contree1$tip.label %in% input$highlight 
      concolors1.2 <- contree1$tip.label %in% input$highlight2
      concolvec1 <- rep("black", length=length(contree1$tip.label))
      concolvec1[which(concolors1.1==TRUE)] <- "#d7191c"
      concolvec1[which(concolors1.2==TRUE)] <- "#2c7bb6"
      plot(ladderize(contree1), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1))
    }
    if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
      plot(ladderize(root(contree1, outgroup=input$outgroup, resolve.root=T)), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1))
    }
  }
  if(length(input$treefile$name)==2){
    par(mfrow=c(1,2))
    trees1 <- trees1()[(input$conburnin+1):length(trees1())]
    contree1 <- consensus.edges(trees1,method="least.squares")
    sv <- prop.clades(contree1, trees1)
    contree1$nodelabel <- round((sv/length(trees1)),2)
    concolors1.1 <- contree1$tip.label %in% input$highlight 
    concolors1.2 <- contree1$tip.label %in% input$highlight2
    concolvec1 <- rep("black", length=length(contree1$tip.label))
    concolvec1[which(concolors1.1==TRUE)] <- "#d7191c"
    concolvec1[which(concolors1.2==TRUE)] <- "#2c7bb6"
    if(input$outgroup=="<None>"){
      plot(ladderize(contree1), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1), cex=input$treecex)
    }
    if(input$outgroup=="<Midpoint>"){
      contree1 <- midpoint.root(contree1)
      concolors1.1 <- contree1$tip.label %in% input$highlight 
      concolors1.2 <- contree1$tip.label %in% input$highlight2
      concolvec1 <- rep("black", length=length(contree1$tip.label))
      concolvec1[which(concolors1.1==TRUE)] <- "#d7191c"
      concolvec1[which(concolors1.2==TRUE)] <- "#2c7bb6"
      plot(ladderize(contree1), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex,  tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1), cex=input$treecex)
    }
    if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
      plot(ladderize(root(contree1, outgroup=input$outgroup, resolve.root=T)), main=treenames()[[1]], cex=input$treecex, edge.width=input$treecex, tip.color = concolvec1)
      nodelabels(text=contree1$nodelabel, frame="none", adj=c(-1), cex=input$treecex)
    }
    trees2 <- trees2()[(input$conburnin+1):length(trees2())]
    contree2 <- consensus.edges(trees2,method="least.squares")
    sv <- prop.clades(contree2, trees2)
    contree2$nodelabel <-  round((sv/length(trees2)), 2)
    concolors2.1 <- contree2$tip.label %in% input$highlight 
    concolors2.2 <- contree2$tip.label %in% input$highlight2
    concolvec2 <- rep("black", length=length(contree2$tip.label))
    concolvec2[which(concolors2.1==TRUE)] <- "#d7191c"
    concolvec2[which(concolors2.2==TRUE)] <- "#2c7bb6"
    if(input$outgroup=="<None>"){
      plot(ladderize(contree2), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, direction = "leftwards", tip.color = concolvec2)
      nodelabels(text=contree2$nodelabel, frame="none", adj=c(1.5), cex=input$treecex)
    }
    if(input$outgroup=="<Midpoint>"){
      contree2 <- midpoint.root(contree2)
      concolors2.1 <- contree2$tip.label %in% input$highlight 
      concolors2.2 <- contree2$tip.label %in% input$highlight2
      concolvec2 <- rep("black", length=length(contree2$tip.label))
      concolvec2[which(concolors2.1==TRUE)] <- "#d7191c"
      concolvec2[which(concolors2.2==TRUE)] <- "#2c7bb6"
      plot(ladderize(contree2), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, direction = "leftwards", tip.color = concolvec2)
      nodelabels(text=contree2$nodelabel, frame="none", adj=c(1.5), cex=input$treecex)
    }
    if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
      plot(ladderize(root(contree2, outgroup=input$outgroup, resolve.root=T)), main=treenames()[[2]], cex=input$treecex, edge.width=input$treecex, direction = "leftwards", tip.color = concolvec2)
      nodelabels(text=contree2$nodelabel, frame="none", adj=c(1.5), cex=input$treecex)
    }
      }
})

output$consensusPlot.ui <- renderUI({
  plotOutput("consensusPlot", height=treeheight(), width=treewidth())
}) 

output$differencePlot <- renderPlot({
  req(input$treefile[[2]])
  validate(
    need(input$outgroup != "<Midpoint>", "Midpoint rooting not possible for consensus cladograms.\nPlease choose different root.")
  )
  trees1 <- trees1()[(input$conburnin+1):length(trees1())]
  contree1 <- consensus(trees1,p=0.5)
  if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
    contree1 <- root(contree1, outgroup=input$outgroup)
  }
  trees2 <- trees2()[(input$conburnin+1):length(trees2())]
  contree2 <- consensus(trees2,p=0.5)
  if(input$outgroup!="<Midpoint>" & input$outgroup!="<None>"){
    contree2 <- root(contree2, outgroup=input$outgroup)
  }
  phylo.diff.new <-function (x, y, main1, main2, ...){
    uniqT1 <- distinct.edges(x, y)
    uniqT2 <- distinct.edges(y, x)
    treeA.cs <- rep("black", dim(x$edge)[1])
    treeA.cs[uniqT1] = "red"
    treeB.cs <- rep("black", dim(y$edge)[1])
    treeB.cs[uniqT2] = "red"
    par(mfrow = c(1, 2))
    plot(x, edge.color = treeA.cs, main=main1, ...)
    plot(y, edge.color = treeB.cs, main=main2, ...)
    invisible()
  }
  phylo.diff.new(contree1, contree2, cex=input$treecex, edge.width=input$treecex, main1=treenames()[[1]], main2=treenames()[[2]])
  })
phylo.diff
output$differencePlot.ui <- renderUI({
  req(input$treefile[[2]])
  plotOutput("differencePlot", height=treeheight(), width=treewidth())
}) 
  } 

# Run the application 
shinyApp(ui = ui, server = server)

