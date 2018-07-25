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
library(shinycssloaders)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),
   
   # Application title
   titlePanel("Phylobayes trace stats"),
   
   # Sidebar with a slider input for number of bins 
   fluidRow(
     column(4,
            wellPanel(
        fileInput("file1", "Select Phylobayes trace files",
                  multiple = TRUE,
                  accept = c(".trace")),
        hr(),
        sliderInput("burnin",
                     "Burnin [% of iterations]:",
                     min = 0,
                     max = 100,
                     value = 20),
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
        conditionalPanel(condition = 'output.tracePlot', downloadButton("downloadPDF", "Download pdf of trace plots")),
        hr(),
        includeMarkdown("background.md")
      )
     ),
      # Show a plot of the generated distribution
      column(8,
        tabsetPanel(type = "tabs",
                    tabPanel("Trace", uiOutput("tracePlot.ui")),
                    tabPanel("Distribution", uiOutput("densePlot.ui")),
                    tabPanel("Summary statistics", withSpinner(color="#0dc5c1", tableOutput("table"))))
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, server) {
  library(ggplot2)
  library(gridExtra)
  
  myPlot <- function(plotlist){
    grid.arrange(grobs=plotlist, ncol=2)
  }
  
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
     burnin <- nrow(trace1)/100*burnin 
     trace1 <- trace1[burnin+1:nrow(trace1), ]
     trace1 <- trace1[seq(0, nrow(trace1), input$prop),]
     
     # same for trace 2–4
     burnin <- input$burnin
     try(burnin <- nrow(trace2)/100*burnin, silent = T) 
     try(trace2 <- trace2[burnin+1:nrow(trace2), ], silent=T)
     try(trace2 <- trace2[seq(0, nrow(trace2), input$prop),], silent=T)
     burnin <- input$burnin
     try(burnin <- nrow(trace3)/100*burnin, silent = T) 
     try(trace3 <- trace3[burnin+1:nrow(trace3), ], silent=T)
     try(trace3 <- trace3[seq(0, nrow(trace3), input$prop),], silent=T)
     burnin <- input$burnin
     try(burnin <- nrow(trace4)/100*burnin, silent = T) 
     try(trace4 <- trace4[burnin+1:nrow(trace4), ], silent=T)
     try(trace4 <- trace4[seq(0, nrow(trace4), input$prop),], silent=T)
     
     traceplots <- list() 
     denseplots <- list()
     
     
     if(length(input$file1$name)==1){
       for (i in 4:length(names)){
         p1 <- ggplot()+
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
           xlab("")+
           ylab(names[i])+   
           theme_light()
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
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
           annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
           xlab("")+
           ylab(names[i])+   
           theme_light()
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
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
           annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
           annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
           xlab("")+
           ylab(names[i])+   
           theme_light()
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
         annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
         annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
         annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
         annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3)+
         xlab("")+
         ylab(names[i])+   
         theme_light()
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
  burnin <- nrow(trace1)/100*burnin 
  trace1 <- trace1[burnin+1:nrow(trace1), ]
  trace1 <- trace1[seq(0, nrow(trace1), input$prop),]
  
  # same for trace 2–4
  burnin <- input$burnin
  try(burnin <- nrow(trace2)/100*burnin, silent = T) 
  try(trace2 <- trace2[burnin+1:nrow(trace2), ], silent=T)
  try(trace2 <- trace2[seq(0, nrow(trace2), input$prop),], silent=T)
  burnin <- input$burnin
  try(burnin <- nrow(trace3)/100*burnin, silent = T) 
  try(trace3 <- trace3[burnin+1:nrow(trace3), ], silent=T)
  try(trace3 <- trace3[seq(0, nrow(trace3), input$prop),], silent=T)
  burnin <- input$burnin
  try(burnin <- nrow(trace4)/100*burnin, silent = T) 
  try(trace4 <- trace4[burnin+1:nrow(trace4), ], silent=T)
  try(trace4 <- trace4[seq(0, nrow(trace4), input$prop),], silent=T)
  
  denseplots <- list()
  
  if(length(input$file1$name)==1){  
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
        ylab("")+
        theme_light()
      denseplots[[i]] <- p2
    }
  }
  if(length(input$file1$name)==2){  
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.4, size=0 )+
        geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.4, size=0)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
        ylab("")+
        theme_light()
      denseplots[[i]] <- p2
    }
  }
  if(length(input$file1$name)==3){
    for (i in 4:length(names)){
      p2 <- ggplot()+
        geom_density(data=trace1,aes_string(names[i]), fill="#377eb8", alpha=0.5, size=0 )+
        geom_density(data=trace2,aes_string(names[i]), fill="#ff7f00", alpha=0.5, size=0)+
        geom_density(data=trace3,aes_string(names[i]), fill="#4daf4a", alpha=0.5, size=0)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
        annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
        ylab("")+
        theme_light()
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
      annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[2]][1], color="#ff7f00", vjust=2.5, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[3]][1], color="#4daf4a", vjust=4, hjust=1, size=3)+
      annotate("text", x = Inf, y = Inf, label = chainnames[[4]][1], color="#984ea3", vjust=5.5, hjust=1, size=3)+
      ylab("")+
      theme_light()
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
  library(markdown)
  library(LaplacesDemon)

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
  burnin <- nrow(trace1)/100*burnin 
  burnin <- burnin+1
  trace1 <- trace1[burnin:nrow(trace1), ]
  
  burnin <- input$burnin
  try( burnin <- nrow(trace2)/100*burnin )
  burnin <- burnin+1
  try( trace2 <- trace2[burnin:nrow(trace2), ] )
  try(trace1 <- rbind(trace1,trace2))
  
  burnin <- input$burnin
  try( burnin <- nrow(trace3)/100*burnin )
  burnin <- burnin+1
  try( trace3 <- trace3[burnin:nrow(trace3), ] )
  try(trace1 <- rbind(trace1,trace3))
  
  burnin <- input$burnin
  try( burnin <- nrow(trace4)/100*burnin )
  burnin <- burnin+1
  try( trace4 <- trace2[burnin:nrow(trace4), ] )
  try(trace4 <- rbind(trace1,trace4))
  
  results <- as.data.frame(ESS(trace1))
  names(results) <- "ESS" 
  results$Mean <- colMeans(trace1)
  results$SD <- sapply(trace1, sd, na.rm = TRUE)
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

} 

# Run the application 
shinyApp(ui = ui, server = server)

