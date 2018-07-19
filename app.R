#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Phylobayes trace stats"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
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
                     value = 10)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Trace", plotOutput("tracePlot", height="800px")),
                    tabPanel("Distribution", plotOutput("densePlot", height="800px")),
                    tabPanel("Summary statistics", tableOutput("table")),
                    tabPanel("About", includeMarkdown("background.md")))
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$tracePlot <- renderPlot({
     # Load packages 
     library(ggplot2)
     library(gridExtra)

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
           geom_point(data=trace1, aes_string(x=names[1], y=names[i]),
                      color='#377eb8',
                      size=0.5)+
           annotate("text", x = Inf, y = Inf, label = chainnames[[1]][1], color="#377eb8", vjust=1, hjust=1, size=3)+
           xlab("")+
           ylab(names[i])+   
           theme_light()
         traceplots[[i]]<- p1
       }
     } 
     
       if(length(input$file1$name)==2){
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
       }
       }
     
     if(length(input$file1$name)==3){
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
       }
     }
     
     if(length(input$file1$name)==4){
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
     }
     }

     newtraceplots <- list()
     for(i in 4:ncol(trace1)){
       newtraceplots[[i-3]] <- traceplots[[i]]  
     } 
     
     grid.arrange(grobs=newtraceplots, ncol = 2)
     
  })

output$densePlot <- renderPlot({
  
  # Load packages 
  library(ggplot2)
  library(gridExtra)
  
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
  
  grid.arrange(grobs=newdenseplots, ncol = 2)

  })

output$table <- renderTable(rownames=TRUE, {
  
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


} 
# Run the application 
shinyApp(ui = ui, server = server)

