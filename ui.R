# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"), # defines theme of app
  # Layout of app: 1 navbar with 3 tab panels. The first 2 tab panels have 1 sidebar + 1 main display (including multiple tabs each)
  # Final tab panel is markdown help file
  navbarPage(
    HTML("postpb"),
    # first tab panel for trace analyses
    tabPanel(
      "Parameters",
      fluidRow(
        column(4, wellPanel( # the following defines the elements of the sidebar
          fileInput("tracefile", "Select trace files",
            multiple = TRUE,
            accept = c(".trace", ".p")
          ),
          div(style = "font-size: 10px; padding: 0px 0px; margin-top:-1em",
            actionGroupButtons(c("exampletrace1", "exampletrace2", "exampletraceinf"), c("Load example 1", "Load example 2", "Info"),
              status = "default",
              size = "xs",
              fullwidth = TRUE,
              direction = "vertical"
            )
          ),
          hr(),
          uiOutput("burnin"),
          uiOutput("whichchain"),
          numericInput("prop", "Consider every Nth iteration",
            value = 10,
            width = "30%"
          ),
          sliderInput("cex", "Scaling factor for points and lines",
            min = 0.1,
            max = 10,
            value = 1,
            step = 0.1
          ),
          sliderInput("height", "Height of plot in pixels",
            min = 100,
            max = 5000,
            value = 1000,
            step = 100
          ),
          sliderInput("width", "Width of plot in pixels",
            min = 100,
            max = 5000,
            value = 1000,
            step = 100
          ),
          hr(),
          conditionalPanel(
            condition = "output.whichchain",
            downloadButton("downloadPDF", "Download pdf of all trace plots")
          )
        )),
        column(
          8, # the following defines the elements of the main display (4 tabs in total)
          tabsetPanel(
            type = "tabs",
            tabPanel(
              "Trace",
              prettyCheckboxGroup(
                inputId = "traceplotstyle",
                label = "",
                choices = c("lines", "points"),
                selected = "lines",
                inline = TRUE,
                status = "primary",
                icon = icon("check")
              ),
              hr(),
              uiOutput("tracePlot.ui")
            ),
            tabPanel(
              "Violin",
              prettyCheckboxGroup(
                inputId = "violinplotstyle",
                label = "",
                choices = c("boxplot", "points"),
                selected = c("boxplot"),
                inline = TRUE,
                status = "primary",
                icon = icon("check")
              ),
              hr(),
              uiOutput("violinPlot.ui")
            ),
            tabPanel(
              "Density",
              uiOutput("densePlot.ui")
            ),
            tabPanel(
              "Summary statistics",
              withSpinner(DT::dataTableOutput("table", width = "80%"),
                color = "#2C4152", size = 0.5
              ),
              conditionalPanel(
                condition = "output.table",
                br(),
                useShinyjs(),
                actionButton("explanation",
                  "Toggle explanations",
                  style = "padding:5px 10px; font-size:90%; background-color:white; color:black"
                ),
                hidden(div(id = "stats", includeMarkdown("stats.md")))
              )
            ),
            br()
          )
        )
      )
    ),
    # second tab panel for tree analyses, set up as "trace" tab panel
    tabPanel(
      "Trees",
      fluidRow(
        column(
          4,
          wellPanel(
            prettyRadioButtons(
              inputId = "treefiletype",
              label = "File format",
              choices = c("Newick (e.g., Phylobayes)", "Nexus (e.g., MrBayes)"),
              selected = "character(0)",
              inline = TRUE,
              shape = "round",
              status = "primary"
            ),
            conditionalPanel(
              "input.treefiletype == 'Newick (e.g., Phylobayes)' || input.treefiletype =='Nexus (e.g., MrBayes)'",
              fileInput("treefile", "Select tree files",
                multiple = TRUE,
                accept = c(".treelist", ".t")
              )
            ),
            div(style = "font-size: 10px; padding: 0px 0px; margin-top:-1em",
                actionGroupButtons(c("exampletree1", "exampletree2", "exampletreeinf"), 
                                   c("Load example 1", "Load example 2", "Info"),
                                   status = "default",
                                   size = "xs",
                                   fullwidth = TRUE,
                                   direction = "vertical"
                )
            ),
            hr(),
            uiOutput("conburnin"),
            uiOutput("whichtree"),
            flowLayout(
              uiOutput("ncores"),
              numericInput("treethin", "Consider every Nth tree",
                           value = 10)
            ),
            uiOutput("outgroups"),
            conditionalPanel(
              condition = "output.outgroups",
              actionGroupButtons(c("reroot", "midpoint", "unroot"),
                c("Reroot with outgroup", "Midpoint rooting", "Unroot"),
                status = "primary",
                size = "sm",
                fullwidth = TRUE
              )
            ),
            hr(),
            uiOutput("highlight"),
            uiOutput("highlight2"),
            conditionalPanel(
              condition = "output.outgroups",
              helpText(HTML("<b>NOTE: Taxa selected twice will be highlighted in <font color=\"#4daf4a\">green</font></b>"))
            ),
            hr(),
            prettyCheckboxGroup(
              inputId = "treeopts",
              label = "Tree plot options",
              choiceNames = c("Align labels", "Ignore branch lengths"),
              choiceValues = c("align", "ignore"),
              inline = TRUE,
              status = "primary",
              icon = icon("check")
            ),
            prettyCheckboxGroup(
              inputId = "treefont",
              label = "Tree labels",
              choices = c("bold", "italic"),
              inline = TRUE,
              status = "primary",
              icon = icon("check")
            ),
            sliderInput("treecex",
              "Scaling factor for labels and lines",
              min = 0.1,
              max = 10,
              value = 1,
              step = 0.1
            ),
            sliderInput("treeheight",
              "Height of plot in pixels",
              min = 100,
              max = 5000,
              value = 800,
              step = 100
            ),
            sliderInput("treewidth",
              "Width of plot in pixels",
              min = 100,
              max = 5000,
              value = 1000,
              step = 100
            ),
            conditionalPanel(
              condition = "output.conburnin",
              downloadButton("downloadtreePDF", "Download pdf of current tree plot")
            )
          )
        ),
        column(
          8,
          tabsetPanel(
            type = "tabs",
            tabPanel(
              "Trees",
              uiOutput("treegens"),
              conditionalPanel(condition = "output.treegens", hr()),
              uiOutput("treePlot.ui")
            ),
            tabPanel(
              "Consensus",
              conditionalPanel(
                condition = "output.consensusPlot",
                sliderInput("postprop",
                  label = "Collapse branches with posterior probabilities lower than",
                  min = 0.5,
                  max = 1,
                  value = 0.5,
                  step = 0.05,
                  ticks = FALSE
                ),
                hr()
              ),
              uiOutput("consensusPlot.ui"),
              br(),
              conditionalPanel(
                condition = "output.consensusPlot",
                downloadButton("newick",
                  "Export tree in newick format",
                  style = "padding:5px 10px; font-size:90%; background-color:white; color:black"
                )
              ),
              # verbatimTextOutput("conroot"),
              br()
            ),
            tabPanel(
              "Difference",
              uiOutput("differencePlot.ui"),
              br(),
              helpText(HTML("Conflicting nodes between consensus trees will be highlighted in <font color=\"#FF1493\"><b>pink</b></font>."))
            ),
            tabPanel(
              "Pairwise Robinson-Foulds",
              uiOutput("rfPlot.ui"),
              br(),
              conditionalPanel(
                condition = "output.rfPlot",
                downloadButton("downloadrfplot",
                  "Download plot as pdf",
                  style = "padding:5px 10px; font-size:90%; background-color:white; color:black"
                ),
                br()
              ),
              helpText("Calculation of pairwise RF distances can be computationally intensive. Reduce the number of considered trees to speed this up.")
            ),
            tabPanel(
              "Bipartition support",
              helpText("For any group of taxa in the dataset, this will determine in how many trees the group was monophyletic."),
              br(),
              div(style = "display: inline-block;vertical-align:top", conditionalPanel(
                "output.bpselect",
                textAreaInput("bptext",
                  "Enter taxon names (1 per line)", "",
                  height = "170px",
                  width = "400px"
                )
              )),
              div(style = "display: inline-block;vertical-align:top; width: 100px;", HTML("<br>")),
              div(style = "display: inline-block; horizontal-align: middle", 
                  uiOutput("bpselect"), HTML("<br><br><br>"),
                  conditionalPanel(
                    "output.bpselect",
                    actionButton("check", "Submit", style = "float: right; padding:8px 10px; font-size:100%; background-color:white; color:black"))),
              htmlOutput("bipart"),
              uiOutput("bipartPlot.ui")
            ),
            tabPanel(
              "RWTY",
              helpText(HTML("All plotting functions are all taken from the RWTY package. Please refer to the <a href='https://cran.r-project.org/web/packages/rwty/vignettes/rwty.html'>package vignette</a> for details on how to interpret the plots.<br><br>")),
              pickerInput(
                inputId = "rwtytype", 
                label = "Select RWTY plot", 
                choices = c("None", "Autocorrelation", "Split frequencies", "Topology trace", "Tree space"), 
                multiple = FALSE,
                selected = NULL
              ),
              uiOutput("rwtyPlot.ui"),
              conditionalPanel(
                condition = "output.rwtyPlot",
                downloadButton("downloadrwtyplot",
                               "Download plot as pdf",
                               style = "padding:5px 10px; font-size:90%; background-color:white; color:black"))
              )
          )
        )
      )
    ),
    tabPanel(
      "About",
      includeMarkdown("README.md")
    )
  )
)
