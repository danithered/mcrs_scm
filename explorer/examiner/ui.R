library(shiny)
library(plotly)

shinyUI(fluidPage(
  # Settings
  div(
    # Menu
    div(fluidRow(
      column(4,
             h4("State"),
             uiOutput("reports"),
             textOutput("state_data"),
             dataTableOutput("table"),
             plotlyOutput("plot")
      ),
      column(6, style="position:static;background-color:lightblue",
             h4("Cell"),
             dataTableOutput("reps"),
             tabsetPanel(
               tabPanel("MFE", plotOutput("hist_mfe")),
               tabPanel("Pfold", plotOutput("hist_Pfold")),
               tabPanel("R", plotOutput("hist_R")),
               tabPanel("active sites", plotOutput("hist_no_sites")),
               tabPanel("A", plotOutput("hist_no_acts")),
               tabPanel("Pdeg", plotOutput("hist_Pdeg"))
             ),
             plotOutput("nice")
      ),
      column(2,
             h4("Replicator"),
             plotOutput("replicator"),
             textOutput("seq"),
             textOutput("str")
      )
    ))
)))
