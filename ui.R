library(shiny)
library(shinythemes)
library(bslib)
library(htmltools)

# Define UI for application that draws a histogram
page_sidebar(
  
  #define the theme
  theme = shinytheme("sandstone"),
  
  # Application title
  title = "Transcriptomic Data Meta Analyst",
  
  #provide choices in the side par
  sidebar = sidebar(
    selectInput(
      inputId = "analysis_type",
      label = "Choose the type of analysis:",
      choices = c("Microarray", "RNA-Seq"),
      selected = "RNA-Seq"
    ),
    
    #add upload counts button for RNA-Seq counts file
    conditionalPanel(
      condition = "input.analysis_type == 'RNA-Seq'",
      #ask for counts file and provide an example
      fileInput(
        inputId = "rnaseq_counts",
        label = "Upload the raw counts file (tab- or comma-delimited)",
        accept = c(".txt", ".csv")
      ),
      
      tags$a("Download example counts file", target = "_blank",
          href = "data/example_counts.csv",
          download = "example_counts.csv"),
      
      #ask for design file and provide an example
      tags$hr(),
      fileInput(
        inputId = "rnaseq_design",
        label = "Upload the experimental design file (tab- or comma-delimited)",
        accept = c(".txt", ".csv")
      ),
      
      p(a("Download example experimental design file", target = "_blank",
          href = "data/experimental_design.txt",
          download = "experimental_design.txt")),
    ),
    
    #add experiment placeholder for microarray
    conditionalPanel(
      condition = "input.analysis_type == 'Microarray'",
      textAreaInput(
        inputId = "array_experiments",
        label = "Type the GEO accessions of the selected experiments (one per line)",
        rows = 12,
        resize = "none"
      )
    )
  )
)
