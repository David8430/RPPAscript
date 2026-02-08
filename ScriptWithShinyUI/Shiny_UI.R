#run the app with
#library(shiny)
#runApp("filepath\\Shiny_UI.R")

library(tidyr)
library(RPPASPACE)
library(mgcv)
library(dplyr)
library(tibble)
library(readxl)
library(ggplot2)

source(file.path("R", "main.R"))
source(file.path("R", "customRPPASPACEfunctions.R"))
source(file.path("R", "PostFittingAdj.R"))
source(file.path("R", "PreFittingDataHandling.R"))
source(file.path("R", "PostFittingOutput.R"))
source(file.path("R", "IterativeProcess.R"))
source(file.path("R", "Visualisation.R"))

UI_element = fluidPage(
  titlePanel("Welcome to RPPA slide processing."),
  
  sidebarLayout(
    sidebarPanel(
      textInput("ProjDir",
                "Enter the Main Directory path:",
                placeholder = "Copy folderpath here"),
      actionButton("startRun", 
                   "Initialise processing"),
      checkboxInput("doFlip",
                    "Invert sample plate well assingments.",
                    value = TRUE),
      checkboxInput("spatCorr",
                    "Correct for spatial deviations.",
                    value = TRUE),
      checkboxInput("protNorm",
                    "Normalise series to the respective FCF slide.",
                    value = FALSE),
      checkboxInput("devFilt",
                    "Remove outlier points.",
                    value = TRUE),
      HTML("Explanation:<br>
           -Inversion means swapping the printed sample assignents as if the plate was inserted the
           other way around.<br>
           i.e. A1 <-> P24; A2 <-> P23<br>
           -Spatial correction means adjusting fluorescence values based on bulk intensity shifts
           correlated with the physical (X,Y) coordinates.<br>
           -Outliers are defined as:<br>
           *symmetric pair dot difference >12000<br>
           NOTE: Requires the slide to have identically printed top and bottom blocks<br>
           *deviation from response curve greater than 2.5 times the respective pair<br>
           *deviation from response curve greater than 10000<br>
           -FCF normalisation will remove any samples lacking either data points.")
      ),
    mainPanel(
      textOutput("textPanel")
      )
    )
  )

server_element = function(input, output, session) {
  textPanel = reactiveValues(text = "Waiting for user input.")
  progressTracker = reactiveValues(is_running = FALSE)
  output$textPanel = renderText(textPanel$text)
    
  observeEvent(input$startRun, {
    #take input
    txt = input$ProjDir
    booleanFlip = input$doFlip
    booleanNorm = input$protNorm
    booleanSpat = input$spatCorr
    booleanFilt = input$devFilt
    
    #prevent multiclick
    if (progressTracker$is_running) return()
    
    #check valid input
    if (!dir.exists(file.path(txt))) {
      textPanel$text = "Error in directory input."
      return()
    }
    if (is.null(booleanFlip)) {
      textPanel$text = "Error in flipping checkbox input."
      return()
    }
    if (is.null(booleanNorm)) {
      textPanel$text = "Error in protein normalization input."
      return()
    }
    if (is.null(booleanSpat)) {
      textPanel$text = "Error in spatial correction input."
      return()
    }
    if (is.null(booleanFilt)) {
      textPanel$text = "Error in outlier removal input."
      return()
    }
    
    #set
    progressTracker$is_running = TRUE
    textPanel$text = "starting" #will not appear, only makes the initial text gray
    showNotification("Processing has started. Keep track of the progress in the R console.",
                     id = "prog_not",
                     type = "message",
                     duration = NULL)

    #run
    tryCatch({runMainProcess(txt, booleanFlip, booleanNorm, booleanSpat, booleanFilt)},
             error=function(e) {
              message(conditionMessage(e))
              print(paste("Error in processing."))
              textPanel$text = "Processing has encountered an unexpected error."
              NULL
    })
    
    #reset
    progressTracker$is_running = FALSE
    
    #give feedback
    removeNotification(id = "prog_not")
    showNotification("Done!",
                     type = "message",
                     duration = 3)
    textPanel$text = "Processing complete."
  })
}

shinyApp(ui = UI_element, server = server_element)