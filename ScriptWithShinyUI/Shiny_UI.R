#run the app with
#library(shiny)
#runApp("filepath\\Shiny_UI.R")

library(RPPASPACE)
library(mgcv)
library(dplyr)
library(tibble)

source(file.path("R", "customRPPASPACEfunctions.R"))
source(file.path("R", "PostFittingAdj.R"))
source(file.path("R", "PreFittingDataHandling.R"))
source(file.path("R", "PostFittingOutput.R"))
source(file.path("R", "IterativeProcess.R"))

UI_element = fluidPage(
  titlePanel("Welcome to RPPA slide processing."),
  
  sidebarLayout(
    sidebarPanel(
      textInput("ProjDir",
                "Enter the Main Directory path:",
                value = "Copy folderpath here"),
      actionButton("startRun", 
                   "Initialise processing"),
      checkboxInput("doFlip",
                    "Flip the plate the other way around?",
                    value = TRUE)
      ),
    mainPanel(
      textOutput("textPanel")
      )
    )
  )

server_element = function(input, output, session) {
  
  progressTracker = reactiveValues(is_running = FALSE)
  
  observeEvent(input$startRun, {
    #take input
    txt = input$ProjDir
    booleanFlip = input$doFlip
    
    #prevent multiclick
    if (progressTracker$is_running) return()
    
    #check valid input
    if (!dir.exists(file.path(txt))) return()
    if (is.null(booleanFlip)) return()
    
    #set
    progressTracker$is_running = TRUE

    #run
    tryCatch({runMainProcess(txt, booleanFlip)},
             error=function(e) {
              message(conditionMessage(e))
              print(paste("Error in slide ", index, " ", antibody))
              NULL
    })
    
    #reset
    progressTracker$is_running = FALSE
    
    #give feedback
    output$textPanel = renderText({
      "Processing complete."
    })
  })
}

shinyApp(ui = UI_element, server = server_element)