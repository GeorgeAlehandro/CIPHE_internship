#Packages used in the server

library(cipheCytoExploreR)
library(flowCore)

# Adding the limits for upload files
options(shiny.maxRequestSize = 100000 * 1024 ^ 2)

server <- function(session, input, output) {
  #Function to fix the names of the uploaded files in server side Shiny
  fixUploadedFilesNames <- function(x) {
    if (is.null(x)) {
      return()
    }
    
    oldNames = x$datapath
    newNames = file.path(dirname(x$datapath),
                         x$name)
    #Renaming files
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
  }
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  #reactiveValues
  values <- reactiveValues(
    files_selected = NULL
  )
  
  
  ###Disabling button on Loading page
  shinyjs::disable('submit')
  #Enabling submit button when files are loaded
  observe({
    if (!is.null(input$input_fcs)) {
      shinyjs::enable('submit')
    }
  })
  
  ##Program behavior when FCS files are loaded
  observeEvent(input$input_fcs, {
    ##Saves on reactive values the "metadata" of input_fcs
    values$files_selected <- input$input_fcs
    ##Display rendering of the name of the files selected by the user
    output$filechosen <- renderText({
      as.character(values$files_selected$name)
    })
    #To extract all the availables choices of columns for the files uploaded
    output$column_choice <-
      renderUI({
        selectInput("column_choice", "Column to equal sample upon:", choices = colnames(read.FCS(input$input_fcs$datapath[1])))
      })
  })
  #Submit command behavior
  observeEvent(input$submit,
               
               {showNotification("Initiated equal sampling...", type = "message")
                 files_uploaded<- fixUploadedFilesNames(input$input_fcs)
                 for (file in files_uploaded$datapath){
                   original <- read.FCS(file)
                   #Extracting all the choices of indexes (factors) in the column
                   indexes <- unique(original@exprs[,input$column_choice])
                   #Minimum sample size of 10^6
                   ground_sample <- 1000000
                   #Calculating the minimum
                   for (index in indexes){
                     subsetting = subset(original@exprs, original@exprs[,input$column_choice] == index)
                     ground_sample <- min(nrow(subsetting), ground_sample)
                   }
                   #Initiating a new matrix to concatenate the subsets
                   new = NULL
                   for (index in indexes){
                     #Selecting the data equal to each index
                     subsetting = subset(original@exprs, original@exprs[,input$column_choice] == index)
                     #Transforming to matrix
                     subsetting = as.matrix(subsetting)
                     if (is.null(new)){
                       new = subsetting[sample(nrow(subsetting), ground_sample), ]
                     }
                     else{
                       new = rbind(subsetting[sample(nrow(subsetting), ground_sample), ],new)
                     }
                   }
                   #Overwriting the expression matrix with the new concatenated
                   original@exprs <- new
                   #Saving the FCS files
                   write.FCS(original, paste0('/media/data/html/OUPUT/equalsampling/',basename(file)))
                 }
                 showNotification("Job done.", type = "message")
               })
 
  #Refresh
  observeEvent(input$refresh_input, {
    #Reset of the gatingSet
    values$files_selected <- NULL
    #Reset of the text render
    output$filechosen <- renderText({
      as.character('Input refreshed')
    })
    showNotification("all_refresh input done", type = "message")
  })

  
  #session$onSessionEnded(stopApp)
}
