##Packages for the server

library(dplyr)
library(flowCore)
library(cipheinfinityFlow)
library(fs)
library(ggplot2)
library(ComplexHeatmap)
library(xlsx)

options(shiny.maxRequestSize = 100000 * 1024 ^ 2)
server <- function(session, input, output) {
  #Used for retrieving server-side volumes
  volumes = getVolumes()
  #Functions to retrieve all the channel markers for the FCS files and label them by user selection
  select_modify <- function(files) {
    requireNamespace("flowCore")
    representative_file <-
      read.FCS(files[1],
               truncate_max_range = FALSE,
               ignore.text.offset = TRUE)
    data_channels <-
      pData(flowCore::parameters(representative_file)[, c("name", "desc")])
    choices <- c("backbone", "exploratory", "discard")
    result <-
      data.frame(data_channels, factor = factor(
        NA,
        levels = c('backbone', 'exploratory', 'discard'),
        ordered = TRUE
      ))
    colnames(result) <- (c('name', 'desc', 'type'))
    return(result)
    
  }
  first.quartile <- function(fr){
    quantile((fr), probs = 0.25)
    
  }
  third.quartile <- function(fr){
    quantile((fr), probs = 0.75)
    
  }
  
  ninth.decile <- function(fr){
    quantile(fr, prob = 0.9, type = 5)
    
  }
  first.decile <- function(fr){
    quantile(fr, prob = 0.1, type = 5)
    
  }
  modefunc <- function(x){
    tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- NA
    return(themode)
  }
  fixUploadedFilesNames <- function(x, plate_nb = 0) {
    if (is.null(x)) {
      return()
    }
    if (plate_nb != 0){
      
      oldNames = x$datapath
      newNames = file.path(dirname(x$datapath), x$name)
    }
    else{
    oldNames = x$datapath
    plate_annot = paste('Plate', plate_nb, '_', x$name, sep = '')
    newNames = file.path(dirname(x$datapath), plate_annot)}
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
  }
  

  ##### HIDE TABS
  hideTab("antibody_selection", "background_exploratory_selection")
  hideTab("antibody_selection", "infinity_markers_selection")
  ##### PATH TO
  ## PATHS TO PRESET
  infinity_legendscreen_plate1_path_to <-
    "/media/data/html/INPUT/markers_panel_presets/infinity_isotypes_LEGENDSCREEN_plate_1.txt"
  infinity_legendscreen_plate2_path_to <-
    "/media/data/html/INPUT/markers_panel_presets/infinity_isotypes_LEGENDSCREEN_plate_2.txt"
  infinity_legendscreen_plate3_path_to <-
    "/media/data/html/INPUT/markers_panel_presets/infinity_isotypes_LEGENDSCREEN_plate_3.txt"
  ## PATH TO OUTPUT
  path_to_output <-
    "/media/data/html/OUTPUT/infinityFlow/shiny_infinity_flow"
  
  ##PATH TO BACKBONE SPECIFICATION
  path_to_backbone_specification <-
    "/media/data/html/OUTPUT/infinityFlow/backbone_exploratory_annotation.csv"
  
  ##PATH TO INTERMEDIARY FOLDER TO MERGE MULTIPLE FILES
  
  ##PATH TO SCRIPT
  path_to_script <- "/media/data/html/source/CipheInfinity/www/"
  # Cleaning the www file of the script by deleting all the files there
  do.call(file.remove, list(list.files(path_to_script, full.names = TRUE)))
  ##### GLOBAL #####################################
  
  values <- reactiveValues(
    umap_to_plot_pop = NULL,
    umap_to_plot_stats = NULL,
    backbone_specification = NULL,
    number_files_1 = NULL,
    number_files_2 = NULL,
    number_files_3 = NULL,
    path_to_infinity_markers_1_1 = NULL,
    infinity_markers_specification_1_1 = NULL,
    infinity_markers_specification_2_1 = NULL,
    infinity_markers_specification_2_2 = NULL,
    infinity_markers_specification_3_1 = NULL,
    infinity_markers_specification_3_2 = NULL,
    infinity_markers_specification_3_3 = NULL,
    infinity_markers_specification_all_wells = NULL,
    flow.frames = NULL,
    flowAI = NULL,
    flowClean = NULL,
    file_selected = NULL,
    average_number_events = NULL,
    plate_selected_1_1 = NULL,
    plate_selected_2_1 = NULL,
    plate_selected_2_2 = NULL,
    plate_selected_3_1 = NULL,
    plate_selected_3_2 = NULL,
    plate_selected_3_3 = NULL
  )
  
  ###Disabling buttons on Loading page
  shinyjs::disable('confirm')
  shinyjs::disable('confirm_infinity')
  shinyjs::disable('plot_umap')
  shinyjs::disable('plot_heatmap')
  
  ##### UPLOAD FCS PANEL ###########################
  volumes = getVolumes()
  
  #####################################
  ### For each value selected by the slider, observe a certain event
  #### Returns the destination  but should be fixed.
  ### The value should be returned, but also showed as print as a good format
  # that can also be taken by infinity flow
  observeEvent(input$number_plates, {
    if (input$number_plates == '1') {
      reset()
      observeEvent(input$file_destination_1_1, {
        files_1_1 <- fixUploadedFilesNames(input$file_destination_1_1, 1)
        values$plate_selected_1_1 <- files_1_1
        output$filechosen_1_1 <- renderText({
          # values$plate_selected_1_1$name
          as.vector(basename(values$plate_selected_1_1$datapath))
        })
      })
    }
  })
  ## values$file_selected sticks
  observeEvent(input$number_plates, {
    if (input$number_plates == '2') {
      reset()
      observeEvent(input$file_destination_2_1, {
        files_2_1 <- fixUploadedFilesNames(input$file_destination_2_1, 1)
        values$plate_selected_2_1 <- files_2_1
        output$filechosen_2_1 <- renderText({
          values$plate_selected_2_1$name
        })
      })
      
      observeEvent(input$file_destination_2_2, {
        files_2_2 <- fixUploadedFilesNames(input$file_destination_2_2, 2)
        values$plate_selected_2_2 <- files_2_2
        output$filechosen_2_2 <- renderText({
          values$plate_selected_2_2$name
        })
      })
    }
  })
  #UI Generation based on the number of input of plates selected by the user
  observeEvent(input$number_plates, {
    if (input$number_plates == '3') {
      reset()
      observeEvent(input$file_destination_3_1, {
        files_3_1 <- fixUploadedFilesNames(input$file_destination_3_1, 1)
        values$plate_selected_3_1 <- files_3_1
        
        output$filechosen_3_1 <- renderText({
          values$plate_selected_3_1$name
        })
      })
      observeEvent(input$file_destination_3_2, {
        files_3_2 <- fixUploadedFilesNames(input$file_destination_3_2, 2)
        values$plate_selected_3_2 <- files_3_2
        output$filechosen_3_2 <- renderText({
          values$plate_selected_3_2$name
        })
      })
      observeEvent(input$file_destination_3_3, {
        files_3_3 <- fixUploadedFilesNames(input$file_destination_3_3, 3)
        values$plate_selected_3_3 <- files_3_3
        output$filechosen_3_3 <- renderText({
          values$plate_selected_3_3$name
        })
      })
    }
  })
  ####Function used to reset all the inputs when moving through slider
  reset <- function() {
    values$file_selected <- NULL
  }
  #####################################SUBMIT
  observeEvent(input$submit,
               {
                 if (input$number_plates == '1') {
                   values$file_selected <- (as.vector(values$plate_selected_1_1))
                   values$number_files_1 <-
                     length(values$plate_selected_1_1)
                 }
                 if (input$number_plates == '2') {
                   values$file_selected <-
                     rbind(values$plate_selected_2_1, values$plate_selected_2_2)
                   values$number_files_1 <-
                     length(values$plate_selected_2_1)
                   values$number_files_2 <-
                     length(values$plate_selected_2_2)

                 }
                 if (input$number_plates == '3') {
                   values$file_selected <-
                     do.call(
                       "rbind",
                       list(
                         values$plate_selected_3_1,
                         values$plate_selected_3_2,
                         values$plate_selected_3_3
                       )
                     )
                   values$number_files_1 <-
                     length(values$plate_selected_3_1)
                   
                   values$number_files_2 <-
                     length(values$plate_selected_3_2)
                   
                   values$number_files_3 <-
                     length(values$plate_selected_3_3)
                   
                   
                 }
                 fcs_files <- values$file_selected$datapath
                 #Calls the infinityFlow package for background specification
                 values$backbone_specification <-
                   select_modify(values$file_selected$datapath)
                 #Render the table
                 output$selection_table <- renderRHandsontable({
                   rhandsontable(values$backbone_specification, width = 950)
                 })
                 
                 
                 #Copy the files in path_to_merge, calculates the sum of events for all the wells
                 #One iteration that does 2 jobs
                 sum_of_all_events <- 0
                 for (file in fcs_files) {
                   #fs::file_copy()
                   # file_copy(file, path_to_merge, overwrite = T)
                   #Reading the fcs file to have access into the description of total events
                   read_fcs_file <- read.FCS(file)
                   sum_of_all_events <-
                     sum_of_all_events + as.numeric(read_fcs_file@description[["$TOT"]])
                 }
                 #Calculates the average of events of all the wells
                 values$average_number_events <-
                   sum_of_all_events / length(fcs_files)

                 ##When all is set and done, enable the buttons
                 showTab("antibody_selection",
                         "background_exploratory_selection")
                 
               })
  observeEvent(input$selection_table, {
    if (!is.null(input$selection_table)) {
      values$backbone_specification <-
        as.data.frame(hot_to_r(input$selection_table))
      output$selection_table <- renderRHandsontable({
        rhandsontable(values$backbone_specification)
      })
      if (any(is.na(values$backbone_specification[, 'type']))) {
        shinyjs::disable('confirm')
      }
      
      if (!any(is.na(values$backbone_specification[, 'type']))) {
        shinyjs::enable('confirm')
      }
      
      
      
    }
  })
  ########################CONFIRM EXPLORATORY AND BACKGROUND SELECTION
  ## Checks if there are NA values --> throws an error if value is NA
  ## Then if all the values are not NA --> checks if at least one of the data
  ## acquisition is 'exploratory'
  observeEvent(input$confirm,
               {
                 if (any(is.na(values$backbone_specification[, 'type']))) {
                   showNotification("One or more type values are empty inside the table.", type = 'error')
                   
                 }
                 else{
                   if (!any(values$backbone_specification[, 'type'] == "exploratory")) {
                     showNotification("At least one measurement must be exploratory.", type = 'error')
                   }
                 }
                 ## path_to_backbone_specification <- "S:\\Mcyto\\Experiments\\R&D\\2022\\George-Alehandro_InfinityFlow\\markers_panel_presets\\test_annotation_11_4.csv"
                 write.csv(values$backbone_specification,
                           file = path_to_backbone_specification,
                           row.names = FALSE)
                 showTab("antibody_selection", "infinity_markers_selection")
               })
  ######################### ADRESSING THE INFINITY MARKERS TABLE
  observeEvent(input$infinity_markers_1_1,
               {
                 values$infinity_markers_specification_1_1 <-
                   read.csv(input$infinity_markers_1_1$datapath, sep = "\t")
                 output$infinity_markers_table_1_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_1_1,
                       width = 950,
                       height = 400
                     )
                   })
                 shinyjs::enable('confirm_infinity')
                 
               })
  
  
  observeEvent(input$infinity_markers_table_1_1, {
    values$infinity_markers_specification_1_1 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_1_1))
  })
  
  observeEvent(input$infinity_markers_preset_1_1,
               {
                 if (input$infinity_markers_preset_1_1 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_1_1 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_1_1 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_1_1 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_1_1 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_1_1 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                   
                 }
                 output$infinity_markers_table_1_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_1_1,
                       width = 950,
                       height = 400
                     )
                   })
                 shinyjs::enable('confirm_infinity')
                 
               })
  
  #################### 2 PLATES, FIRST PLATE
  
  observeEvent(input$infinity_markers_2_1,
               {
                 values$infinity_markers_specification_2_1 <-
                   read.csv(input$infinity_markers_2_1$datapath, sep = "\t")
                 # -> dataframe
                 output$infinity_markers_table_2_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_2_1,
                       width = 950,
                       height = 400
                     )
                   })
                 #Render the datatable of background and exploratory specification
                 ##LOADED
                 shinyjs::enable('confirm_infinity')
                 
               })
  observeEvent(input$infinity_markers_table_2_1, {
    values$infinity_markers_specification_2_1 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_2_1))
  })
  
  observeEvent(input$infinity_markers_preset_2_1,
               {
                 if (input$infinity_markers_preset_2_1 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_2_1 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                 }
                 if (input$infinity_markers_preset_2_1 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_2_1 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                 }
                 if (input$infinity_markers_preset_2_1 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_2_1 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                 }
                 output$infinity_markers_table_2_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_2_1,
                       width = 950,
                       height = 400
                     )
                   })
                 shinyjs::enable('confirm_infinity')
                 
               })
  #################### 2 PLATES, SECOND PLATE
  
  observeEvent(input$infinity_markers_2_2,
               {
                 values$infinity_markers_specification_2_2 <-
                   read.csv(input$infinity_markers_2_2$datapath, sep = "\t")
                 output$infinity_markers_table_2_2 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_2_2,
                       width = 950,
                       height = 400
                     )
                   })
                 
               })
  observeEvent(input$infinity_markers_table_2_2, {
    values$infinity_markers_specification_2_2 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_2_2))
  })
  
  observeEvent(input$infinity_markers_preset_2_2,
               {
                 if (input$infinity_markers_preset_2_2 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_2_2 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_2_2 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_2_2 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_2_2 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_2_2 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                 }
                 output$infinity_markers_table_2_2 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_2_2,
                       width = 950,
                       height = 400
                     )
                   })
                 
               })
  #################### 3 PLATES, FIRST PLATE
  
  observeEvent(input$infinity_markers_3_1,
               {
                 values$infinity_markers_specification_3_1 <-
                   read.csv(input$infinity_markers_3_1$datapath, sep = "\t")
                 output$infinity_markers_table_3_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_3_1,
                       width = 250,
                       height = 400
                     )
                   })
                 shinyjs::enable('confirm_infinity')
               })
  observeEvent(input$infinity_markers_table_3_1, {
    values$infinity_markers_specification_3_1 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_3_1))
  })
  
  observeEvent(input$infinity_markers_preset_3_1,
               {
                 if (input$infinity_markers_preset_3_1 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_3_1 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_1 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_3_1 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_1 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_3_1 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                 }
                 output$infinity_markers_table_3_1 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_3_1,
                       width = 250,
                       height = 400
                     )
                   })
                 shinyjs::enable('confirm_infinity')
                 
               })
  #################### 3 PLATES, SECOND PLATE
  
  observeEvent(input$infinity_markers_3_2,
               {
                 values$infinity_markers_specification_3_2 <-
                   read.csv(input$infinity_markers_3_2$datapath, sep = "\t")
                 
                 
               })
  observeEvent(input$infinity_markers_table_3_2, {
    values$infinity_markers_specification_3_2 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_3_2))
  })
  
  observeEvent(input$infinity_markers_preset_3_2,
               {
                 if (input$infinity_markers_preset_3_2 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_3_2 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_2 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_3_2 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_2 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_3_2 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                 }
                 output$infinity_markers_table_3_2 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_3_2,
                       width = 250,
                       height = 400
                     )
                   })
                 
               })
  
  #################### 3 PLATES, THIRD PLATE
  
  observeEvent(input$infinity_markers_3_3,
               {
                 values$infinity_markers_specification_3_3 <-
                   read.csv(input$infinity_markers_3_3$datapath, sep = "\t")
                 output$infinity_markers_table_3_3 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_3_3,
                       width = 250,
                       height = 400
                     )
                   })
                 
               })
  observeEvent(input$infinity_markers_table_3_3, {
    values$infinity_markers_specification_3_3 <-
      as.data.frame(hot_to_r(input$infinity_markers_table_3_3))
  })
  
  observeEvent(input$infinity_markers_preset_3_3,
               {
                 if (input$infinity_markers_preset_3_3 == "infinity_isotypes_LEGENDSCREEN_plate_1") {
                   values$infinity_markers_specification_3_3 <-
                     read.csv(infinity_legendscreen_plate1_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_2 == "infinity_isotypes_LEGENDSCREEN_plate_2") {
                   values$infinity_markers_specification_3_3 <-
                     read.csv(infinity_legendscreen_plate2_path_to, sep = "\t")
                   
                 }
                 if (input$infinity_markers_preset_3_3 == "infinity_isotypes_LEGENDSCREEN_plate_3") {
                   values$infinity_markers_specification_3_3 <-
                     read.csv(infinity_legendscreen_plate3_path_to, sep = "\t")
                   
                 }
                 output$infinity_markers_table_3_3 <-
                   renderRHandsontable({
                     rhandsontable(
                       values$infinity_markers_specification_3_3,
                       width = 250,
                       height = 400
                     )
                   })
                 
                 
               })
  
  #Functionality of the Confirm Infinity Markers button, that takes the dataframe
  #made by the different antibodies in order to infer each exploratory antibody to each fcs file
  #Resulting value is stored within the reactive value of values$infinity_markers_specification_all_wells
  observeEvent(input$confirm_infinity,
               {
                 if (input$number_plates == '1' &&
                     !is.null(values$infinity_markers_specification_1_1)) {
                   values$infinity_markers_specification_all_wells <-
                     values$infinity_markers_specification_1_1
                 }
                 if (input$number_plates == '2' &&
                     !is.null(values$infinity_markers_specification_2_1) &&
                     !is.null(values$infinity_markers_specification_2_2)) {
                   values$infinity_markers_specification_all_wells <-
                     rbind(
                       values$infinity_markers_specification_2_1,
                       values$infinity_markers_specification_2_2
                     )
                 }
                 if (input$number_plates == '3' &&
                     !is.null(values$infinity_markers_specification_3_1) &&
                     !is.null(values$infinity_markers_specification_3_2) &&
                     !is.null(values$infinity_markers_specification_3_3)) {
                   values$infinity_markers_specification_all_wells <-
                     rbind(
                       values$infinity_markers_specification_3_1,
                       values$infinity_markers_specification_3_2,
                       values$infinity_markers_specification_3_3
                     )
                   
                 }
                 #<!>
                 #Sorting files in respect to OS selection order
                 files <- gtools::mixedsort(basename(values$file_selected$datapath))
                 row.names(values$infinity_markers_specification_all_wells) = files
               })
  
  #Function to create an empty dataframe with number of rows equal to the number of files
  create_empty_df <- function(number_of_rows) {
    df <- data.frame(matrix(NA, nrow = number_of_rows, ncol = 2))
    
    colnames(df) <- c('Infinity_target', 'Infinity_isotype')
    df$Infinity_target <- as.character(df$Infinity_target)
    df$Infinity_isotype <- as.character(df$Infinity_isotype)
    return(df)
    
  }
  #Initialized customized table button behavior
  observeEvent(input$initialize_empty,
               {
                 if (input$number_plates == '1') {
                   empty_table_1 <- create_empty_df(values$number_files_1)
                   output$infinity_markers_table_1_1 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_1, width = 950, height = 400)
                     })
                   
                 }
                 if (input$number_plates == '2') {
                   empty_table_1 <- create_empty_df(values$number_files_1)
                   empty_table_2 <-
                     create_empty_df(values$number_files_2)
                   output$infinity_markers_table_2_1 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_1, width = 950, height = 400)
                     })
                   output$infinity_markers_table_2_2 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_2, width = 950, height = 400)
                     })
                   
                 }
                 if (input$number_plates == '3') {
                   empty_table_1 <- create_empty_df(values$number_files_1)
                   empty_table_2 <-
                     create_empty_df(values$number_files_2)
                   empty_table_3 <-
                     create_empty_df(values$number_files_3)
                   output$infinity_markers_table_3_1 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_1, width = 250, height = 400)
                     })
                   output$infinity_markers_table_3_2 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_2, width = 250, height = 400)
                     })
                   output$infinity_markers_table_3_3 <-
                     renderRHandsontable({
                       rhandsontable(empty_table_3, width = 250, height = 400)
                     })
                   
                 }
               })
  
  ######################### Begin pipeline
  observeEvent(input$begin_pipeline, {
    if (file.exists(path_to_output)) {
      unlink(path_to_output, recursive = TRUE)
      cat("The directory is deleted")
    }
    #Fetching the value of specified sliders
    input_events_downsampling <-
      (input$input_events_downsampling / 100) * values$average_number_events
    prediction_events_downsampling <-
      (input$prediction_events_downsampling / 100) * values$average_number_events
    cores <- input$cores_used
    if (input$transform == 'Yes') {
      transform <- F
    }
    else{
      transform <- T
    }
    targets <-
      values$infinity_markers_specification_all_wells$Infinity_target
    names(targets) <-
      rownames(values$infinity_markers_specification_all_wells)
    isotypes <-
      values$infinity_markers_specification_all_wells$Infinity_isotype
    names(isotypes) <-
      rownames(values$infinity_markers_specification_all_wells)
    
    
    ####Real beginning of the pipeline
    #to_concatenate <- c(input$cores_used, 'L')
    # cores <- paste(to_concatenate, collapse = '')
    cores <- as.numeric(input$cores_used)
    files <- values$file_selected$datapath
    files <- files[match(gtools::mixedsort(basename(files)), basename(files))]
    
    imputed_data <- infinity_flow(
      vector_of_file_absolute_paths = files,
      path_to_output = path_to_output,
      #   path_to_intermediary_results = path_to_intermediary_results,
      backbone_selection_file = path_to_backbone_specification,
      annotation = targets,
      isotype = isotypes,
      input_events_downsampling = input_events_downsampling,
      prediction_events_downsampling = prediction_events_downsampling,
      verbose = TRUE,
      cores = cores,
      transform = transform
    )
  })
  
  
  
  
  observeEvent(input$files_to_view, {
    values$choices_to_add <- input$files_to_view$name
    values$choices_to_add <-
      gsub(".*target_", 'Target ', values$choices_to_add)
    values$choices_to_add <- gsub(".fcs", '', values$choices_to_add)
    showTab("view_tab", "view_heatmap")
    showTab("view_tab", "view_umap")
  })
  
  
  
  observeEvent(input$plot_umap, {
    values$umap_to_plot_pops = NULL
    values$umap_to_plot_stats = NULL
    index_pops = NULL
    eve = NULL
    index_plots <-
      match(input$selectize_files, values$choices_to_add)
    concatenated_file <- NULL
    #If user didn't specify any particular populations to plot
    for (index in index_plots) {
      a <-
        (read.FCS(input$files_to_view$datapath[index])@exprs)
      
      if (is.null(concatenated_file)) {
        concatenated_file <- a
      }
      else{
        concatenated_file <- rbind(concatenated_file, a)
      }
      
      
    }
    if (!is.null(input$selectize_pops))
      ({
        values$pop_annot = read.csv(input$load_annot$datapath)
        index_pops <-
          match(input$selectize_pops, values$pop_annot[, 'Population'])
        concatenated_file <-
          subset(concatenated_file, subset = concatenated_file[, input$selectize_col_factor] %in% index_pops)
      })
    concatenated_file = as.data.frame(concatenated_file)
    concatenated_file[, input$selectize_col_factor] = as.factor(concatenated_file[, input$selectize_col_factor])
    values$umap_to_plot_pops <-
      ggplot(concatenated_file, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = eval(as.symbol(input$selectize_col_factor))))
    #Ploting stats on UMAP
    values$umap_to_plot_stats <-
      ggplot(concatenated_file, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(col = eval(
        as.symbol(input$selectize_cols_umap)
      ))) +
      scale_colour_gradient2(
        low = ("blue"),
        mid = "yellow",
        high = ("red"),
        midpoint = quantile(concatenated_file[, input$selectize_cols_umap], probs = 0.85) /
          2
      ) +
      geom_point(data = concatenated_file[concatenated_file[, input$selectize_cols_umap] >
                                            quantile(concatenated_file[, input$selectize_cols_umap], probs = 0.85), ], col = "red")
    #Loading annotation if inserted by the user
    if (!is.null(input$load_annot)) {
      if (!is.null(index_pops)) {
        values$pop_annot = subset(values$pop_annot, subset = values$pop_annot[, 'Index'] %in% index_pops)
      }
      legend = values$pop_annot[, 'Population']
      values$umap_to_plot_pops <- values$umap_to_plot_pops +
        scale_color_discrete(name = input$selectize_col_factor, labels = legend)
    }
    output$plot_umap_pops <-
      renderPlot({
        values$umap_to_plot_pops
      })
    output$plot_umap_stats <-
      renderPlot({
        values$umap_to_plot_stats
      })
    
    
  })
  
  #Heatmap not background corrected
  observeEvent(input$files_to_view, {
    values$choices_to_add <- input$files_to_view$name
    values$choices_to_add <-
      gsub(".*target_", 'Target ', values$choices_to_add)
    values$choices_to_add <-
      gsub(".fcs", '', values$choices_to_add)
    output$ui_files <- renderUI({
      pickerInput(
        inputId = "selectize_files",
        label = "Select files for plot",
        choices = values$choices_to_add,
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE
      )
    })
    
    
    output$ui_cols_heatmap <- renderUI({
      pickerInput(
        inputId = "selectize_cols_heatmap",
        label = 'Select markers for heatmap plot',
        choices = as.vector(colnames(
          read.FCS(input$files_to_view$datapath[1])@exprs
        )),
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE
      )
    })
    output$ui_factor <- renderUI({
      pickerInput(
        inputId = "selectize_col_factor",
        label = 'Select column of annotation',
        choices = as.vector(colnames(
          read.FCS(input$files_to_view$datapath[1])@exprs
        )),
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = F
      )
    })
    output$ui_cols_umap <- renderUI({
      pickerInput(
        inputId = "selectize_cols_umap",
        label = 'Select marker for UMAP plot',
        choices = as.vector(colnames(
          read.FCS(input$files_to_view$datapath[1])@exprs
        )),
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = F
      )
    })

  })
  
  
  
  observeEvent(input$plot_heatmap, {
    concatenated = NULL
    index_plots <-
      match(input$selectize_files, values$choices_to_add)
    #If user didn't specify any particular populations to plot
    
    for (index in index_plots) {
      a <-
        (read.FCS(input$files_to_view$datapath[index])@exprs)
      if (is.null(concatenated)) {
        concatenated = a
      }
      else{
        concatenated = rbind(concatenated, a)
      }
      
      
    }
    
    
    d <- as.data.frame(concatenated, row.names = NULL)
    #d <- lapply(d, as.numeric)
    input_is <- input$selectize_col_factor
    d[,input$selectize_col_factor] <- factor(d[,input$selectize_col_factor]) # make year a factor.
    
    if (!(is.null(input$selectize_cols_heatmap))) {
      #creating the aggregate mean per population dataframe
      df_for_heatmap <-
        aggregate(cbind(sapply(input$selectize_cols_heatmap, function(x)
          eval(as.symbol(x)))) ~ eval(as.symbol(input$selectize_col_factor)), d, mean)
      colnames(df_for_heatmap) <-
        c(input$selectize_col_factor, input$selectize_cols_heatmap)
      if (!(is.null(values$pop_annot))) {
        df_for_heatmap[, input$selectize_col_factor] <-
          values$pop_annot[, "Population"]
        if (!(is.null(input$selectize_pops))) {
          df_for_heatmap <-
            df_for_heatmap[df_for_heatmap[,input$selectize_col_factor] %in% as.vector(input$selectize_pops), ]
        }
      }
      #Rendering the heatmap
      output$plot_heatmap_map <-
        renderPlot({
          Heatmap(
            scale(df_for_heatmap[, input$selectize_cols_heatmap]),
            name = "Markers' expression",
            #title of legend
            column_title = "Markers",
            split = df_for_heatmap[,input$selectize_col_factor],
            #Split factor
            show_row_names = T,
            row_names_gp = gpar(fontsize = 1),
            row_title_rot = 0,
            row_title_side = c("right")
          )
          
        })
    }
    #Done plotting
  })

  #########################
  observeEvent(input$load_annot, {
    values$pop_annot = read.csv(input$load_annot$datapath)
    updateSelectizeInput(
      session,
      'selectize_pops',
      choices = values$pop_annot[, 'Population'],
      server = TRUE,
      options = list(maxOptions = 5)
    )
    if (!is.null(values$umap_to_plot_pops)) {
      values$pop_annot = read.csv(input$load_annot$datapath)
      legend = values$pop_annot[, 'Population']
      values$umap_to_plot_pops <- values$umap_to_plot_pops +
        scale_color_discrete(name = input$selectize_col_factor, labels = legend)
      output$plot_umap <-
        renderPlot({
          values$umap_to_plot_pops
          
        })
      
    }
    
    #Done plotting
  })
  
  observeEvent(input$load_annot, {
    values$pop_annot = read.csv(input$load_annot$datapath)
    updateCheckboxGroupInput(
      session,
      'selectize_pops',
      label = 'Select populations for heatmap plot',
      choices = values$pop_annot[, 'Population'],
      inline = T
    )
    updateCheckboxGroupInput(
      session,
      'selectize_pops',
      label = 'Select populations for heatmap plot',
      choices = values$pop_annot[, 'Population'],
      inline = T
    )
    output$plot_legend_heatmap <-
      renderPlot({
        plot(
          NULL ,
          xaxt = 'n',
          yaxt = 'n',
          bty = 'n',
          ylab = '',
          xlab = '',
          xlim = 0:1,
          ylim = 0:1
        )
        if (!is.null(values$pop_annot_heatmap)) {
          legend("topleft",
                 legend = values$pop_annot_heatmap[, 'Population'],
                 fill = myColor)
        }
      })
  })
  #Enabling heatmap and umap buttons accordingly
  observe({
    if (!is.null(input$selectize_files)) {
      if (!is.null(input$selectize_cols_heatmap)) {
        shinyjs::enable('plot_heatmap')
      }
      if (!is.null(input$selectize_cols_umap)) {
        shinyjs::enable('plot_umap')
      }
    }
  })
  observeEvent(input$refresh_input, {
    
    values$step.ff = list()
    values$step.gates = list()
    values$polygon.temp = NULL
    values$gate.strat = NULL
    values$log = matrix(
      nrow = 1,
      ncol = 1,
      data = "Upload",
      dimnames = list(NULL, "Log")
    )
    values$package.table = NULL
    showNotification("all_refresh input done", type = "message")
  })
  shinyDirChoose(input, 'dir_to_stat', roots=volumes(), session=session)
  
  
  #################

  observeEvent(input$calculate, {
    unlink(paste(
      getwd(),
      paste0('stats_calculate', '_file.xlsx'),
      sep = '/'))
    files_stat <- fixUploadedFilesNames(input$files_to_stat)
    all_files_to_pick <- sample(files_stat$datapath,10, replace = T)
    con = NULL
    file_id = 1
    first_file <- all_files_to_pick[1]
    input$selectize_stat_factor
    indexes = unique(read.FCS(first_file)@exprs[,input$selectize_stat_factor])
    #Initializing the iteration
    for (file in all_files_to_pick){
      #Extracting the matrix of expression while itterating over the files
      a<-read.FCS(file)@exprs
      #Transforming to data frame
      a<-as.data.frame(a)
      #Creating a data.frame of FILE ID and the expression matrix of the FCS files
      a<-data.frame(`File ID` = file_id, a, check.names=FALSE)
      new = NULL
      for (index in indexes){
        subsetting = subset(a, a[,input$selectize_stat_factor] == index)
        subsetting = as.matrix(subsetting)
        if (is.null(new)){
          #Sampling 20 events of each population annotation
          new = subsetting[sample(nrow(subsetting), 20, replace = T), ]
        }
        else{
          new = rbind(subsetting[sample(nrow(subsetting), 20, replace = T), ],new)
        }
      }
      if (is.null(con)){
        con = new
      }
      else{
        con = rbind(con, a)
      }
      file_id = file_id + 1
    }
    a<-read.FCS(file)
    #Dropping the last column
    a@exprs <- as.matrix(con[,-1])
    #Generating a download handler for the concatenated FCS file
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".fcs", sep="")
      },
      content = function(file) {
        write.FCS(a, file)
      }
    )
    con = con[,!sapply(con, function(x) mean((x)))==1]
    #Extracting only the data having XGBoost in their names
    data_xgboost = con[,grepl("*XGBoost",names(con))]
    file_ka = as.symbol(input$selectize_stat_factor)
    mydata = data.frame(`File ID`=con[,'File ID'], `Annotation_ID`=con[,input$selectize_stat_factor], UMAP1 = con[,'UMAP1'],UMAP2=con[,'UMAP2'],data_xgboost, check.names=FALSE)
    if ("Mean" %in% input$selectize_stats){
      mean_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, mean)
      mean_stats<-mean_stats[,-c(3,4)]
      write.xlsx(mean_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="Average", append=TRUE)
    }
    if ("Median" %in% input$selectize_stats){
      median_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, median)
      median_stats<-median_stats[,-c(3,4)]
      write.xlsx(median_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="Median", append=TRUE)
    }
    if ("Standard Deviation" %in% input$selectize_stats){
      standard_deviation_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, sd)
      standard_deviation_stats<-standard_deviation_stats[,-c(3,4)]
      write.xlsx(standard_deviation_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="SD", append=TRUE)
    }
    if ("First Decile" %in% input$selectize_stats){
      first_decile_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, first.decile)
      first_decile_stats<-first_decile_stats[,-c(3,4)]
      write.xlsx(first_decile_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="First Decile", append=TRUE)
    }
    if ("Ninth Decile" %in% input$selectize_stats){
      ninth_decile_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, ninth.decile)
      ninth_decile_stats<-ninth_decile_stats[,-c(3,4)]
      write.xlsx(ninth_decile_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="Ninth Decile", append=TRUE)
    }
    if ("First Quartile" %in% input$selectize_stats){
      first_quartile_stats <- aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, first.quartile)
      first_quartile_stats<-first_quartile_stats[,-c(3,4)]
      
      write.xlsx(first_quartile_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="First Quartile", append=TRUE)
    }
    if ("Third Quartile" %in% input$selectize_stats){
      third_quartile_stats <- aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, third.quartile)
      third_quartile_stats<-third_quartile_stats[,-c(3,4)]
      
      write.xlsx(third_quartile_stats, file = paste(
        getwd(),
        paste0('stats_calculate', '_file.xlsx'),
        sep = '/'), 
        sheetName="Third Quartile", append=TRUE)
    }
    if ("Mode" %in% input$selectize_stats){
      mode_stats<-aggregate(cbind(sapply(colnames(mydata), function(x) eval(as.symbol(x)))) ~ `Annotation_ID` + `File ID`, mydata, modefunc)
      mode_stats<-mode_stats[,-c(3,4)]
      
    }
    showNotification("Stats generation done!", type = 'message')
    output$download_fcs <- renderUI({downloadLink("downloadData", "Download concatenated FCS")})
    output$download_stats <-renderUI({downloadLink("downloadstatsData", "Download calculated stats")})
    
    output$downloadstatsData <- downloadHandler(
      filename = function() {
        paste("stats-", Sys.Date(), ".xlsx", sep="")
      },
      #Generating the text the sheets of stats that will be found inside the excel
      content = function(file) {
        if ("Mean" %in% input$selectize_stats){
          write.xlsx(mean_stats, file, 
            sheetName="Average", append=TRUE)
        }
        if ("Median" %in% input$selectize_stats){
          write.xlsx(median_stats, file, 
            sheetName="Median", append=TRUE)
        }
        if ("Standard Deviation" %in% input$selectize_stats){
          write.xlsx(standard_deviation_stats, file, 
            sheetName="SD", append=TRUE)
        }
        if ("First Decile" %in% input$selectize_stats){
          write.xlsx(first_decile_stats, file, 
            sheetName="First Decile", append=TRUE)
        }
        if ("Ninth Decile" %in% input$selectize_stats){
          write.xlsx(ninth_decile_stats, file, 
            sheetName="Ninth Decile", append=TRUE)
        }
        if ("First Quartile" %in% input$selectize_stats){
          write.xlsx(first_quartile_stats, file,
            sheetName="First Quartile", append=TRUE)
        }
        if ("Third Quartile" %in% input$selectize_stats){
        write.xlsx(third_quartile_stats, file, 
          sheetName="Third Quartile", append=TRUE)}
        if ("Mode" %in% input$selectize_stats){
        write.xlsx(mode_stats, file, 
          sheetName="Mode", append=TRUE)}
      }
    )
  }
  )
  #When loading files to generate stats on
  observeEvent(input$files_to_stat, {
    #
    output$ui_factor_stat <- renderUI({
      pickerInput(
        inputId = "selectize_stat_factor",
        label = 'Select column of annotation',
        choices = as.vector(colnames(
          read.FCS(input$files_to_stat$datapath[1])@exprs
        )),
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = F
      )
    })
    
  }
  )
}
