#Packages used in the tool

library(dplyr)
library(cipheCytoExploreR)
library(spsComps)
library(tools)
library(CytoML)
library(tidyverse)
library(xlsx)

#Increasing the limits for uploaded files
options(shiny.maxRequestSize = 100000 * 1024 ^ 2)
#Increasing the rJava memory for writing the excel file
options(java.parameters = "-Xmx8000m")


server <- function(session, input, output) {
  #Function to rename the files in the server-side shiny tool
  fixUploadedFilesNames <- function(x) {
    if (is.null(x)) {
      return()
    }
    
    oldNames = x$datapath
    newNames = file.path(dirname(x$datapath),
                         x$name)
    file.rename(from = oldNames, to = newNames)
    x$datapath <- newNames
    x
  }
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  #The package will not check for the default wd
  options("CytoExploreR_wd_check" = FALSE)
  #Global reactiveValues
  
  values <- reactiveValues(
    excel_file = NULL,
    choices_update = 'root',
    parent_child = NULL,
    markers_file = NULL,
    details_file = NULL,
    gating_template = NULL,
    channels = NULL,
    markers = NULL,
    marker1 = NULL,
    marker2 = NULL,
    marker3 = NULL,
    gs = NULL,
    backbone_specification = NULL,
    fcs_files = NULL,
    filechosen = NULL,
    experiment_name = NULL,
    experiment_markers = NULL,
    experiment_details = NULL,
    file_selected = NULL,
    merged_fs = NULL,
    flowjo_selected = NULL,
    files_flowjo_fcs = NULL,
    index_down = 0,
    experiment = NULL,
    experiment_1 = NULL,
    messageData = NULL,
    taskData = NULL,
    bar_values = 0,
    trans = NULL
  )
  ###Functions
  ###Function to render gating menu (to be called when files are loaded and processed)
  render_gating_menu <- function() {
    output$y.axis <-
      renderUI({
        selectInput("y.axis", "Y axis", choices = values$channels)
      })
    output$x.axis <-
      renderUI({
        selectInput("x.axis", "X axis", choices = values$channels)
      })
    output$invert_gating <-
      renderUI({
        radioButtons("invert_gating", "Invert gating", choices = c('No', 'Yes'))
      })
    progress$set(message = "Done successfully", value = 1)
  }
  
  update_view_options <- function() {
    #Extracting the unique name the nodes
    values$choices_update <-
      cyto_nodes(values$gs, path = 1)
    #Updating the select input of the figures
    updateSelectInput(session,
                      'parent',
                      label = 'parent',
                      choices = values$choices_update)
    #Updating the choices of populations to gate on
    updateSelectizeInput(
      session,
      'select_pops',
      choices = (gs_get_pop_paths(values$gs, path = 1)),
      server = TRUE,
      options = list(maxOptions = 8)
    )
    full_nodes_path <-
      cyto_nodes(values$gs, path = 1)
    updateSelectInput(session,
                      'download_select',
                      label = 'Download',
                      choices = full_nodes_path)
    progress$set(message = "Drawing plots...", value = 0.9)
    output$gating_scheme <- renderPlot({
      cyto_plot_gating_scheme(values$gs[1], header = NA)
    })
    output$output_view_plots_all_files <-
      renderMenu({
        menuItem("View Plots",
                 icon = icon("eye"),
                 tabName = "view")
      })
    
    updateSelectInput(
      session,
      'view_plot',
      label = 'View Plot',
      choices = input$input_fcs_flowjo$name
    )
    progress$set(message = "Done successfully", value = 1)
  }
  
  ###Disabling buttons on Loading page
  shinyjs::disable('confirm')
  shinyjs::disable('submit_flowjo')
  shinyjs::disable('submit')
  menu.df <-
    read.csv("/media/data/html/source/CIPHEGate/dropDownMenus2.csv")
  #partition menu data into message, notification, and task
  values$messageData <- subset(menu.df, type == "message")
  notificationData <- subset(menu.df, type == "notification")
  values$taskData <- subset(menu.df, type == "task")
  #Messages extracted from the csv file
  observeEvent(values$messageData, {
    output$messageMenu <- renderMenu({
      msgs <- apply(values$messageData, 1, function(row) {
        messageItem(from = row[["from"]], message = row[["message"]])
      })
      
      dropdownMenu(type = "messages",
                   badgeStatus = "success",
                   .list = msgs)
    })
  })
  
  
  output$notificationMenu <- renderMenu({
    nots <- apply(notificationData, 1, function(row) {
      notificationItem(text = row[["message"]], status = row[["status"]])
    })
    dropdownMenu(type = "notifications",
                 badgeStatus = "warning",
                 .list = nots)
  })
  
  output$taskMenu <- renderMenu({
    taks <- apply(values$taskData, 1, function(row) {
      taskItem(text = row[["message"]],
               color = row[["color"]],
               value = row[["value"]])
    })
    dropdownMenu(
      type = "tasks",
      badgeStatus = "danger",
      .list = taks,
      headerText = 'Ongoing tasks:'
    )
  })
  
  #FlowJo submission
  observeEvent(input$submit_flowjo, {
    #Values$bar_values is a reactive value to indicate the progress
    values$bar_values = values$bar_values + 30
    #Deletes any data saved before (server side)
    unlink('folder_download', recursive = TRUE)
    ##Saves on reactive values the "metadata" of input_fcs
    values$files_flowjo_fcs = dirname(input$input_fcs_flowjo$datapath[1])
    
    values$flowjo_selected = input$input_flowjo_template$datapath
    ##Display rendering of the name of the files selected by the user
    output$filechosen2 <- renderText({
      as.character(input$input_flowjo_template$name)
    })
    #Processing the FlowJo template
    ws <- open_flowjo_xml(values$flowjo_selected)
    flowjo_workspace <-
      flowjo_to_gatingset(ws, path = values$files_flowjo_fcs, name = 1)
    file.remove('Spillover-Matrix.csv')
    test =  cyto_spillover_extract(flowjo_workspace)[[1]]
    gh <- flowjo_workspace[[1]]
    values$experiment_name = input$name1
    ######################
    #Compensate
    rownames(test) <- colnames(test)
    write.csv(test, "Spillover-Matrix.csv")
    ######################
    ######################
    ######################
    ######################
    progress <- shiny::Progress$new()
    progress$set(message = "Submitting files", value = 0)
    #Setting up the variables of the filenames of details, markers
    #And gating template
    values$gating_template = paste0(values$experiment_name, "_gating-template.csv")
    cyto_gatingTemplate_generate(gh, values$gating_template)
    #Eliminating 'Comp-' from the column names of the files if available
    a <- read.csv(values$gating_template)
    a$gating_args <- gsub('Comp-', '', a$gating_args)
    a$dims <- gsub('Comp-', '', a$dims)
    file.remove(values$gating_template)
    
    
    write.csv(a, values$gating_template)
    ##FROM HERE
    values$markers_file = paste0(values$experiment_name, '_markers.csv')
    values$details_file = paste0(values$experiment_name, '_details.csv')
    #Setting up the experiments, taking only the first file no need
    #To itterate through the others too
    progress$set(message = "Setting up GatingSet", value = 0.2)
    #Fixing from 0.fcs to original names
    flowjo_new_names <-
      fixUploadedFilesNames(input$input_fcs_flowjo)
    values$gs <-
      cyto_setup(
        flowjo_new_names$datapath[1],
        gatingTemplate = 'initial.csv',
        markers = values$markers_file,
        details = values$details_file
      )
    #sampleNames(values$gs) <- values$files_selected$name
    #Reading compensate
    progress$set(message = "Setting up Compensation...", value = 0.3)
    cyto_compensate(values$gs, "Spillover-Matrix.csv")
    
    #Reading transformation
    extract_transformations <-
      gh_get_transformations(flowjo_workspace[[1]], only.function = F)
    names_channels_used <- names(extract_transformations)
    names_channels_used <- gsub('Comp-', '', names_channels_used)
    names(extract_transformations) <- names_channels_used
    progress$set(message = "Setting up the transformation...", value = 0.35)
    values$gs <-
      cyto_transform(values$gs,
                     transformerList(names_channels_used, extract_transformations))
    progress$set(message = "Applying the gating template...", value = 0.4)
    cyto_gatingTemplate_apply(values$gs,
                              gatingTemplate = values$gating_template)
    #Rendering the table of experiment markers
    values$experiment_markers <-
      read.csv(values$markers_file, sep = ",")
    progress$set(message = "Generating table for markers", value = 0.45)
    output$table_markers <- renderRHandsontable({
      rhandsontable(values$experiment_markers,
                    width = 950,
                    height = 600)
    })
    #Rendering the table of experiment details
    values$experiment_details <-
      read.csv(values$details_file, sep = ",")
    progress$set(message = "Generating table for details", value = 0.5)
    output$table_details <- renderRHandsontable({
      rhandsontable(values$experiment_details,
                    width = 950,
                    height = 600)
    })
    progress$set(message = "Rendering setup menu tab", value = 0.6)
    output$setup_menu <-
      renderMenu({
        menuItem("View Setup",
                 tabName = "setup",
                 icon = icon("wrench"))
      })
    
    
    output$download_populations <-
      renderMenu({
        menuItem("Download files",
                 tabName = "download",
                 icon = icon("download"))
      })
    
    
    progress$set(message = "Rendering cleaning menu tab", value = 0.65)
    output$cleaning_menu <-
      renderMenu({
        menuItem(
          "Cleaning/Gating",
          tabName = "clean",
          icon = icon("filter", lib = "glyphicon")
        )
      })
    values$channels <- cyto_channels(values$gs)
    values$markers <-
      as.character(cyto_markers(values$gs))
    values$markers <- append("NULL", values$markers)
    #Renders Y and X axis selection, also invert gating choice
    progress$set(message = "Rendering gating options", value = 0.7)
    #Calls for the function that renders the gating menu
    render_gating_menu()
    progress$set(message = "Updating View plots options", value = 0.75)
    #Calls for the function of rendering the viewing options
    update_view_options()
    
    
    progress$set(message = "Done processing files", value = 0.8)
    Sys.sleep(0.1)
    progress$close()
    # #############
    notif_to_add = c('',
                     'FlowJo file succesfully processed',
                     'green',
                     values$bar_values,
                     '',
                     '')
    #values$taskData <- values$taskData[-1,]
    values$taskData[2, ] <- notif_to_add
  })
  progress$set(message = "Updating options...", value = 0.9)
  
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  ##### Enables and disables gating button depending on user_entry
  ### <!> Gating names should be unique
  observe({
    if (!is.null(input$marker1)) {
      if (input$marker1  %in% values$choices_update) {
        shinyjs::disable('gating')
        shinyjs::enable('edit')
      } else{
        if (!is.null(input$marker2)) {
          if (input$marker2  %in% values$choices_update) {
            shinyjs::disable('gating')
            shinyjs::enable('edit')
            return
          } else{
            if (!is.null(input$marker3)) {
              if (input$marker3  %in% values$choices_update) {
                shinyjs::disable('gating')
                shinyjs::enable('edit')
                return
              } else{
                shinyjs::enable('gating')
                shinyjs::disable('edit')
              }
            }
            
          }
        }
        
      }
    }
    
    
  })
  
  
  
  observeEvent(input$gating, {
    #Setting up the cyto_gate_draw
    if (input$invert_gating == 'No') {
      values$marker1 = input$marker1
      values$marker2 = input$marker2
      values$marker3 = input$marker3
      if (input$marker2 == "") {
        values$marker2 = NULL
      }
      
      if (input$marker3 == "") {
        values$marker3 = NULL
      }
      
      #Function to draw the gates
      cyto_gate_draw(
        values$gs,
        parent = input$parent,
        alias = c(values$marker1, values$marker2, values$marker3),
        channels = c(input$x.axis, input$y.axis),
        type = input$type,
        negate = FALSE,
        popup = TRUE
      )
      #Retrieving the new nodes after gate plotting
      values$choices_update <- cyto_nodes(values$gs, path = 1)
      #<!> Not used right now
      # values$parent_child <-
      #   strsplit(cyto_nodes(values$gs, path = 2),
      #            split = '/',
      #            fixed = T)
      #Updating the choices
      updateSelectInput(session,
                        'parent',
                        label = 'parent',
                        choices = values$choices_update)
      full_nodes_path <- cyto_nodes(values$gs, path = 1)
      updateSelectInput(session,
                        'download_select',
                        label = 'Download',
                        choices = full_nodes_path)
    }
    #If user wants to invert gate, takes into consideration the invert gating factors
    if (input$invert_gating == 'Yes') {
      values$marker1 = input$marker1
      values$marker2 = paste0('NOT', '(', input$marker1, ')')
      cyto_gate_draw(
        values$gs,
        parent = input$parent,
        alias = c(values$marker1, values$marker2, values$marker3),
        channels = c(input$x.axis, input$y.axis),
        type = input$type,
        #Invert gating argument
        negate = TRUE
      )
      values$choices_update <- cyto_nodes(values$gs, path = 1)
      updateSelectInput(session,
                        'parent',
                        label = 'parent',
                        choices = values$choices_update)
      full_nodes_path <- cyto_nodes(values$gs, path = 1)
      updateSelectInput(session,
                        'download_select',
                        label = 'Download',
                        choices = full_nodes_path)
    }
    
    
    output$gating_scheme <- renderPlot({
      cyto_plot_gating_scheme(values$gs[1], header = NA)
    })
    
    output$output_view_plots_all_files <-
      renderMenu({
        menuItem("View Plots",
                 icon = icon("eye"),
                 tabName = "view")
      })
    
    updateSelectizeInput(
      session,
      'select_pops',
      choices = (gs_get_pop_paths(values$gs, path = 1)),
      server = TRUE,
      options = list(maxOptions = 5)
    )
  })
  
  ##Edit Button
  observeEvent(input$edit, {
    #values$gs <- cyto_compensate(values$gs)
    if (input$invert_gating == 'No') {
      values$marker1 = input$marker1
      values$marker2 = input$marker2
      values$marker3 = input$marker3
      if (input$marker2 == "") {
        values$marker2 = NULL
      }
      
      if (input$marker3 == "") {
        values$marker3 = NULL
      }
      cyto_gate_edit(
        values$gs,
        parent = input$parent,
        alias = c(values$marker1, values$marker2, values$marker3),
        channels = c(input$x.axis, input$y.axis),
        type = input$type,
        negate = FALSE
      )
      output$gating_scheme <- renderPlot({
        cyto_plot_gating_scheme(values$gs[1], header = NA)
      })
      
    }
    
    
  })
  #Enable flowJo submission buttons on data entry
  observeEvent(input$name1, {
    values$experiment_1 = input$name1
  })
  #Reactive values for the notifiations to be added to the UI
  observeEvent(input$input, {
    values$bar_values = 10
    if (req(input$input) == 'Upload data') {
      initializing_message = c('',
                               'Starting fcs gating analysis...',
                               'red',
                               values$bar_values,
                               '',
                               '')
    } else{
      #if (is.null(input$input_fcs_flowjo) &
      #    is.null(input$input_flowjo_template)){
      initializing_message = c('',
                               'Starting FlowJo template analysis...',
                               'red',
                               values$bar_values,
                               '',
                               '')
      # }
    }
    values$taskData <-
      unique(rbind(values$taskData, initializing_message))
    values$taskData[2, ] <- initializing_message
    if (nrow(values$taskData) > 2) {
      values$taskData <- values$taskData[1:2,]
    }
    
  })
  observeEvent(input$input_fcs_flowjo, {
    values$bar_values = values$bar_values + 30
    if (!is.null(input$input_fcs_flowjo)) {
      notif_to_add = c('',
                       'FCS files uploaded successfully',
                       'yellow',
                       values$bar_values,
                       '',
                       '')
      values$taskData[2, ] <- notif_to_add
    }
  })
  
  observeEvent(input$input_flowjo_template, {
    values$bar_values = values$bar_values + 30
    if (!is.null(input$input_flowjo_template)) {
      notif_to_add = c('',
                       'FlowJo template uploaded successfully',
                       'yellow',
                       values$bar_values,
                       '',
                       '')
      #values$taskData <- values$taskData[-1,]
      values$taskData[2, ] <- notif_to_add
    }
  })
  
  observeEvent(input$input_fcs, {
    values$bar_values = values$bar_values + 30
    if (!is.null(input$input_fcs)) {
      notif_to_add = c('',
                       'FCS files uploaded successfully',
                       'yellow',
                       values$bar_values,
                       '',
                       '')
      values$taskData[2, ] <- notif_to_add
    }
  })
  
  observeEvent(input$input_gating_template, {
    values$bar_values = values$bar_values + 30
    if (!is.null(input$input_gating_template)) {
      notif_to_add = c('',
                       'Template uploaded successfully',
                       'yellow',
                       values$bar_values,
                       '',
                       '')
      #values$taskData <- values$taskData[-1,]
      values$taskData[2, ] <- notif_to_add
    }
  })
  
  
  observe({
    if (!is.null(input$input_fcs_flowjo) &
        !is.null(input$input_flowjo_template) &
        (values$experiment_1 != "")) {
      shinyjs::enable('submit_flowjo')
    }
  })
  #Enable submission buttons on data entry
  observeEvent(input$name, {
    notif_to_add = c('', 'Detected...', 'yellow', '40', '', '')
    values$experiment = input$name
  })
  observe({
    if (!is.null(input$input_fcs) & (values$experiment != "")) {
      shinyjs::enable('submit')
    }
  })
  #####################################input_fcs
  ##Program behavior when FCS files are loaded
  
  
  observeEvent(input$input_fcs_flowjo, {
    ##Saves on reactive values the "metadata" of input_fcs_flowjo
    values$files_selected <- input$input_fcs_flowjo
    ##Display rendering of the name of the files selected by the user
    output$filechosen <- renderText({
      as.character(values$files_selected$name)
    })
  })
  observeEvent(input$input_fcs, {
    ##Saves on reactive values the "metadata" of input_fcs
    values$files_selected <- input$input_fcs
    ##Display rendering of the name of the files selected by the user
    output$filechosen <- renderText({
      as.character(values$files_selected$name)
    })
  })
  #Submit ordinary CytoExploreR 
  observeEvent(input$submit,
               {
                 values$bar_values = 100
                 values$experiment_name = input$name
                 progress <- shiny::Progress$new()
                 progress$set(message = "Submitting files", value = 0)
                 #Setting up the variables of the filenames of details, markers
                 #And gating template
                 values$gating_template = paste0(values$experiment_name, "_gating-template")
                 values$markers_file = paste0(values$experiment_name, '_markers.csv')
                 values$details_file = paste0(values$experiment_name, '_details.csv')
                 #Setting up the experiments, taking only the first file no need
                 #To itterate through the others too
                 progress$set(message = "Setting up GatingSet", value = 0.2)
                 test_new_names <-
                   fixUploadedFilesNames(input$input_fcs)
                 values$gs <-
                   cyto_setup(
                     dirname(test_new_names$datapath[1]),
                     gatingTemplate = values$gating_template,
                     markers = values$markers_file,
                     details = values$details_file
                   )
                 #sampleNames(values$gs) <- values$files_selected$name
                 #`cyto_names<-`(values$gs, values$files_selected$name)
                 progress$set(message = "Checking for compensation argument", value = 0.2)
                 #Check for the compensation argument
                 if (input$compensate == 'Yes') {
                   progress$set(message = "Applying compensation", value = 0.3)
                   values$gs <- cyto_compensate(values$gs)
                 }
                 #Check for the transformation argument
                 if (input$transform == 'Yes') {
                   progress$set(message = "Applying logicle transformation", value = 0.4)
                   values$trans <-
                     cyto_transformer_logicle(values$gs)
                   values$gs <- cyto_transform(values$gs,
                                               trans = values$trans)
                 }
                 #logicle transformation
                 progress$set(message = "Applying logicle transformation", value = 0.4)
                 
                 values$experiment_markers <-
                   read.csv(values$markers_file, sep = ",")
                 progress$set(message = "Generating table for markers", value = 0.45)
                 output$table_markers <- renderRHandsontable({
                   rhandsontable(values$experiment_markers,
                                 width = 950,
                                 height = 600)
                 })
                 values$experiment_details <-
                   read.csv(values$details_file, sep = ",")
                 progress$set(message = "Generating table for details", value = 0.5)
                 output$table_details <- renderRHandsontable({
                   rhandsontable(values$experiment_details,
                                 width = 950,
                                 height = 600)
                 })
                 progress$set(message = "Rendering setup menu tab", value = 0.6)
                 output$setup_menu <-
                   renderMenu({
                     menuItem("View Setup",
                              tabName = "setup",
                              icon = icon("wrench"))
                   })
                 progress$set(message = "Rendering cleaning menu tab", value = 0.65)
                 #Rendering gating tab
                 output$cleaning_menu <-
                   renderMenu({
                     menuItem(
                       "Cleaning/Gating",
                       tabName = "clean",
                       icon = icon("filter", lib = "glyphicon")
                     )
                   })
                 #Rendering Download tab
                 output$download_populations <-
                   renderMenu({
                     menuItem("Download files",
                              tabName = "download",
                              icon = icon("download"))
                   })
                 
                 values$channels <- cyto_channels(values$gs)
                 values$markers <-
                   as.character(cyto_markers(values$gs))
                 #Adding the possibility of the markers of being "NULL"
                 values$markers <- append("NULL", values$markers)
                 #Renders Y and X axis selection, also invert gating choice
                 progress$set(message = "Rendering gating options", value = 0.7)
                 render_gating_menu()
                 progress$set(message = "Updating View plots options", value = 0.75)
                 
                 updateSelectInput(
                   session,
                   'view_plot',
                   label = 'View Plot',
                   choices = values$files_selected$name
                 )
                 
                 observeEvent(input$input_gating_template, {
                   progress$set(message = "Copying and applying the provided gating template...", value = 0.8)
                   #Writing the provided gating template
                   write.csv(
                     read.csv(input$input_gating_template$datapath),
                     paste0(values$gating_template, '.csv'),
                     row.names = FALSE
                   )
                   #Function to apply it
                   cyto_gatingTemplate_apply(values$gs,
                                             gatingTemplate = input$input_gating_template$datapath)
                   values$choices_update <-
                     cyto_nodes(values$gs, path = 1)
                   #Updating the choices of parents for the gates
                   updateSelectInput(session,
                                     'parent',
                                     label = 'parent',
                                     choices = values$choices_update)
                   updateSelectizeInput(
                     session,
                     'select_pops',
                     choices = (gs_get_pop_paths(values$gs, path = 1)),
                     server = TRUE,
                     options = list(maxOptions = 5)
                   )
                   full_nodes_path <-
                     cyto_nodes(values$gs, path = 1)
                   updateSelectInput(session,
                                     'download_select',
                                     label = 'Download',
                                     choices = full_nodes_path)
                   output$gating_scheme <- renderPlot({
                     cyto_plot_gating_scheme(values$gs[1], header = NA)
                   })
                   #Render the view tab
                   output$output_view_plots_all_files <-
                     renderMenu({
                       menuItem("View Plots",
                                icon = icon("eye"),
                                tabName = "view")
                     })
                   
                 })
                 
                 progress$set(message = "Done processing files", value = 0.8)
                 Sys.sleep(0.1)
                 progress$close()
                 notif_to_add = c('',
                                  'FlowJo file succesfully processed',
                                  'green',
                                  values$bar_values,
                                  '',
                                  '')
                 #values$taskData <- values$taskData[-1,]
                 values$taskData[2, ] <- notif_to_add
               })
  progress$set(message = "Rendering gating menu", value = 0.9)
  
  observeEvent(input$invert_gating, {
    if (input$invert_gating == 'Yes') {
      output$marker1 <-
        renderUI({
          textInput("marker1", "Population", value = "")
        })
      output$marker2 <- renderUI({
        
      })
      output$marker3 <- renderUI({
        
      })
      values$marker2 <- NULL
      values$marker3 <- NULL
    }
    if (input$invert_gating == 'No') {
      output$marker1 <-
        renderUI({
          textInput("marker1", "Population 1", value = "")
        })
      output$marker2 <-
        renderUI({
          textInput("marker2", "Population 2", value = "")
        })
      output$marker3 <-
        renderUI({
          textInput("marker3", "Population 3", value = "")
        })
      
    }
  })
  ###############################################################################
  #Allows the observation and modification of experiment markers and details
  observeEvent(input$table_markers, {
    values$experiment_markers <-
      as.data.frame(hot_to_r(input$table_markers))
  })
  observeEvent(input$table_details, {
    values$experiment_details <-
      as.data.frame(hot_to_r(input$table_details))
  })
  ###############################################################################
  #Viewing plot tab
  observeEvent(input$view_plot, {
    index_fcs <- match(input$view_plot, values$files_selected$name)
    output$view_gating_scheme <- renderPlot({
      #TODO
      #Should be interactive depending on the number of gates
      cyto_plot_gating_scheme(
        values$gs[index_fcs],
        header = NA,
        contour_lines = 15,
        layout = c(5, 5)
      )
      
    })
    
  })
  
  #########################REFRESH
  observeEvent(input$refresh_input, {
    #Reset of the gatingSet
    values$gs <- NULL
    #Reset of the text render
    output$filechosen <- renderText({
      as.character('Input refreshed')
    })
    #Reset of
    updateTextInput(session, 'name', label = 'Name of the experiment', value = '')
    updateTextInput(session, 'name1', label = 'Name of the experiment', value = '')
    showNotification("all_refresh input done", type = "message")
  })
  ###Download FCS files based on gating selection
  observeEvent(input$download_button, {
    fs_gated = gs_pop_get_data(values$gs, input$download_select)
    directory_to_save <-
      paste(input$download_select,
            values$experiment_name,
            sep = '_',
            values$index_down)
    cyto_save(fs_gated, save_as = directory_to_save)
    files_downloaded <-
      list.files(path = directory_to_save, full.names = T)
    original_names_no_ext <-
      lapply(values$files_selected$name, file_path_sans_ext)
    node_with_extension <- paste0(directory_to_save, '.fcs')
    values$index_down = values$index_down + 1
    file.rename(files_downloaded,
                paste(
                  directory_to_save,
                  paste(original_names_no_ext, node_with_extension, sep = '_'),
                  sep = '/'
                ))
    showNotification(paste0("All files saved at ", directory_to_save), type = "message")
    
  })
  ##Download indexed
  observeEvent(input$download_indexed, {
    #Creating dataframe to be as a memo for the indexing process
    df_for_index <-
      data.frame(Index = character(), Population = character())
    ##Adding gating info reference to check for bugs
    percentage_for_gates <-
      data.frame(File = character())
    #Files on the left
    #Gate on TOP
    #Then percentages in between
    #Indexing variable definition
    index_event = 0
    #Initating merged_fs
    
    merged_fs = NULL
    #Stats functions
    pop.quantiles <- function(fr) {
      chnls <- colnames(fr)
      res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
      names(res) <- chnls
      res
    }
    
    #Itterating over all the selected populations
    for (pop in input$select_pops) {
      #Index
      index_event = index_event + 1
      #Creates gates of populations selected for each file
      fs_gated = gs_pop_get_data(values$gs, pop)
      #Barcoding the populations on a column of Event ID
      fs_gated = cyto_barcode(fs_gated, type = "events")
      #To keep track of the populations gated for the Log File (see below)
      df_for_index[nrow(df_for_index) + 1, ] <- c(index_event, pop)
      #Done for only the first population
      if (index_event == 1) {
        merged_fs = fs_gated
        #Initializing a first matrix concatenating all the first population
        #selection from all the samples
        
        for (i in sampleNames(fs_gated)) {
          #Replace the numbers on Event ID with '1'
          merged_fs@frames[[i]]@exprs[, 'Event ID'] <- index_event
          #  percentage_for_gates[t,pop]<- nrow(merged_fs@frames[[i]]@exprs)/50000
        }
        
      }
      else{
        for (i in sampleNames(fs_gated)) {
          #Replace the numbers on Event ID with the index of the population gated
          fs_gated@frames[[i]]@exprs[, 'Event ID'] <- index_event
          #Merging using rbind for each of the samples
          merged_fs@frames[[i]]@exprs <-
            rbind(merged_fs@frames[[i]]@exprs, fs_gated@frames[[i]]@exprs)
        }
      }
    }
    #Name of the dir to save to
    directory_to_save <-
      paste('selection',
            values$experiment_name,
            sep = '_',
            values$index_down)
    
    #Using cyto_save to save the fcs files
    cyto_save(merged_fs, save_as = directory_to_save)
    files_downloaded <-
      list.files(path = directory_to_save, full.names = T)
    #Sorting the files_downloaded vector in respect to the operating system sort
    #Exp: A10 should come after A9
    files_downloaded <- gtools::mixedsort(files_downloaded)
    #To retrieve back the original files names (0.fcs -> A1_specimen_A1.fcs)
    
    original_names_no_ext <-
      lapply(values$files_selected$name, file_path_sans_ext)
    original_names_no_ext <-
      gtools::mixedsort(unlist(original_names_no_ext))
    node_with_extension <- paste0(directory_to_save, '.fcs')
    values$index_down = values$index_down + 1
    #Renaming back to original names
    file.rename(files_downloaded,
                paste(
                  directory_to_save,
                  paste(original_names_no_ext, node_with_extension, sep = '_'),
                  sep = '/'
                ))
    #Getting the stats of gating
    quantiles_stats = gs_pop_get_stats(values$gs, type = pop.quantiles)
    percent_stats = gs_pop_get_stats(values$gs, type = 'percent')
    count_stats = gs_pop_get_stats(values$gs)
    mfi_stats = gs_pop_get_stats(values$gs, type = pop.MFI)
    #Writing the excel file with each stat being a sheet
    ##Indexing for annotation IDs and population
    write.xlsx(
      df_for_index,
      file = paste(
        directory_to_save,
        paste0(values$experiment_name, '_log_file.xlsx'),
        sep = '/'
      ),
      sheetName = "population-index",
      append = FALSE
    )
    ##Quantiles
    write.xlsx(
      quantiles_stats,
      file = paste(
        directory_to_save,
        paste0(values$experiment_name, '_log_file.xlsx'),
        sep = '/'
      ),
      sheetName = "75_quantiles_stats",
      append = TRUE
    )
    ##Percents
    write.xlsx(
      percent_stats,
      file = paste(
        directory_to_save,
        paste0(values$experiment_name, '_log_file.xlsx'),
        sep = '/'
      ),
      sheetName = "percent_stats",
      append = TRUE
    )
    ##Count
    write.xlsx(
      count_stats,
      file = paste(
        directory_to_save,
        paste0(values$experiment_name, '_log_file.xlsx'),
        sep = '/'
      ),
      sheetName = "count_stats",
      append = TRUE
    )
    ##MFI
    write.xlsx(
      mfi_stats,
      file = paste(
        directory_to_save,
        paste0(values$experiment_name, '_log_file.xlsx'),
        sep = '/'
      ),
      sheetName = "mfi_stats",
      append = TRUE
    )
    #Rendering the view menu
    output$view_populations <-
      renderMenu({
        menuItem("View Populations",
                 icon = icon("eye"),
                 tabName = "view_data")
      })
    #To render the stats inside the R shiny
    excel_read <- read.xlsx(paste(
      directory_to_save,
      paste0(values$experiment_name, '_log_file.xlsx'),
      sep = '/'
    ), 1, header = TRUE)
    #RHandsontable to render the data in shiny
    output$view_data_table <- renderRHandsontable({
      rhandsontable(excel_read,
                    width = 950,
                    height = 600)
    })
    observeEvent(input$view_data, {
      if (input$view_data == "Index") {
        excel_read <- read.xlsx(paste(
          directory_to_save,
          paste0(values$experiment_name, '_log_file.xlsx'),
          sep = '/'
        ), 1, header = TRUE)
        output$view_data_table <- renderRHandsontable({
          rhandsontable(excel_read,
                        width = 950,
                        height = 600)
        })
      }
      if (input$view_data == "Quantile .75") {
        excel_read <- read.xlsx(paste(
          directory_to_save,
          paste0(values$experiment_name, '_log_file.xlsx'),
          sep = '/'
        ), 2, header = TRUE)
        output$view_data_table <- renderRHandsontable({
          rhandsontable(excel_read,
                        width = 950,
                        height = 600)
        })
      }
      if (input$view_data == "Percent") {
        excel_read <- read.xlsx(paste(
          directory_to_save,
          paste0(values$experiment_name, '_log_file.xlsx'),
          sep = '/'
        ), 3, header = TRUE)
        output$view_data_table <- renderRHandsontable({
          rhandsontable(excel_read,
                        width = 950,
                        height = 600)
        })
      }
      if (input$view_data == "Count") {
        excel_read <- read.xlsx(paste(
          directory_to_save,
          paste0(values$experiment_name, '_log_file.xlsx'),
          sep = '/'
        ), 4, header = TRUE)
        output$view_data_table <- renderRHandsontable({
          rhandsontable(excel_read,
                        width = 950,
                        height = 600)
        })
      }
      if (input$view_data == "MFI") {
        excel_read <- read.xlsx(paste(
          directory_to_save,
          paste0(values$experiment_name, '_log_file.xlsx'),
          sep = '/'
        ), 5, header = TRUE)
        output$view_data_table <- renderRHandsontable({
          rhandsontable(excel_read,
                        width = 950,
                        height = 600)
        })
      }
    })
    #Saving done notification
    showNotification(paste0("All files saved at ", directory_to_save), type = "message")
  })
  
  #Dictating the behavior of the program on export
  observeEvent(input$export_button, {
    #Could face issues with the NOT annotation
    if (input$export_choices == 'cytoML') {
      cyto_export(values$gs,
                  paste0(values$experiment_name, "_cytoML.xml"))
    }
    # <!> For FlowJo format the docker should have an image of FlowJo for compatibility
    if (input$export_choices == 'FlowJo') {
      cyto_export(values$gs,
                  paste0(values$experiment_name, "_FlowJo.wsp"))
    }
  })
  #session$onSessionEnded(stopApp)
}
