#Packages for the UI

library(shiny)
library(shinydashboard)
library(shinyFiles)
library(rhandsontable)
library(shinyWidgets)
header <- dashboardHeader(title = "Ciphe infinityFlow",
                          uiOutput("test"),
                          dropdownMenuOutput("messageMenu"))
###THIS IS THE ONLINE VERSION test

sidebar <- dashboardSidebar(
  tags$head(tags$style(
    HTML('.content-wrapper { height: 2000px !important;}')
  )),
  ##### TAGS HEAD #####
  tags$head(
    tags$style(
      ##For shiny notification popup design
      HTML(
        ".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             "
      )
    ),
    HTML(
      "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 150000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
    )
  ),
  shinyjs::useShinyjs(),
  
  ##### MENU ITEM ####
  sidebarMenu(
    id = "tabs",
    menuItem(
      "Guide Tutorial",
      tabName = "guide",
      icon = icon("dashboard")
    ),
    menuItem(
      "Upload Data",
      icon = icon("folder-open"),
      tabName = "upload"
    ),
    menuItem(
      "Process by infinityFlow",
      icon = icon("digital-ocean"),
      tabName = "process"
    ),
    menuItem("Stats", icon = icon("calculator"), tabName = "stats"),
    menuItem("View Plots", icon = icon("eye"), tabName = "view")
  )
)
my_css <- "
.bs-select-all {
  display: none;
}
.bs-deselect-all {
  width: 100%;
}
"
body <- dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    
    ##### GUIDE TAB ####
    tabItem(
      tabName = "guide",
      shinydashboard::box(
        title = "Guide Line",
        status = "success",
        solidHeader = T,
        width = 7,
        shinydashboard::box(
          title = "Upload DATA",
          status = "primary",
          collapsed = T,
          solidHeader = T,
          width = 12,
          collapsible = TRUE,
          tags$ul(
            tags$p("Begin by uploading data."),
            tags$li(
              tags$b("Upload FCS:"),
              "Start by choosing the number of plates of your experiment."
            ),
            tags$li(
              tags$b("Upload CSV :"),
              "For each file input box, select the files that belong to the specific plate in your experiment"
            ),
            tags$li("Select the role of each acquistion channel and then the infinityMarkers."),
            tags$li("Proceed by modifiying the parameters of your experiment."),
            tags$li("Make us of the machine's multithreaded system.")
          )
        ),
        shinydashboard::box(
          title = "Preprocess",
          status = "primary",
          collapsed = T,
          solidHeader = T,
          width = 12,
          collapsible = TRUE,
          tags$ul(
            tags$li(
              "Careful that the number of CPUs shouldn't exceed the number of files."
            ),
            tags$li(
              "You can run the experiment on already transformed data or let the pipeline take care of the transformation."
            )
          )
        )
      )
    ),
    ##### UPLOAD TAB ####
    # Begins by the user setting up the number of plates to analyze then the
    # number of files input boxes would be equal to that
    tabItem(
      tabName = "upload",
      tabBox(
        id = "input",
        width = 12,
        height = 400,
        tabPanel(
          title = "Upload plates data",
          column(
            5,
            sliderInput(
              "number_plates",
              "Number of plates for the experiment:",
              min = 1,
              max = 3,
              value = 2
            )
          ),
          #####################################################################
          ###Depending on selected number of plates, the UI will be modified###
          #####################################################################
          
          conditionalPanel(condition = "input.number_plates == '1'",
                           column(
                             5,
                             fileInput(
                               "file_destination_1_1",
                               "Load .FCS/.TXT/.CSV File(s) from your desktop",
                               multiple = TRUE,
                               accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
                                 ".fcs"
                               )
                             ),
                             verbatimTextOutput("filechosen_1_1")
                           )),
          conditionalPanel(
            condition = "input.number_plates == '2'",
            column(
              5,
              fileInput(
                "file_destination_2_1",
                "Load .FCS/.TXT/.CSV File(s) from your desktop",
                multiple = TRUE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".fcs"
                )
              ),
              verbatimTextOutput("filechosen_2_1"),
              fileInput(
                "file_destination_2_2",
                "Load .FCS/.TXT/.CSV File(s) from your desktop",
                multiple = TRUE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".fcs"
                )
              ),
              verbatimTextOutput("filechosen_2_2")
            )
          ),
          
          conditionalPanel(
            condition = "input.number_plates == '3'",
            column(
              5,
              fileInput(
                "file_destination_3_1",
                "Load .FCS/.TXT/.CSV File(s) from your desktop",
                multiple = TRUE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".fcs"
                )
              ),
              verbatimTextOutput("filechosen_3_1"),
              fileInput(
                "file_destination_3_2",
                "Load .FCS/.TXT/.CSV File(s) from your desktop",
                multiple = TRUE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".fcs"
                )
              ),
              verbatimTextOutput("filechosen_3_2"),
              fileInput(
                "file_destination_3_3",
                "Load .FCS/.TXT/.CSV File(s) from your desktop",
                multiple = TRUE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".fcs"
                )
              ),
              verbatimTextOutput("filechosen_3_3")
            )
          ),
          #Creating the refresh button
          column(2, actionButton("refresh_input", "Refresh")),
          #Creating the button submit
          column(
            6,
            align = "center",
            offset = 3,
            actionButton("submit", label = "Submit", style = "margin-top:20px;")
          )
        )
      ),
      
      
      #Down TabBox for Further analysis selection of parameters
      tabBox(
        id = "antibody_selection",
        width = 12,
        height = 780,
        ##########CHECK FOR tabsetPanel
        tabPanel(
          title = "Background/Exploratory Selection",
          value = 'background_exploratory_selection',
          column(8, rHandsontableOutput('selection_table')),
          column(
            12,
            align = 'center',
            actionButton("confirm", "Confirm Background Exploratory Selection", style =
                           "margin-top:10px;")
          )
        ),
        tabPanel(
          title = "Infinity Markers selection",
          value = 'infinity_markers_selection',
          ##################################################################################################
          ####Number of tables to be generated also depends on the number of plates selected by the user####
          ##################################################################################################
          
          conditionalPanel(condition = "input.number_plates == '1'",
                           column(
                             12,
                             fileInput(
                               "infinity_markers_1_1",
                               "Choose .xl or CSV file for Infinity panel (1)",
                               accept =
                                 c('.txt',
                                   '.csv')
                             ),
                             selectInput(
                               "infinity_markers_preset_1_1",
                               "Or Choose a preset Infinity panel:",
                               choices = c(
                                 'infinity_isotypes_LEGENDSCREEN_plate_1',
                                 'infinity_isotypes_LEGENDSCREEN_plate_2',
                                 'infinity_isotypes_LEGENDSCREEN_plate_3'
                               )
                             ),
                             
                             column(8, rHandsontableOutput('infinity_markers_table_1_1'))
                           ), ),
          conditionalPanel(
            condition = "input.number_plates == '2'",
            column(
              6,
              fileInput(
                "infinity_markers_2_1",
                "Choose .xl or CSV file for Infinity panel (1)",
                accept =
                  c('.txt',
                    '.csv')
              ),
              #Nested Row
              selectInput(
                "infinity_markers_preset_2_1",
                "Or Choose a preset Infinity panel:",
                choices = c(
                  'infinity_isotypes_LEGENDSCREEN_plate_1',
                  'infinity_isotypes_LEGENDSCREEN_plate_2',
                  'infinity_isotypes_LEGENDSCREEN_plate_3'
                )
              ),
              column(6, rHandsontableOutput('infinity_markers_table_2_1'))
            ),
            column(
              6,
              fileInput(
                "infinity_markers_2_2",
                "Choose .xl or CSV files for Infinity panel (2)",
                accept =
                  c('.txt',
                    '.csv')
              ),
              selectInput(
                "infinity_markers_preset_2_2",
                "Or Choose a preset Infinity panel:",
                choices = c(
                  'infinity_isotypes_LEGENDSCREEN_plate_1',
                  'infinity_isotypes_LEGENDSCREEN_plate_2',
                  'infinity_isotypes_LEGENDSCREEN_plate_3'
                )
              ),
              column(6, rHandsontableOutput('infinity_markers_table_2_2'))
            )
            
            
            
            
          ),
          conditionalPanel(
            condition = "input.number_plates == '3'",
            column(
              4,
              fileInput(
                "infinity_markers_3_1",
                "Choose .xl or CSV file for Infinity panel (1)",
                accept =
                  c('.txt',
                    '.csv')
              ),
              #Nested Row
              selectInput(
                "infinity_markers_preset_3_1",
                "Or Choose a preset Infinity panel:",
                choices = c(
                  'infinity_isotypes_LEGENDSCREEN_plate_1',
                  'infinity_isotypes_LEGENDSCREEN_plate_2',
                  'infinity_isotypes_LEGENDSCREEN_plate_3'
                ),
                width = 400
              ),
              
              column(6, rHandsontableOutput('infinity_markers_table_3_1'))
            ),
            
            
            
            column(
              4,
              fileInput(
                "infinity_markers_3_2",
                "Choose .xl or CSV file for Infinity panel (2)",
                accept =
                  c('.txt',
                    '.csv')
              ),
              selectInput(
                "infinity_markers_preset_3_2",
                "Or Choose a preset Infinity panel:",
                choices = c(
                  'infinity_isotypes_LEGENDSCREEN_plate_1',
                  'infinity_isotypes_LEGENDSCREEN_plate_2',
                  'infinity_isotypes_LEGENDSCREEN_plate_3'
                ),
                width = 400
              ),
              
              #Nested Row
              column(6, rHandsontableOutput('infinity_markers_table_3_2'))
            ),
            
            
            column(
              4,
              fileInput(
                "infinity_markers_3_3",
                "Choose .xl or CSV file for Infinity panel (3)",
                accept =
                  c('.txt',
                    '.csv')
              ),
              selectInput(
                "infinity_markers_preset_3_3",
                "Or Choose a preset Infinity panel:",
                choices = c(
                  'infinity_isotypes_LEGENDSCREEN_plate_1',
                  'infinity_isotypes_LEGENDSCREEN_plate_2',
                  'infinity_isotypes_LEGENDSCREEN_plate_3'
                ),
                width = 400
              ),
              
              #Nested Row
              column(6, rHandsontableOutput('infinity_markers_table_3_3'))
            ),
            
            
            
            
          ),
          #To allow for ocustomized panel insertion
          column(
            12,
            align = 'center',
            actionButton("initialize_empty", "Create custom input",
                         
                         style =
                           "margin-top:10px;")
          ),
          column(
            12,
            align = 'center',
            actionButton("confirm_infinity", "Confirm Infinity markers",
                         
                         style =
                           "margin-top:10px;")
          )
          
          
        )
      )
    ),
    #Pipeline parameters tab
    tabItem(
      tabName = "process",
      tabBox(
        id = "input",
        width = 12,
        height = 400,
        tabPanel(
          title = "infinityFlow pipeline",
          column(
            5,
            sliderInput(
              "input_events_downsampling",
              "Input Events Downsampling %:",
              min = 0,
              max = 100,
              value = 50
            )
          ),
          column(
            5,
            sliderInput(
              "prediction_events_downsampling",
              "Prediction Events Downsampling %:",
              min = 0,
              max = 100,
              value = 50
            )
          ),
          column(
            6,
            align = "center",
            offset = 3,
            selectInput(
              "cores_used",
              "Cores/Threads to be used:",
              choices = c(1, 2, 4, 8, 16, 32)
            ),
            radioButtons(
              "transform",
              "Are the fcs files transformed?",
              choices = c('No', 'Yes')
            ),
            
            actionButton("begin_pipeline", label =
                           "Submit pipeline"),
            
            
            verbatimTextOutput("infinity_flow_dialog")
          ),
          
        ),
      ),
      
      
    ),
    tabItem(
      tabName = "stats",
      tabBox(
        id = "stat_tab",
        width = 12,
        height = 300,
        column(
          6,
          fileInput(
            "files_to_stat",
            "Load .FCS from your desktop",
            multiple = TRUE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".fcs"
            )
          ),
          shinyDirButton('dir_to_stat', title='Choose directory from server',label = 'Choose directory from server'),
          uiOutput('ui_factor_stat')
        ),
        column(
          6,

          checkboxGroupInput("selectize_stats", "Select stats to calculate:",
                             choices = c("Mean", "Median", "Standard Deviation", "First Decile", "Ninth Decile", "First Quartile","Third Quartile","Mode"), width = '100%', inline = T),
          actionButton("calculate", "Calculate"),
          uiOutput('download_fcs'),
          uiOutput('download_stats')
          ),
        
        
      )
      
    ),
    tabItem(
      tabName = "view",
      tabBox(
        id = "test",
        width = 12,
        height = 300,
        column(
          6,
          fileInput(
            "files_to_view",
            "Load .FCS from your desktop",
            multiple = TRUE,
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv",
              ".fcs"
            )
          ),
          uiOutput("ui_files"),
          uiOutput('ui_factor')
        ),
        column(
          6,
          fileInput(
            "load_annot",
            "Load annotation from your desktop",
            multiple = F,
            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
          ),
          checkboxGroupInput("selectize_pops", "",
                             choices = NULL, width = '100%')
        ),
        
        
      ),
      #Visualization by UMAP and heatmap possibilities inside the tool
      tabBox(
        id = "view_tab",
        width = 12,
        height = 1950,
        
        
        tabPanel(
          title = "UMAP",
          value = "view_umap",
          
          column(6,
                 uiOutput('ui_cols_umap'), ),
          
          
          column(6,
                 actionButton("plot_umap",
                              "Plot UMAP"))
          
          ,
          fluidRow(column(
            12,
            plotOutput("plot_umap_pops"),
            plotOutput("plot_umap_stats")
          ))
        ),
        
        
        tabPanel(
          title = "Heatmap",
          value = "view_heatmap",
          
          column(6,
                 
                 uiOutput('ui_cols_heatmap')),
          column(6,
                 actionButton("plot_heatmap",
                              "Plot Heatmap"))
          ,
          #Select columns for heatmap
          
          tags$head(tags$style(HTML(my_css))),
          
          
          fluidRow(column(12, plotOutput("plot_heatmap_map")), )
        )
      ),
      
      
    )
    
  )
)
#UI initialization
ui <- dashboardPage(header,
                    sidebar,
                    body)
