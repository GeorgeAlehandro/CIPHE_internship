#Packages used for the UI
library(shinydashboardPlus)
library(shinydashboard)
library(shinyFiles)
library(rhandsontable)

#Setting up the header
##Notifications
header <- dashboardHeader(
  title = "Ciphe Gate",
  dropdownMenuOutput("messageMenu"),
  dropdownMenuOutput("notificationMenu"),
  dropdownMenuOutput("taskMenu")
)

#Setting up the sidebar
sidebar <- shinydashboardPlus::dashboardSidebar(
  minified = F,
  #Tags head
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
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
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
  
  #Sidebar Menu Items
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
    #Items that will show interactively
    menuItemOutput("setup_menu"),
    menuItemOutput("cleaning_menu"),
    menuItemOutput("output_view_plots_all_files"),
    menuItemOutput("download_populations"),
    menuItemOutput("view_populations")
  ),
  collapsed = F
)

#Setting up the body
body <- dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    #Library tab
    tabItem(
      tabName = "library",
      shinydashboard::box(
        title = "Library statues",
        width = 12,
        collapsible = F,
      )
    ),
    #Guide tab with interactive boxes
    tabItem(
      tabName = "guide",
      shinydashboard::box(
        title = "News in Ciphe Gating",
        status = "primary",
        solidHeader = T,
        width = 5,
        shinydashboard::box(
          width = 12,
          tags$h4("Ciphe Gating", style =
                    "font-style:'bold'"),
          tags$p("Go to upload data to begin", style =
                   "font-size:16px")
        ),
        shinydashboard::box(
          width = 12,
          tags$h4("This tool offers:"),
          tags$p("Manual gating of fcs files", style =
                   "font-size:16px;"),
          tags$p("Usage of cytoexploreR R package", style =
                   "font-size:16px;"),
          tags$p("Compatibility with FlowJo files", style =
                   "font-size:16px;")
        )
      ),
      shinydashboard::box(
        title = "Guide Line",
        
        solidHeader = T,
        width = 7,
        shinydashboard::box(
          title = "How to begin?",
          status = "primary",
          collapsed = T,
          solidHeader = T,
          width = 12,
          collapsible = TRUE,
          tags$ul(
            tags$p(
              "You can directly go to the 'Upload Data' tab on the left and based on the gating template, you will find the following choices:"
            ),
            tags$li(
              tags$b("Upload FCS:"),
              "Upload FCS files on their own without any gating template that dictates the gates"
            ),
            tags$li(
              tags$b("Bring your gating template:"),
              "Or with the data uploaded you can apply a gating template that has been generated by openCyto format."
            ),
            tags$li(
              tags$b("Using a FlowJo template:"),
              "Upload a FlowJo template and the FCS files that correspond to it, this tool offers the possibility for cross-software analysis."
            )
            
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
              "You can compensate the data in the first tab by ticking the compensation."
            ),
            tags$li(
              "Please note: If the file is in FlowJo format, it should be compensated originally by FlowJo software."
            )
          )
        )
      )
    ),
    #Upload tab
    ##Begins by the user setting up the number of plates to analyze then the
    ##number of files the input boxes would be equal to that
    tabItem(
      tabName = "upload",
      tabBox(
        id = "input",
        width = 12,
        height = 600,
        tabPanel(
          id = 'tab_fcs_1',
          title = "Upload data",
          # Accept works only on google chrome? or all browser
          # However not working on shiny browser (normal thing)
          column(
            12,
            align = 'center',
            textInput("name", "Name of the experiment", value = ""),
            fluidRow(
              column(
                12,
                align = 'center',
                fileInput(
                  "input_fcs",
                  "Load .FCS File(s) from your desktop",
                  multiple = TRUE,
                  accept = c(".fcs")
                )
              ),
              column(
                12,
                align = 'center',
                fileInput(
                  "input_gating_template",
                  "(Optional) Load gatingTemplate",
                  multiple = F,
                  accept = c(".csv")
                )
              ),
            ),
            verbatimTextOutput("filechosen"),
            fluidRow(
              column(6,
                     radioButtons("compensate", "Apply compensation?", choices = c('No', 'Yes'))),
              column(6,
                     radioButtons("transform", "Apply logicle transformation?", choices = c('No', 'Yes')))
            )),
          
          column(
            12,
            align = 'center',
            actionButton("submit", label = "Submit", width = '10%')
          ),
          #Creating the refresh button
        ),
        tabPanel(
          title = "Load FlowJo Workspace",
          
          column(
            12,
            align = 'center',
            textInput("name1", "Name of the experiment", value = ""),
            fileInput(
              "input_fcs_flowjo",
              "Load .FCS File(s) from your desktop",
              multiple = TRUE,
              accept = c(".fcs")
            )
          ),
          #flowjo_input
          column(
            12,
            align = 'center',
            fileInput("input_flowjo_template",
                      "FlowJo Template Input",
                      multiple = FALSE,
                      accept = c(".wsp")),
            verbatimTextOutput("filechosen2")
          ),
          column(
            12,
            align = 'center',
            actionButton("submit_flowjo", "Submit FlowJo Template")
          )
        ),
        column(
          12,
          align = 'center',actionButton("refresh_input", "Reset"),style = "margin-top:50px")
      ),
      
    ),
    tabItem(
      tabName = "setup",
      tabBox(
        id = "antibody_selection",
        width = 12,
        height = 680,
        tabPanel(
          title = "Channels and markers",
          value = 'channels_and_markers',
          column(8, rHandsontableOutput('table_markers')),
        ),
        tabPanel(
          title = "Files",
          value = 'files_tab',
          column(8, rHandsontableOutput('table_details')),
          
          
          
        ),
        
        
        
      ),
    ),
    
    
    tabItem(
      tabName = "clean",
      tabBox(
        id = "input",
        width = 12,
        height = 800,
        column(12,  plotOutput("gating_scheme", height = '700'))
        
      ),
      #To render all the choices and buttons of gating
      fluidRow(
        column(4, uiOutput("y.axis")),
        column(4, uiOutput("x.axis")),
        column(4, uiOutput("invert_gating"))
      ),
      fluidRow(
        column(4, uiOutput("marker1"))
        ,
        column(4, uiOutput("marker2")),
        column(4, uiOutput("marker3"))
      ),
      selectInput("parent", "parent", choices = c("root")),
      selectInput(
        "type",
        "type",
        #Different choices of gating
        choices = c(
          "rectangle",
          "polygon",
          "interval",
          "ellipse",
          "boundary",
          "quadrant"
        )
      ),
      fluidRow(column(
        6, actionButton("gating", label = "Gating")
      ),
      column(
        6, actionButton("edit", label = "Edit Gate")
      ))
    )
    ,
    
    tabItem(
      tabName = "view",
      tabBox(
        id = "input",
        width = 12,
        height = 800,
        column(12, selectInput("view_plot", "View Plot", choices = c())),
        column(12, plotOutput("view_gating_scheme", height = '1800px',width='1800px'))
      )
      
      
    ),
    tabItem(
      tabName = "view_data",
      tabBox(
        id = "input",
        width = 12,
        height = 800,
        column(12, selectInput("view_data", "View Plot", choices = c("Index", "Quantile .75", "Percent", "Count", "MFI"))),
        column(12, rHandsontableOutput("view_data_table", height = '700px'))
      )
      
      
    ),
    tabItem(
      #tabName = "download_one_pop",
      tabName = "download",
      tabBox(
        id = "files_download",
        width = 12,
        height = 400,
        tabPanel(
          title = 'Download one population',
          value = 'download_one_pop',
          column(12, selectInput(
            "download_select", "Download FCS", choices = c()
          )),
          column(12, actionButton(
            "download_button", "Download", choices = c()
          ))
        ),
        tabPanel(
          title = 'Download indexed populations',
          value = 'download_many_pops',
          column(
            12,
            selectizeInput(
              'select_pops',
              'Select populations to index fcs files',
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(placeholder = 'Select populations to index')
            )
          ),
          column(12, actionButton("download_indexed", "Download"))
        )
      ),
      tabBox(
        id = "export_as",
        width = 12,
        height = 400,
        tabPanel(
          title = 'Export Gating Set',
          value = 'export_gs',
          column(12, selectInput(
            "export_choices", "Download FCS", choices = c('cytoML','FlowJo')
          )),
          column(12, actionButton(
            "export_button", "Export", choices = c()
          ))
        )
      )
      
      
    )
    
  )
)
#UI final setup
ui <- dashboardPage(header,
                    sidebar,
                    body,
                    skin = "blue")
