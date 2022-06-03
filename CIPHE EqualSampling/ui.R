#Packages used for the UI

library(shinydashboardPlus)
library(shinydashboard)
library(shinyFiles)

#Setting up the header
header <- dashboardHeader(
  title = "Ciphe EqualSampling"
)

#Setting up the sidebar
sidebar <- shinydashboardPlus::dashboardSidebar(
  minified = F,
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
    #Prevents disconnection in phases of inactivity
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
    )
  ),
  collapsed = F
)

#Setting up the body
body <- dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    
    ##### GUIDE TAB ####
    tabItem(
      tabName = "guide",
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
              "You can directly go to the 'Upload Data' tab on the left and begin uploading your files, this tool offers:"
            ),
            tags$li(
              tags$b("Upload FCS:"),
              "Equally sampling a file based on the lowest number of events observed for a certain index"
            ),
            tags$li(
              tags$b("One click process:"),
              "This minitool will take care of the steps."
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
        height = 600,
        tabPanel(
          id = 'tab_fcs_1',
          title = "Upload data",
          column(
            12,
            align = 'center',
            
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
              column(12, align = 'center', uiOutput("column_choice"))
              ,
              
            ),
            verbatimTextOutput("filechosen"),
            ),
          
          column(
            12,
            align = 'center',
            actionButton("submit", label = "Submit", width = '10%')
          ),
          #Creating the refresh button
        ),
        column(
          12,
          align = 'center',actionButton("refresh_input", "Reset"),style = "margin-top:50px")
      ),
      
    )
    
  )
)
#UI initialization
ui <- dashboardPage(header,
                    sidebar,
                    body,
                    skin = "black")
