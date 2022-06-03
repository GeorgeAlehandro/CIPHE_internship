library(shiny)
library(shinydashboard)
library(shinyFiles)
library(rhandsontable)
header <- dashboardHeader(
  title="Ciphe Gating",
  uiOutput("test"),
  dropdownMenuOutput("messageMenu")
)

sidebar <- dashboardSidebar(

  ##### TAGS HEAD #####
  tags$head(
    tags$style(
      ##For shiny notification popup design
      HTML(".shiny-notification {
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
  
  ##### MENU ITEM ####
  sidebarMenu(id = "tabs",
              menuItem("Guide Tutorial", tabName = "guide", icon = icon("dashboard")),
              menuItem("Upload Data", icon = icon("folder-open"), tabName = "upload"),
              menuItem("Cleaning",tabName="clean", icon = icon("filter",lib="glyphicon")),
              menuItem("View Plots",icon=icon("eye"),tabName = "view"),
              menuItem("Clustering",icon=icon("braille"),
                       menuSubItem("Run Clustering",icon=icon("braille"),tabName="clusters"),
                       menuSubItem("Clusters Analysis",icon=icon("braille"),tabName="clustersAnalysis")
              ),
              menuItem("BackGating",icon=icon("undo-alt"),tabName = "backgating"),
              # menuItem("OneSENSE",icon=icon("bolt"),tabName="onesens"),
              menuItem("Annotation",icon=icon("magic"),
                       menuSubItem("Manual",icon=icon("hand-spock"),tabName="manuannot"),
                       menuSubItem("Automatique",icon=icon("robot"),tabName="autoannot")
              ),
                  column(12,actionButton("bigDDl","Download FCS")),tags$br(),tags$br(),
                  column(12,uiOutput("lnkBigDDl")),tags$br(),tags$br()

                ,
                column(12,verbatimTextOutput('ex_out'))
              
  )
)

body <- dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    ##### LIBRARY #####
    tabItem(tabName="library",
            shinydashboard::box(title = "Library statues", width=12,collapsible = F,
            )
    ),
    ##### GUIDE TAB ####
    tabItem(tabName = "guide",
            shinydashboard::box(title="News in Ciphe infinityFlow", status="success", solidHeader = T,width=5,
                                shinydashboard::box(width=12,
                                                    tags$h4("infinityFlow Pipeline",style="font-style:'bold'"),
                                                    tags$p("Go to upload data to begin",style="font-size:16px")
                                ),
                                shinydashboard::box(width=12,
                                                    tags$h4("To do list"),
                                                    tags$p("Create one tab with label edit and keywords edit", style="font-size:16px;"),
                                                    tags$p("Change plot view of backgating with hypergate", style="font-size:16px;"),
                                                    tags$p("Add select annotation in View FCS tab",style="font-size:16px;")
                                )
            ),
            shinydashboard::box(title="Guide Line",status="success",solidHeader = T,width=7,
                                shinydashboard::box(title="Upload DATA", status = "primary",collapsed = T,  solidHeader = T, width=12,collapsible=TRUE,
                                                    tags$ul(
                                                      tags$p("Ce premier onglet vous permet de charger vos donn??es,il est obligatoire de passer par ce premiere onglet
            et de charg?? des fichier pour d??bloqu?? la possibilit?? d'utiliser les autres onglet et algorithme disponible dans
            la CIPHE-Box. VOus avez la possibilit?? d'entrer des donn??es diff??rentes maniere et de diff??rents formats"),
                                                      tags$li(tags$b("Upload FCS:"),"Vous pouvez uploader un ou plusieurs fichier FCS (3.0) il seront alors lue dans l'ordre d'uploade dans une liste"),
                                                      tags$li(tags$b("Upload CSV :"), "A la palce de FCS vous pouvez loader plusieurs CSV qui serons alors convertit en FCS (Attention ne m??langer pas FCS et CSV)"),
                                                      tags$li("FLowRepository : Vous pouvez charg?? des fichier FCS depuis un repertoire pr??sent sur FlowRepository. Cela est un peu plus long car il y a le temps de t??l??chargement a rpendre en copmpte et une connection internet est requise."),
                                                      tags$li("From Server : C'est en cours de developpement mais il seras possible d'entrer une addresse IP et un chemin pour permettre a l'outil de scanner les sous dossier pr??sent et proposer le chargement des FCS dans l'outil"),
                                                      tags$li("RData : En reflection, proposer un format de sauvegarde du travail en cours au format RData et en relecture en upload.")
                                                    )
                                ),
                                shinydashboard::box(title="Preprocess", status = "primary",collapsed = T,   solidHeader = T, width=12,collapsible=TRUE,
                                                    tags$ul(
                                                      tags$li("You can Compensate, Transform, Divide , Multiple, Center and download result under menu"),
                                                      tags$li("Its possible to explore your parameters witxh overview of each label for one file or each file for one label")
                                                    )
                                )
            )
    ),
    ##### UPLOAD TAB ####
    # Begins by the user setting up the number of plates to analyze then the
    # number of files input boxes would be equal to that
    tabItem(tabName = "upload",
            tabBox(id = "input",width=12,height = 400,
                   tabPanel(title="Upload data",
              # Accept works only on google chrome? or all browser
              # However not working on shiny browser (normal thing) 

                              column(12, align = 'center',
                                     fileInput("input_fcs", "Load .FCS/.TXT/.CSV File(s) from your desktop", multiple = TRUE,
                                               accept = c("text/csv","text/comma-separated-values,text/plain",".csv",".fcs")
                                     ),
                                     verbatimTextOutput("filechosen")),                
              #Creating the refresh button
                            column(2,actionButton("refresh_input","Refresh")),
              #Creating the button submit
              column(6, align="center", offset = 3,actionButton("submit", label="Submit",style="margin-top:20px;"))
                   ),
                   tabPanel(title="Load RData",
                            column(2,fileInput("values_input","RData Input",multiple=FALSE, accept=c(".RData")))
                   )
            ),
            
          
    ),
    tabItem( tabName="clean",
             column(9,
                    uiOutput("selectViewPlot"),
                    plotOutput("cleaningPlot")
             ),
             column(3,
                    div(class="live",fluidRow(
                      column(6,uiOutput("clean.live1")),
                      column(6,uiOutput("clean.live2"))
                    )),
                    div(class = "live_number", fluidRow(
                      column(6,sliderInput("clean.gate.live1",
                                           "Live Gate X",
                                           min = 0, max = 250000, value = c(0,5000), step = 0.1
                      )),
                      column(6,sliderInput("clean.gate.live2",
                                           "Live Gate Y",
                                           min = 0, max = 250000, value = c(0,5000), step = 0.1
                      )) 
                    )),
                    div(class="size",fluidRow(
                      column(6,uiOutput("clean.size1")),
                      column(6,uiOutput("clean.size2"))
                    )),
                    div(class = "size_number", fluidRow(
                      column(6,sliderInput(
                        "clean.gate.size1",
                        "Size Gate X",
                        min = 0, max = 250000, value = c(45000, 55000), step = 5000
                      )),
                      column(6,sliderInput(
                        "clean.gate.size2",
                        "Size Gate Y",
                        min = 0, max = 250000, value = c(45000, 55000), step = 5000
                      ))
                    )),
                    div(class= "sSSC",fluidRow(
                      column(6,uiOutput("clean.singletsFCS1")),
                      column(6,uiOutput("clean.singletsFSC2"))
                    )),
                    div(class="sFSC",fluidRow(
                      column(6,
                             sliderInput("x.border.left",
                                         "Singlet FSC: bot left X",
                                         min = 0, max = 250000, value = 0, step = 5000)
                      ),
                      column(6,
                             sliderInput("y.border.left",
                                         "Singlet FSC: bot left Y",
                                         min = 0, max = 250000, value = 0, step = 5000)
                      ),
                      column(6,
                             sliderInput("x.border.right",
                                         "Singlet FSC: top right X",
                                         min = 0, max = 250000, value = 100000, step = 5000)
                      ),
                      column(6,
                             sliderInput("y.border.right",
                                         "Singlet FSC: top right Y",
                                         min = 0, max = 250000, value = 250000, step = 5000)
                      )
                    )),
                    div(class="sSSC",fluidRow(
                      column(6,uiOutput("clean.singletsSSC1")),
                      column(6,uiOutput("clean.singletsSSC2"))
                    )),
                    actionButton("gating", label = "Gating")
             )            
    ),
    
    tabItem(tabName = "view",
            tabBox(id = "input",width=12,height = 2000,
                     
                   
                   tabPanel(title="UMAP Annotated",
                            
                            column(5,actionButton("generate_umap", "Generate PDF"),uiOutput("plot_umap"))
                   ),
                   tabPanel(title="UMAP Annotated Background Corrected",
                            
                            column(5,actionButton("generate_umap_corrected", "Generate PDF"),uiOutput("plot_umap_corrected"))
                   )
            ),
            
            
    )
    
  )
)
##### UI ######
ui <- dashboardPage(
  header,
  sidebar,
  body
)










