require("shiny")
require("shinydashboard")
require("rhandsontable")
require("plotly")
require("DT")
require("shinyjs")

header <- dashboardHeader(
  title="Ciphe-Box",
  uiOutput("test"),
  dropdownMenuOutput("messageMenu")
)

sidebar <- dashboardSidebar(
  ##### TAGS HEAD #####
  tags$head(
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
    menuItem("Preprocess",icon=icon("digital-ocean"),tabName="preprocess"),
    menuItem("View FCS",icon=icon("eye"),tabName = "view"),
    menuItem("Manage FCS", icon=icon("edit"),
      menuSubItem("Configuration",icon=icon("code"),tabName="conf"),
      menuSubItem("Concatenate", icon=icon("boxes"),tabName="concat"),
      menuSubItem("Descriptive value",icon=icon("sort-numeric-up"),tabName="mfi"),
      menuSubItem("Diferential value",icon=icon("sort-numeric-up"),tabName="difstats")
    ),
    menuItem("Clean FCS",icon=icon("broom"),
      menuSubItem("FlowAI", icon=icon("digital-ocean"),tabName="flowai"),
      menuSubItem("FlowClean",icon=icon("digital-ocean"),tabName="flowclean")
    ),
    menuItem("Sampling", icon=icon("angle-double-down"),
      menuSubItem("Random", icon=icon("compress"),tabName="rand"),
      menuSubItem("Percentile ",icon=icon("fingerprint"),tabName="percent")
    ),
    menuItem("Reduction Dimension", icon = icon("cuttlefish"),
      menuSubItem("PCA (linear)",icon=icon("cut"),tabName="pca"),
      menuSubItem("Unlinear",icon=icon("wrench"),tabName="umap"),
      menuSubItem("oneSENSE",icon=icon("bolt"),tabName="onesens")
    ),
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
    menuItem("Manual Gating",icon=icon("hand-spock"),
      menuSubItem("Gate",icon=icon("tree"),tabName="gate"),
      menuSubItem("View",icon=icon("road"),tabName="gateview")
    ),
    conditionalPanel(
      condition = "output.fileUploaded",
      # column(6,actionButton("save","Save RData")),
      # column(6,uiOutput("lnkRdata")),tags$br(),tags$br(),
      conditionalPanel(
        condition = "output.bigFileUploaded",
        column(12,actionButton("bigDDl","Download FCS")),tags$br(),tags$br(),
        column(12,uiOutput("lnkBigDDl")),tags$br(),tags$br()
      ),
      conditionalPanel(
        condition = "output.smallFileUploaded",
        column(12,downloadButton("downloadData","Download FCS")),tags$br(),tags$br()
      ),
      column(12,verbatimTextOutput('ex_out'))
    )
  )
)

body <- dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    ##### LIBRARY #####
    tabItem(tabName="library",
      box(title = "Library statues", width=12,collapsible = F,
        column(6,dataTableOutput("lib.bio.table",height = 800)),
        column(3,uiOutput("selectUpdate")),
        column(3,actionButton("updatePackage","Update/Install"))
      )
    ),
    ##### GUIDE TAB ####
    tabItem(tabName = "guide",
      shinydashboard::box(title="News in CIPHE-Box", status="success", solidHeader = T,width=5,
        shinydashboard::box(width=12,
          tags$h4("OneSENSE Pipeline",style="font-style:'bold'"),
          tags$p("You can use OneSense but the result its just one files with 4 new dimensions",style="font-size:16px")
        ),
        shinydashboard::box(width=12,
          tags$h4("Add down Sampling by Density",style="font-style:'bold';"),
          tags$p("We add a down sampling method with a kernel density estimation",style="font-size:16px;")
        ),
        shinydashboard::box(width=12,
          tags$h4("Change tab organisation"),
          tags$p("Change organisation of first tab, upload Data is for upload data, check number of events, dimensions,
                 and edit name of file (in progress)",style="font-size:16px;"),
          tags$p("Preprocess part is for transform, compensate, scale and divide data with two overview,
                 one with value of 10 first events, and density plot (works just if you have the same number and names of parameters"
                 ,style="font-size:16px;")
          
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
        ),
        shinydashboard::box(title="View FCS", status = "primary",collapsed = T,   solidHeader = T, width=12,collapsible=TRUE,
          tags$ul(
            tags$li("You can Compensate, Transform, Divide , Multiple, Center and download result under menu"),
            tags$li("Its possible to explore your parameters witxh overview of each label for one file or each file for one label")
          )
        ),
        shinydashboard::box(title="Manage FCS", solidHeader = T,collapsed = T,  width=12, collapsible=T, status="primary",
          shinydashboard::box(title="Configuration FCS", status = "primary",  solidHeader = F,collapsed = T, width=12, collapsible=TRUE,
            tags$p("Update labels and names in FCS"),
            tags$p("Update keywords and download with right-click")
          ),
          shinydashboard::box(title="Compute MFI", status = "primary",  solidHeader = F, width=12,collapsed = T, collapsible=TRUE,
            tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
              labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
              laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
              voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
              cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
            )
          )
        ),
        shinydashboard::box(title="Clean FCS", status = "primary",collapsed = T,   solidHeader = T, width=12,collapsible=TRUE,
          shinydashboard::box(title="FlowAI",status="primary", collapsed=T,width=12,collapsible = T,
            tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
              labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
              laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
              voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
              cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
            )
          ),
          shinydashboard::box(title="FlowClean",status="primary", collapsed=T,width=12,collapsible = T,
            tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
              labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
              laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
              voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
              cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
            )
          ),
          shinydashboard::box(title="FlowCut",status="primary", collapsed=T,width=12,collapsible = T,
            tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
              labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
              laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
              voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
              cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
            )
          )
        ),
        shinydashboard::box(title="Sampling", status = "primary",collapsed = T,   solidHeader = T, width=12,collapsible=TRUE,
          tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
            labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
            laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
            voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
            cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
          )
        ),
        shinydashboard::box(title="Dimension Reduction", status = "primary",collapsed = T,   solidHeader = T, width=12,collapsible=TRUE,
          tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
            labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
            laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
            voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
            cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
          )
        ),
        shinydashboard::box(title="Clustering", status = "primary",collapsed = T, solidHeader = T, width=12, collapsible=TRUE,
          tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
            labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
            laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
            voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
            cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
          )
        ),
        shinydashboard::box(title="BackGating", status = "primary",collapsed = T, solidHeader = T, width=12, collapsible=TRUE,
          tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut 
            labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco 
            laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in 
            voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat 
            cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
          )
        )
      )
    ),
    ##### UPLOAD TAB ####
    tabItem(tabName = "upload",
      tabBox(id = "input",width=12,height = 150,
        tabPanel(title="Upload FCS",
          column(4,
            fileInput("fcs_input","Files Input (FCS or CSV to convert)",
            multiple=TRUE, accept=c(".fcs",".csv"))
          ),
          column(2,actionButton("refresh_input","Refresh")),
          column(2,selectInput("sep","Separateur", choice=c(";",",","tab")))
        ),
        tabPanel(title="Flow Repository",
          column(3,textInput("searchFR","Research")),
          column(3,uiOutput("listFR")),
          column(1,uiOutput("uiNbrFilesFR")),
          column(1,actionButton("searchFR","Search")),
          column(1,actionButton("loadFR","Load FCS"))
        ),
        tabPanel(title="From Server",
          column(2,textInput("IP","IP","10.71.1.22")),
          column(2,textInput("path","Path","/INPUT/DATA")),
          column(3,uiOutput("uiListDirs"))
        ),
        tabPanel(title="Load RData",
          column(2,fileInput("values_input","RData Input",multiple=FALSE, accept=c(".RData")))
        )
      ),
      box(width=12,
        conditionalPanel(
          condition = "output.fileUploaded",
          column(2,actionButton("renameFCS","Rename FCS")),
          column(2,actionButton("fcsToCSV","Convert to CSV")),
          column(2,uiOutput("lnkBigDDlCSV")),
          tags$br(),tags$br(),
          column(12,dataTableOutput("myTable"))
        )
      )
    ),
    ##### PREPROCESS TAB ####
    tabItem(tabName="preprocess",
      conditionalPanel(condition = "output.fileUploaded",
        box(width=3, status="warning", solidHeader=TRUE,title="Preprocess", collapsible=F,
          box(width = 12,
            actionButton("comp","Compensate"),
            actionButton("decomp","Decompensate")
          ),
          box(width = 12,
            uiOutput("transMarkerOutput")
          ),
          box(width=12, collapsible=TRUE, title="Transformation",
            selectInput("trans_meth","Transformation",choices=c("logicle","arcsinh","flowVS"), multiple=FALSE),
            numericInput("args_trans","Args",value=NULL),
            actionButton("addAllTransParams","All"),
            actionButton("transform","Transform"),
            actionButton("detransform","Detransform")
          ),
          box(width = 12, collapsible = TRUE, title="Change value",collapsed=TRUE,
            column(4,actionButton("divide","Divide")),
            column(5,actionButton("multiple","Multiple")),
            column(3,numericInput("mult_value","",value=1))
          ),
          box(width = 12,title="Normalize", collapsible = TRUE, collapsed= TRUE,
            column(6,actionButton("center","Center")),
            column(6,actionButton("scale","Scale")),
            column(6,numericInput("quantConst","Probs",value=97.5,min = 10,max=99)),
            column(6,actionButton("norm","Constan Norm")),
            column(12,actionButton("clr_fn","Centred Log-Ratio"))
          )
        ),
        box(title="Values OverView", status = "warning",  solidHeader = TRUE,width=9,collapsible=TRUE,collapsed = TRUE,
          dataTableOutput("textPreview"),style = "width:100%;overflow-x: scroll;font-size:12px;"
        ),
        box(title="Density OverView",status = "warning", solidHeader = TRUE, width = 9, collapsible = TRUE,
          column(2,selectInput("flowVizMod","Mode",choices = c("FlowViz","ggridge"),multiple=F,selected=1)),
          column(2,uiOutput("selectFlowViz")),
          column(2,numericInput("jpxmin","Xmin",value=0)),
          column(2,numericInput("jpxmax","Xmax",value=100)),
          column(2,numericInput("height","Height",value=600,min = 60,step=60)),
          column(2,numericInput("width", "Width", value=1200)),
          # column(1,downloadButton("ddljp","Download PNG")),
          plotOutput("overView",height="auto")
        )
      )
    ),
    ##### VIEW TABS ####
    tabItem(tabName = "view",
      conditionalPanel(condition="output.fileUploaded",
        tabBox(id="viewTabBox",width=12,
          tabPanel("viewSimpleFile",
            fluidRow(
              box(width=4,
                column(6,selectInput("view_mod","Select View Mod",choices=c("Scatterplot","Heatmap"), multiple=FALSE)),
                column(6,sliderInput("samplingView","Freq Cells",min=0.1,max=1,value=0.3))
              ),
              box(width=8,
                conditionalPanel(condition="input.view_mod == 'Heatmap'",
                  column(12,uiOutput("selectHeatmapDim"))
                ),
                conditionalPanel(condition="input.view_mod == 'Scatterplot'",
                  column(6,uiOutput("selectXPreview")),
                  column(6,uiOutput("selectYPreview"))
                )
              )
            ),
            fluidRow(
              box(width=3,
                column(12,uiOutput("selectFilePreview")),
                column(6,uiOutput("uiAnnotParams")),
                column(6,uiOutput("uiAnnotID")),
                conditionalPanel(condition="input.view_mod == 'Scatterplot'",
                  column(6,uiOutput("xlimOutput")),
                  column(6,uiOutput("ylimOutput")),
                  column(6,sliderInput("cex","Cex",min=0.1,max=10,value=1))
                ),
                column(12,uiOutput("colorMods"))
              ),
              box(width=9,
                fluidRow(
                  plotOutput("plotPreview",height="auto")
                )
              )
            )          
          ) 
        )
      )
    ),
    ##### CONFIGURATION TAB ####
    tabItem(tabName = "conf",
      conditionalPanel(
        condition = "output.fileUploaded",
        tabBox(id = "tabConf",width = 12,title = "",
          tabPanel("Edit Labels and Names",
            fluidRow(
              column(2,selectInput("selectViewLabels","Select Parameters (Names/Labels)",choices=c("names","labels"))),
              column(2,actionButton("saveParams","Save Parameters")),
              column(2,textInput("markerPattern","Pattern")),
              column(2,radioButtons("wherePattern","Where",choices=c("post","pre"),selected="post",inline = TRUE)),
              column(2,actionButton("addPattern","Add Pattern"))
            ),
            rHandsontableOutput('paramsTable')
          ),
          tabPanel("Edit Keywords",
            actionButton("applyKeywords","Apply Keywords"),
            rHandsontableOutput("keywordsTable")
          ),
          tabPanel("Edit Expression Matrice",
            fluidRow(
              box("Delete Row", width = 8,
                uiOutput("uiSelectEditAnnot"),
                actionButton("deleteRow","Delete Rows"),
                rHandsontableOutput("deleteRowExprs")
              ),
              box(width=4, collapsible = T, collapsed = T,title="How to use",
                tags$ul(
                  tags$li("Selectionner un parametre de votre table d'expression"),
                  tags$li("Si et seulement si il y a moin de 100 valeur unique diff??rentes dans cette dimension la table va apparaitre.
                    Cette table affiche le nombre d'evenement correspondant ac ette valeur dans l'ensemble de vos fichiers"),
                  tags$li("Cocher les valeur pour lesquels vous souhaitez retirer les evenements sur l'ensemble des fichier "),
                  tags$li("Cliqu?? sur le bouton Delete Rows et cela va supprimer de vos FCS les evenements")
                )
              )
            )
          )
        )
      )
    ),
    ##### CONCA TAB ####
    tabItem(tabName="concat",
      conditionalPanel(
        condition = "output.fileUploaded",
        tabBox(title = "", id="tab1",width=12, height="900px",
          tabPanel(title="Concatenante",
            shinydashboard::box(width=6,
                box(width=12,collapsible = T,title="How To use",collapsed = T,
                  tags$ul(
                    tags$li("Pour concaterner des fichier, selectionner le ou les fichier que vous souhaitez rassembler en un
                      puis sur le bouton add clustering group."),
                    tags$li("A droite va apparaitre alors un seul logo 'trash' qui sert a annuler
                      votre selection et vous indique bien que un seul fichier seras cr??er a partir de ceux choisi.") ,
                    tags$li("Vous pouvez cr??er plusieurs concatenation d'un seul coup et meme utiliser plusieurs fois le meme fichiers"),
                    tags$li("Une fois votre selection faites choisissez le nom du parametre permettant de retrouver quelle evenement correspond a quelle fichier (par ordre numerique)
                         et cliqu?? sur Concat Files"),
                    tags$li("Pour s'assurer que votre concatenation a bien march??, la liste de gauche seras remplacer par le ou les fichiers concatener et vous pouvez utiliser ces meme ficheir dans d'autre analyses")
                )),
                uiOutput("clusteringui2"),
                column(6,actionButton("clusteringui_add_clustering_group", "Add clustering group")),
                column(6,actionButton("clusteringui_add_all_groups", "Add all into groups"))
            ),
            shinydashboard::box(width=6,
              uiOutput("clusteringui3"),
              textInput("concatParams","Concatenante Params",value="Flag"),
              actionButton("addConcat","Concat Files")
            )
          ),
          tabPanel(title="Split FCS",
            box(width=8,
              uiOutput("selectSepMarkers"),
              tableOutput("theoricNbrEvents"),
              actionButton("selectSep","Select Sep params")
            ),
            box(width=4, collapsible = T, collapsed = T, title="How to use")
          )
        )
      )
    ),
    ##### FLOWAI TAB ####
    tabItem(tabName = "flowai",
      conditionalPanel(
        condition="output.fileUploaded",
        column(3,
          shinydashboard::box(title="Clean with FlowAI",collapsible = FALSE, solidHeader = TRUE,width=12,
            uiOutput("selectFlowAIStep"),
            uiOutput("selectChFM"),
            uiOutput("selectChFS"),
            actionButton("clean","Apply Clean")
          )
        ),
        conditionalPanel(
          condition="output.fileCleanAI",
          column(9,
            tabBox(title="FlowAI Results",width=12,id="tabBox2",height="900px",
              tabPanel("Overview",
                column(3,
                  uiOutput("IDFlowAIPreviewxOutput"),
                  column(6,uiOutput("selectXFlowAI")),
                  column(6,uiOutput("selectYFlowAI")),
                  textInput("flagQC","Flag QC param",value="QC"),
                  actionButton("addFlowAI","Add Quality Params"),
                  actionButton("selectHQ","Select High QC")
                ),
                column(7,plotOutput("plotPreviewFlowAI",height=500,width=500)),
                shinydashboard::box(width=12,uiOutput("tableQC_ui"))
              ),
              tabPanel("FlowRate",uiOutput("selectFlowRate"),plotOutput("plotFlowRate")),
              tabPanel("Signal Acquisition",uiOutput("selectSignal"),plotOutput("plotSignal",height=700)),
              tabPanel("Dynamic Range",uiOutput("selectDynamic"),plotOutput("plotDynamic"))
            )
          )
        )
      )
    ),
    ##### FLOWCLEAN TAB ####
    tabItem(tabName = "flowclean",
      conditionalPanel(condition="output.fileUploaded",
        box(width=12,title="Clean with FlowClean package",
          fluidRow(column(6,selectInput("vectMarkers","Select Params",choices=NULL, multiple=TRUE)),
          column(2,textInput("grepVectMarkers","Pattern")),
          column(1,actionButton("addGrepVectMarkers","Add Pattern")),
          column(1,actionButton("allVectMarkers","Add All")),
          column(1,actionButton("clearVectMarkers","Clear All"))),
          box(width=3,
            actionButton("runFlowClean","Run"),
            sliderInput("binSize","binSize",value = 0.01,min=0,max=1,step=0.01),
            sliderInput("nCellCutoff","nCellCutoff",value=500,min=100,max=1000,step=100),
            uiOutput("uiSelectXflowClean"),
            uiOutput("uiSelectYflowClean"),
            uiOutput("uiSelectIDflowClean"),
            actionButton("enrichFlowClean","Enrich FCS")
          ),
          box(width=9,
            column(6,plotOutput("plotFlowClean",width = 550, height = 550)),
            column(6,tableOutput("tablResultFlowCLEAN"))
          )
        )
      )
    ),
    ##### BACKGATING ####
    tabItem(tabName = "backgating",
      conditionalPanel(condition="output.fileUploaded",
      box(title="Preview your data", solidHeader = TRUE, status = "success", collapsible = F, width = 4,
        column(10,uiOutput("selectFilesBackGating")),
        fluidRow(column(6,uiOutput("UIannotation")),
                 column(6,selectInput("annotation_id","Select Population", choices = c("")))
        ),
        fluidRow(column(6, uiOutput("UImarkerX")),
                 column(6, uiOutput("UImarkerY"))
        ),
        column(10, plotOutput("plotPreviewBG", height="400px"))
      ),
      box(title="HyperGate Settings", solidHeader = TRUE, status = "warning", collapsible = F, width = 8,
        shinydashboard::box(width=12, collapsible = T,
          column(8,uiOutput("selectMarkersBG")),
          column(4,selectInput("methBG","Algorithme",choice=c("HyperGate","GateFinder"))),
          column(1,actionButton("clearAllBGDim","Clear List")),
          column(2,textInput("grepBG","Pattern")),
          column(1,actionButton("addGrepBG","Add Pattern")),
          column(4,actionButton("addAllMarkersBG","Add All")),
          fluidRow(
            column(4,
              conditionalPanel(condition="input.methBG=='HyperGate'",
                actionButton("runHypergate","Run Hypergate")
              ),
              conditionalPanel(condition="input.methBG=='GateFinder'",
                actionButton("runGateFinder","Run GateFinder")
              )
            )
          )
        ),
        conditionalPanel(condition = "input.methBG == 'HyperGate'",
          shinydashboard::box(width=12, collapsible = T,
            column(12,tableOutput("tableOutput"))),
          shinydashboard::box(width=12, collapsible = T,
            column(12,uiOutput("gatingStrat"))),
          shinydashboard::box(width=12, collapsible = T,
            column(12,plotlyOutput("barPlot")))
        ),
        conditionalPanel(condition = "input.methBG == 'GateFinder'",
          column(12,plotOutput("gatingStrat2")),
          column(12,plotOutput("linePlot"))
        )
      )
    )),
    ##### RANDOM TAB ####
    tabItem(tabName="rand",
      conditionalPanel(
        condition = "output.fileUploaded",
        shinydashboard::box(title="Random Downsamling" ,collapsible= TRUE, width = 8,status="danger", solidHeader=TRUE,
          shinydashboard::box(title="For all files",width=12,collapsible = T,
            column(3,radioButtons("sampleMode","Sampling",c("percent","events"),inline=TRUE)),
            column(6,
              conditionalPanel(
                condition="input.sampleMode=='percent'",
                sliderInput("percentile","%",1,100,value=100,sep=5)
              ),
              conditionalPanel(
                condition = "input.sampleMode=='events'",
                numericInput("events","Nbr events",value=1000)
              )
            ),
            column(3,
              actionButton("applyAllSampling","Apply All Sampling")
            )
          ),
          shinydashboard::box(title="Individual Sampling",width=12,collapsible = FALSE,
            actionButton("applyIndSampling","Apply Individual Sampling"),
            conditionalPanel(
              condition="input.sampleMode=='percent'",
              uiOutput("sliderInputByFiles")
            ),
            conditionalPanel(
              condition = "input.sampleMode=='events'",
              uiOutput("numberInputByFiles")
            )
          )
        ),
        shinydashboard::box(title="Files and size", collapsible=TRUE, width=4,status="danger",solidHeader = TRUE,
          dataTableOutput("sizeSampling")
        )
      )
    ), 
    ##### DENSITY TAB #####
    tabItem(tabName="percent",
      conditionalPanel(condition="output.fileUploaded",
        shinydashboard::box(width=12,
          shinydashboard::box(width=3,
            uiOutput("selectFilesDens"),
            column(6,uiOutput("selectXdens")),
            column(6,uiOutput("selectYdens")),
            column(6,sliderInput("probDensity1","Density Bins",min=0.5,max=10,sep=0.1,value=1)),
            column(6,uiOutput("uiSliderProdDensity2")),
            actionButton("activeDens","Create Density"),
            tags$br(),tags$br(),
            uiOutput("densityTableSample"),
            actionButton("applyDensSample","Apply Sampling"),
            actionButton("enrichDensSample","Enrich Sampling")
          ),
          shinydashboard::box(width=9,
            plotOutput("densityPlotSample",height=500,width = 500),
            plotOutput("histPlotX"),
            plotOutput("histPlotY")
          )
        )
      )
    ),
    ##### CLUSTERING TAB ####
    tabItem(tabName = "clusters",
      conditionalPanel(condition = "output.fileUploaded",
        box(id="tabBoxCluster",width = 12,height = 900,
          box(title="Run Clustering", width = 12,collapsible = TRUE,
            column(6,uiOutput("clusterParamsOutput")),
            column(2,textInput("grepParamsClust","Pattern")),
            column(1,actionButton("addGrepParamsClust","Add Pattern")),
            column(1,actionButton("allParamsClust","Add All")),
            column(1,actionButton("clearParamsClust","Clear All"))
          ),
          box(width=3,
            box(width=12,collapsible = TRUE,
              column(12,selectInput("clustering_meth","Clustering Algorithme", choices = c("k-means","FlowSOM","RPhenograph","CLARA","HClust","flowPeaks"))),
              conditionalPanel(
                condition = "input.clustering_meth == 'flowPeaks'",
                column(12, sliderInput("tol","Tolerance",min=0, max=1, value=0.1))
              ),
              conditionalPanel(
                condition="input.clustering_meth != 'FlowSOM'&&input.clustering_meth != 'RPhenograph'&&input.clustering_meth",
                column(12,numericInput("ncluster","Cluster Number",value=20))
              ),
              conditionalPanel(
                condition="input.clustering_meth == 'RPhenograph'",
                column(12,numericInput("knn","knn",value=30,min = 5,max=100))
              ),
              conditionalPanel(
                condition="input.clustering_meth == 'FlowSOM'",
                column(6,numericInput("Xgrid","X grid",value=20)),
                column(6,numericInput("Ygrid","Y grid",value=20)),
                column(4,numericInput("maxCluster","Clusters",value=20)),
                column(8,selectInput("methodsMeta","Methods",
                  choices=c("None","metaClustering_consensus","metaClustering_hclust","metaClustering_kmeans","metaClustering_som"),selected="None"))
              ),
              column(12, actionButton("runClustering","Run Clustering"))
            ),
            conditionalPanel(
              condition = "output.fileClustered",
              box(width=12,collapsible = TRUE,
                checkboxInput("compareClust","Compare two files", value=FALSE),
                uiOutput("IDClusteredPreviewxOutput"),
                conditionalPanel(condition = "input.compareClust",uiOutput("IDClusteredPreviewxOutput2")),
                fluidRow(
                  column(6,uiOutput("selectXCluster")),
                  column(6,uiOutput("selectYCluster"))
                ),
                # conditionalPanel(condition= "input.clustering_meth == 'HClust'")
                actionButton("addClusterID","Add Cluster ID in FCS")
              )
            )
          ),
          conditionalPanel(
            condition = "output.fileClustered",
            box(title="View Clustering Result",width=9,
              column(6,
                plotOutput("plotPreviewCluster",height = 500,width = 500)
                # downloadButton("clusterCenters")
              ),
              conditionalPanel(
                condition = "input.compareClust",
                column(6,
                  plotOutput("plotPreviewCluster2",height = 500,width = 500)
                )
              )
            )
          )
        )
      )
    ),
    ##### CLUSTERING ANALYSIS #####
    tabItem(tabName="clustersAnalysis",
      conditionalPanel(condition = "output.fileUploaded",
        box(width=9,
          box(width=12,collapsible = T,title="Cluster settings",
            column(2,uiOutput("uiSelectClusterId")),
            column(10,uiOutput("uiSelectClusterIdUnique"))
          ),
          box(width=12, collapsible = T,title="Parameters settings",
            column(6,uiOutput("uiSelectParamsHeatmap")),
            column(1,textInput("grepParamsClustAnal","Pattern")),
            column(2,actionButton("addGrepParamsClustAnal","Add Pattern")),
            column(1,actionButton("allParamsClustAnal","Add All")),
            column(1,actionButton("clearParamsClustAnal","Clear All"))
          ),
          tabBox(id="analyse_cluster",width = 12,
            tabPanel("HeatMap",
              column(2,actionButton("runHcl","Run")),    
              uiOutput("heatmapCluster",height = "auto",width = "100%")
            ),
            tabPanel("Markers Density",
              # column(2,uiOutput("uiMultiPlotXParam")),
              # column(2,uiOutput("uiMultiPlotYParam")),
              column(12,uiOutput("uifileViewMDP")),
              uiOutput("markerDensityPlot")
            )
          )
        ),
        box(width = 3,
          column(6,uiOutput("uiClAnalXparam")),
          column(6,uiOutput("uiClAnalYparam")),
          # column(4,uiOutput("uiClAnalClParam")),
          column(12,uiOutput("useViewPlot"))#,width = "450px",height = "450px"))
        )
      )
    ),
    ##### DESCRIPTIVE STAT TAB #####
    tabItem(tabName = "mfi",
      column(12,
        conditionalPanel(
          condition = "output.fileUploaded",
          tabBox(id="mfi2",width=12,title="",
            tabPanel("Compute Descriptive Stats",
              fluidRow(shinydashboard::box( width=12, collapsible=FALSE,
                column(4,checkboxGroupInput("stats","Export stats",choices = c("Mean","Median","Mode","SD","Min","Max"),
                  selected = c("Mean","Median","Mode","SD","Min","Max"),inline = TRUE
                )),
                column(4,selectInput("quantiles","Quantiles Export",choices=c(5,10,15,20,25,30,35,40,60,65,70,75,80,85,90,95),selected=c(5,25,75,95),multiple=TRUE)),
                column(4,radioButtons("bindType","Bind Type",choices=c("Row","Column"),selected="Row",inline = TRUE)),
                uiOutput("mfiMarkerUI"),
                actionButton("computeMFI","Compute"),
                downloadButton("ddlMFItable"),
                box(width=12,dataTableOutput("mfiTable"),style = "width:100%;overflow-x: scroll;font-size:12px;")
              )
            )),
            tabPanel("Create FCS from statistique",
              fluidRow(box(width=12, collapsible=FALSE,
                fluidRow(column(2,uiOutput("uiParamsSelect")),
                  column(8,uiOutput("uiMarkerSelect")),
                  column(1,actionButton("addAllCreateFCS","Add All")),
                  column(1,actionButton("clearAllCreateFCS","Clear List"))),
                fluidRow(column(6,checkboxGroupInput("statsCreate","Export stats",choices = c("Mean","Median","Mode","SD","Min","Max","size"),
                    selected = c("Mean","size"),inline = TRUE)),actionButton("computeCreate","Compute")
                ),
                fluidRow(
                  downloadButton("ddlCSVCreateFromClsuter","Download CSV"),
                  column(12, rHandsontableOutput("previewTableCreateFCS"))
                )
              )
            ))
          )
        )
      )
    ),
    ##### DIFERENTIAL STAT TAB #####
    tabItem(tabName="difstats",
      conditionalPanel(condition = "output.fileUploaded",
      box(width=12,
        box(width=3,
          column(12,selectInput("annotGroup","Annot",choices = NULL,multiple = F)),
          column(6,uiOutput("uiGroup1")),
          column(6,uiOutput('uiGroup2')),
          conditionalPanel(condition= "input.group1",
            column(12,uiOutput("uiAnnotPoint")),
            column(12,uiOutput("uiVariable")),
            column(12,textInput("grepVolcano","Pattern")),
            column(4,actionButton("addAllVolcano","Add All")),
            column(4,actionButton("clearAllVolcano","Clear All")),
            column(4,actionButton("addGrepVolcano","Add PAttern"))
          ),
          tags$br(),tags$br(),tags$br(),tags$br(),
          column(12,actionButton("plotVolcano","Plot"))
        ),
        box(width = 9,
          plotlyOutput("volcanoPlot",width=700,height = 700)
        )
      ))
    ),
    ##### PCA TAB ####
    tabItem(tabName = "pca",
      conditionalPanel(
        condition = "output.fileUploaded",
        shinydashboard::box(title="Dimension Reduction in FCS by PCA" ,collapsible= FALSE, width = 12,status="success", solidHeader=TRUE,
          fluidRow(
            column(6,uiOutput("selectPCAParams")),
            column(1,actionButton("runPCA","Run")),
            column(1,actionButton("addAllPCADim","Add All")),
            column(1,actionButton("clearAllPCADim","Clear List")),
            column(2,textInput("grepPCA","Pattern")),
            column(1,actionButton("addGrepPCA","Add Pattern"))
          ),
            conditionalPanel(condition = "output.fileReducPCA",
            shinydashboard::box(width=2,
              uiOutput("selectPCA"),
              checkboxInput("comparePCA","Compare two files", value=FALSE),
              conditionalPanel(
                condition = "input.comparePCA",
                uiOutput("selectPCA2")
              ),
              uiOutput("ZPCAPlot"),
              uiOutput("XPCAPlot"),
              uiOutput("YPCAPlot"),
              uiOutput("selectPCAToAdd"),
              actionButton("addPCA","Add PCA dim")
            ),
            shinydashboard::box(collapsible= TRUE, width = 10,
              column(6,plotOutput("PCAPlot",width=600, height = 600)),
              conditionalPanel(
                condition = "input.comparePCA",
                column(6, plotOutput("PCAPlot2",width=600, height = 600))
              )
            ),
            shinydashboard::box(collapsible= TRUE, width = 10,
              box(width=12,column(4,uiOutput("uiInfoPCA")),
              conditionalPanel(condition="input.infoPCAID == 'Contribution'",
              column(3,uiOutput("uiContribAxes")),column(3,uiOutput("uiContribTop")))),                
              column(6,plotOutput("plotInfo1",width=600, height = 400)),
              conditionalPanel(
                condition = "input.comparePCA",
                column(6,plotOutput("plotInfo2",width=600, height = 400))
              )
            )
          )
        )
      )
    ),
    ##### UNLINEAR TAB #####
    tabItem(tabName = "umap",
      conditionalPanel(
        condition = "output.fileUploaded",
        shinydashboard::box(title="Dimension Reduction in FCS by UMAP" ,collapsible= FALSE, width = 12,status="success", solidHeader=TRUE,
          fluidRow(
            column(1,selectInput("unlineartReduc","Select Meth",choices=c("UMAP","R-TSNE","IsoMap","DiffusionMap","EmbedSOM"))),
            column(5,uiOutput("selectUMAPParams")),
            column(1,actionButton("runUMAP","Run")),
            column(1,actionButton("addAllUMAP","Add All")),
            column(1,actionButton("clearAllUMAP","Clear")),
            column(2,textInput("grepUMAP","Pattern")),
            column(1,actionButton("addGrepUMAP","Add Pattern")),
            conditionalPanel(condition="input.unlineartReduc == 'R-TSNE'",
              column(2,numericInput("iter","Iter",1000)),
              column(2,numericInput("theta","Theta",0.4)),
              column(2,numericInput("knn","Perplexity",30)),
              column(2,numericInput("eta","Learning Rate",200))
            ),
            conditionalPanel(condition="input.unlineartReduc == 'UMAP'",
              column(2,checkboxInput("advOptUmap","Adv Set.",value=FALSE)),
              conditionalPanel(condition="input.advOptUmap",
                column(2,numericInput("n_neighbors","Num_neighbot",value=15)),
                column(2,selectInput("metric","Metric",choices=c("euclidean","cosine","manhattan","hamming")))
              )
            )
          ),
          conditionalPanel(condition = "output.fileReducUMAP",
            shinydashboard::box(width=2,
              uiOutput("selectUMAP"),
              checkboxInput("compareUMAP","Compare two files", value=FALSE),
              conditionalPanel(
                condition = "input.compareUMAP",
                uiOutput("selectUMAP2")
              ),
              uiOutput("ZumapPlot"),
              uiOutput("XumapPlot"),
              uiOutput("YumapPlot"),
              # downloadButton("ddlUMAP","Download CSV"),tags$br(),tags$br(),
              uiOutput("selectUMAPToAdd"),
              actionButton("addUMAP","Add New Dim")
            ),
            shinydashboard::box(collapsible= TRUE, width = 10,
              column(6,plotOutput("UMAPPlot",width=700, height = 500)),
              conditionalPanel(
                condition = "input.compareUMAP",
                column(6, plotOutput("UMAPPlot2",width=700, height = 500))
              )
            )
          )
        )
      )
    ),
    ##### ONESENSE #####
    tabItem(tabName="onesens",
      conditionalPanel(condition = "output.fileUploaded",
        tabBox(width=12,height=900,
          tabPanel("Input",
            fluidRow(
              column(2,uiOutput("uiSelectOnsenseFile")),
              column(2,selectInput("oneSENSEMeth","Method",choices = c("Rtsne","UMAP"))),
              column(2,actionButton("firstMarkers","1. Display Markers")),
              column(1,actionButton("updatenamescsv","2. Confirm")),
              column(1,actionButton("submit", "3. Run OneSENSE")),
              column(1,numericInput("num", "Ceiling Input", value = 500)),
              column(1,numericInput("bins", "Bins", value = 250))
            ),
            fluidRow(
              column(8,selectInput("markers1","1st Category", NULL,multiple=TRUE,width = 1200)),
              column(1,actionButton("addAllMarkers1","Add All")),
              column(1,actionButton("clearMarkers1","Clear")),
              column(1,textInput("grepMarkers1","Pattern")),
              column(1,actionButton("addAllGrepMarkers1","Add Pattern"))
            ),
            fluidRow(
              column(8,selectInput("markers2","2st Category", NULL,multiple=TRUE,width = 1200)),
              column(1,actionButton("addAllMarkers2","Add All")),
              column(1,actionButton("clearMarkers2","Clear")),
              column(1,textInput("grepMarkers2","Pattern")),
              column(1,actionButton("addAllGrepMarkers2","Add Pattern"))
            )
            #checkboxGroupInput("markers3","Select 3rd category markers (Optional)", "")
          ),
          tabPanel("Summary Plot",
            column(12,actionButton("enrichOneSENSE","Add OneSENSE dim")),
            column(6,plotlyOutput("testOneSense2", width = 450,height = 400)),
            column(6,plotlyOutput("testOneSense3",width=400,height = 400)),
            column(6,plotlyOutput("testOneSense4", width = 450,height = 450))
            # column(6,plotOutput("testOneSense1",width=400,height = 400))
          ),
          tabPanel("Dot Plot",
            plotlyOutput("dotPlot", width = 800,height = 800)),
          tabPanel('Heatmap 1st category',
            column(3,downloadButton("downloadHeatMap1","Download Heatmap1")),
            plotlyOutput("heatmap1", width = 1200,height = 800)),
          tabPanel('Heatmap 2st Category',
            column(3,downloadButton("downloadHeatMap2","Download Heatmap2")),
            plotlyOutput("heatmap2", width = 1200,height = 800)),
          tabPanel("Coordinate Selection",
            box(width=12,column(6,actionButton("submit2", "Coordinate Selection")),
            uiOutput("markerXdropdown"),
            uiOutput("markerYdropdown"),
            tableOutput("table"),
            plotOutput("plot", click = "plot_click")
          ),
          tabPanel("Frequency Heatplot",
            column(6,actionButton("submit3", "Frequency heatmap")))
          )
        )
      )
    ),
    ##### MANUAL ANNOTATION #####
    tabItem(tabName="manuannot",
      conditionalPanel(condition = "output.fileUploaded",
      box(width=12,
        box(width=6,
          column(4,selectInput("xAnnotPlot","X Param",choices=NULL)),
          column(4,selectInput("yAnnotPlot","Y Param",choices=NULL)),
          column(4,selectInput("zAnnotPlot","Z Param",choices=NULL)),
          box(width=12,plotlyOutput("annotPlot",width=750,height=650))
        ),
        box(width=4,
          column(3,selectInput("xAnnotPlotSub","X Param",choices=NULL)),
          column(3,selectInput("yAnnotPlotSub","Y Param",choices=NULL)),
          column(6,selectInput("zAnnotPlotSub","Z Param",choices=NULL)),
          column(12,plotlyOutput("subAnnotPlot",width=450,height = 400)),
          column(6,textInput("annotName","Annotation Name",value="manual_Annot")),
          column(3,actionButton("addAnnotName","Add")),
          column(3,actionButton("enrichAnnotName","Enrich FCS"))
        ),
        box(width=2,
          column(12,uiOutput("uiAddPop")),
          tableOutput("annotTable")
        )
      ))
    ),
    ##### MANUAL GATING #####
    tabItem(tabName="gate",
      conditionalPanel(condition="output.fileUploaded",
        box(width=3,
            box(width=12, collapsible = T,title="How to use",collapsed = T,
              tags$ul(
                tags$li("Pour r??aliser un petit gating manuel (maximum 5 ??tapes) selectionner le fichier sur lequels 
                        vous allez travailler puis les param??tre X et Y et par defaut l'??tape Root"),
                tags$li("Desinner la gate en r??alisant des clique sur le plot a droite. Attention en fonction de la taille
                        de votre jeux de donn??es les point peuvent prendre du temps a apparaitres. Chaque point seras relier
                        au point pr??c??dent pour dessiner un polygon"),
                tags$li("Une fois terminer, remplisser votre Gate Name et cliquer sur Gate : Cela va automatiquement relier 
                  le premier et le dernier point afin d'avoir une figure terminer"),
                tags$li("Votre sous population apparait dans la liste Select Step et vous pourrez de nouveau travailler dessus
                         ou bien refaire une gate sur ??tape precedente"),
                tags$li("A droite vous avez un r??sumer des population en taille et en frequence par rapport au all events. Pour appliquer 
                        les polygons a l'ensemble des fichier aller dans l'onglet View et laisser l'outil vous g??n??rer les plot par fichier. Une ligne ??gale un fichier")
            )),
          actionButton("refresh","Refresh"),
          uiOutput("selectGateFile"),
          uiOutput("selectGateStep"),
          uiOutput("selectXGate"),
          uiOutput("selectYGate"),
          textInput("stepName","Gate Name"),
          actionButton("createGate","Gate"),
          actionButton("deleteGate","Delete")
        ),
        box(width = 9,
          column(6,plotOutput("plotGate",click = "point",dblclick = "end",width = "600px",height = "600px")),
          column(6,dataTableOutput("tableGateManual"))
        )
      )
    ),
    tabItem(tabName="gateview",
      conditionalPanel(condition="output.fileUploaded",
        box(width=12,
          column(4,uiOutput("selectGateStepView")),
          column(2,actionButton("loadSubPop","Load")),
          column(2,actionButton("ddlSubPop","Download")),
          column(2,textOutput("linkBigSubPop")),
          column(12,uiOutput("gateViewStep"))
        
        )
      )
    ),
    ##### AUTO ANNOTATION #####
    tabItem(tabName="autoannot",
      conditionalPanel(condition="output.fileUploaded",
        box(width=9,title="Select Your enrichment"),
        box(width=3,title="Upload your Reference",
          fileInput("ref_file","Upload FCS reference",accept = c(".fcs"),multiple = TRUE),
          column(2,checkboxInput("ref_com","Comp.",value=TRUE)),
          column(5,selectInput("trans_ref",label = "Transf.",choices =c("Arcsinh","Logicle","None"))),
          column(5,numericInput("trans_ref_args",label="Args",value=5)),
          uiOutput("uiMarkersTransRef"),
          tableOutput("listRef"),
          actionButton("ref_prepross","Preprocess")
        )
      )
    )
  )
)
##### UI ######
ui <- dashboardPage(
  header,
  sidebar,
  body
)
  









