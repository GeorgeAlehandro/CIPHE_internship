##### shiny package #####
library(dplyr)
library(cipheinfinityFlow)
options(shiny.maxRequestSize = 100000*1024^2)
#### PRE-REQ
source("C:\\Users\\gsaad\\Desktop\\libraryCIPHE\\help.func.R")
source("C:\\Users\\gsaad\\Desktop\\libraryCIPHE\\cleaning.func.R")
source("C:\\Users\\gsaad\\Desktop\\libraryCIPHE\\gating.func.R")

roots1 <- c(w="D:/DATA/")
roots2 <- c(w="D:/METADATA/")
roots3 <- c(w="D:/CLEAN/")

no.transform <- c("FSC-A","SSC-A","FSC-H","SSC-H","FSC-W","SSC-W","Time","Flag","FlagDens")
server <- function(session, input, output) {


##### GLOBAL #####################################

  values <- reactiveValues(
    backbone_specification=NULL,
    fcs_files = NULL,
    number_files = NULL,
    path_to_raw_files = NULL,
    infinity_markers_specification_all_wells = NULL,
    filechosen = NULL,
    flow.frames = NULL,
    flowAI = NULL,
    flowClean = NULL,
    file_selected = NULL,
    average_number_events = NULL,
    exp.name = NULL,
    raw = list(),
    flow.frames=list(),
    flow.frames.quality = list(),
    clean.flow.frames = list(),
    metadata = NULL,
    params = NULL,
    markers = NULL,
    trans = "none",
    spill=NULL,
    cleaning.data = list(),
    concatenate.fcs= NULL,
    random.chosen = FALSE,
    buttons.created = FALSE,
    fcs.output = NULL,
    plot.object = NULL,
    list.param = NULL,
    buttons = NULL,
    list.arg = NULL,
    new.value.param = NULL,
    threshold = NULL,
    list.func = NULL,
    logical = FALSE,
    logical2 = FALSE,
    gating.all = NULL,
    validation = FALSE,
    logical.clean = FALSE,
    validation2 = FALSE
  )

  checkDimName <- function(dim,fcs,button_id){
    output$dimName <- renderUI({
      objectDimName <- lapply(c(1:length(dim)),function(i){
        textInput(paste0("dim_",dim[[i]]),"Name of dim",value=dim[[i]])
      })
      return(do.call(tagList,objectDimName))
    })
    modalDialog(
      column(6,uiOutput("dimName")),
      actionButton(button_id,"Enrich")
    )
  }
  ###Disabling buttons on Loading page
   shinyjs::disable('confirm')
   shinyjs::disable('confirm_infinity')
##### DOWNLOAD IN MENUITEM #######################

  output$downloadData <- downloadHandler (
    filename = function() {
      paste("data_", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      progress <- Progress$new()
      flow.frames <- values$flow.frames
      root <- getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
      fs <- c()
      progress$set(message="Zip files", value=0.5)
      for(i in c(1:length(flow.frames))){
        name <- names(values$flow.frames)[i]
        fcs <- flow.frames[[i]]
        path <- gsub(".fcs$",".fcs",name)
        fs <- c(fs, path)
       # fcs <- updateKeywords(fcs)
        write.FCS(fcs, path)
      }
      zip(zipfile=file, files=fs)
      setwd(root)
      progress$close()
    }
  )
  output$ex_out <- renderPrint({
    print(values$log)
  })

  output$bigFileUploaded <- reactive({
    return(object.size(values$flow.frames)/1000000>10)
  })

  output$smallFileUploaded <- reactive({
    return(object.size(values$flow.frames)/1000000<10)
  })

  outputOptions(output, "bigFileUploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "smallFileUploaded", suspendWhenHidden = FALSE)
  ##TRANSFORMED TO clean.flow.frames
  observeEvent(input$bigDDl,{
  
    if(Sys.info()[1] == "Linux"){
      progress <- Progress$new()
      flow.frames <- values$clean.flow.frames
      root <- getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
      fs <- c()

      list.files.remove <- list.files(".",pattern=".fcs",recursive = F,full.names=F)
      for(t in list.files.remove){
        print(t)
        file.remove(t)
      }
      progress$set(message="Zip files", value=0.5)
      for(i in c(1:length(flow.frames))){
        name <- names(values$clean.flow.frames)[i]
        fcs <- flow.frames[[i]]
        path <- gsub(".fcs$",".fcs",name)
        fs <- c(fs, path)
        #fcs <- updateKeywords(fcs)
        write.FCS(fcs, path)
      }

      zip(zipfile="output.zip", files=fs)
      setwd(root)
      values$path.big.file <- paste0(tmpdir,"/output.zip")
      values$path.big.file <- str_replace(values$path.big.file,"/media/data/html","http://10.71.1.22/")
      progress$close()

      output$lnkBigDDl <- renderText({
        return(paste0("<a href='",values$path.big.file,"'>Download</a>"))
      })
    }else if(Sys.info()[1]=="Windows"){
      print('FOUND THE CURRENT FLOW FRAMES')
      print(values$flow.frames)
      print('FOUND THE Clean FLOW FRAMES')
      print(values$clean.flow.frames)
      #dir <- chooseDir()
      dir <- choose.dir()
      dir <- gsub("\\\\","/",dir)
      progress <- Progress$new()
      flow.frames <- values$clean.flow.frames
      fs <- c()
      for(i in c(1:length(flow.frames))){
        name <- paste0(dir,"/",names(values$clean.flow.frames)[i])
        fcs <- flow.frames[[i]]
        path <- gsub(".fcs$",".fcs",name)
        fs <- c(fs, path)
       # fcs <- updateKeywords(fcs)
        write.FCS(fcs, path)
      }
      progress$close()
      showNotification(ui=paste0("Your FCS is write at : ",dir),type="message")
    }
  })
##### CHOICE of FCS files ######
  observeEvent(input$input_fcs,{
    if(!is.null(input$input_fcs) && input$input_fcs != ""){
      # shinyjs::disable("input_fcs")
      progress <- Progress$new()
      progress$set(message="Read Select file : ",value=1)
      print('input_fcs')
      print(input$input_fcs)
      values$exp.name <- basename(input$input_fcs$datapath)
      #listFCS <- list.files(paste0(input$input_fcs), recursive=FALSE, full.names = TRUE)
      listFCS <- input$input_fcs$datapath
      print('listFCS')
      print(listFCS)
    #  values$flow.frames <- NULL
      A <- read.FCS(listFCS[1])
      print('A')
      print(A)
      print('values$flow.frames BEFORE lapply')
      print(values$flow.frames)
      listFCS <- as.vector(listFCS)
      print('new list fcs')
      print(listFCS)
      print(typeof(listFCS))
      print(class(listFCS))
      # progress$set(message=listFCS,value=1)
      #values$flow.frames <- lapply(listFCS, function(i){a <- read.FCS(i); values$flow.frames <- list(a); print(values$flow.frames)})
      for (file in listFCS){
        one_read <- read.FCS(file);
       # values$flow.frames <- list(values$flow.frames, one_read)
        values$flow.frames <-do.call(c, list(values$flow.frames, one_read))
        
      }
      # 
      # if(file.exists(paste0(roots2,basename(input$input_fcs$datapath),"_MetaData.csv"))) {
      #   values$metadata <- read.csv(paste0(roots2,basename(input$input_fcs$datapath),"_MetaData.csv"),header = TRUE,check.names = TRUE)
      # } else {
      #   progress <- Progress$new()
      #   progress$set(message = paste0("Verify your data", paste0(roots2,basename(input$input_fcs$datapath),"_MetaData.csv")," doesn't exist!"))
      # }
      print('here basename')
      print(basename(listFCS))
      print('input_fcs')
      print(input$input_fcs$name)
      names(values$flow.frames) <- input$input_fcs$name
      print('1')
      print(values$flow.frames)
     # values$flow.frames <- values$flow.frames[values$metadata[,"filesnames"]]
      print('2')
      print(values$flow.frames)
      View(values$flow.frames)
      # names(values$flow.frames) <- values$metadata[,"filesnames"]
      names <- as.vector(pData(values$flow.frames[[1]]@parameters)[,2])
      params <- as.vector(pData(values$flow.frames[[1]]@parameters)[,1])
      names(params) <- params
      names(params)[]
      
      progress$set(message=1)
      values$params <- NULL
      values$params <- as.vector(get.names.CIPHE(values$flow.frames[[1]]))
      values$name.params <- NULL
      values$names.params <- names
     # values$spill <- values$flow.frames[[1]]@description[[found.spill.CIPHE(flow.frames)[[1]]]]
      values$list.params <- NULL
      values$list.params <- as.vector(get.markers.CIPHE(values$flow.frames[[1]]))
      values$list.params2 <- NULL
      values$list.params2 <- as.vector(get.markers.2.CIPHE(values$flow.frames[[1]]))
      progress$set(message=2)
      output$select1 <- renderUI({selectInput("transMethod", "Transform Method",choices=c("logicle","estimate","none"))})
      output$select2 <- renderUI({selectInput("markersY", "Y markers", choices = values$list.params2)})
      output$select3 <- renderUI({selectInput("viewFiles", "Files", choices = names(values$flow.frames), selected = 1)})
      output$button1 <- renderUI({actionButton("hide1", "Hide")}) 
      output$button2 <- renderUI({actionButton("show1", "Show")}) 
      output$button3 <- renderUI({actionButton("hide2", "Hide")}) 
      output$button4 <- renderUI({actionButton("show2", "Show")}) 
      values$raw <- values$flow.frames
      
      progress$close()
      
      
      # shinyjs::enable("input_fcs")
      # Initialize if empty
    } else {
      values$exp.name <- NULL
      values$flow.frames <- list()
      values$metadata <- NULL
    }
  })
##### CLEAN ####################
  observe({
    if(is.null(values$params)) return(NULL)
    if(length(values$flow.frames)<1) return(NULL)
    list.param <- values$params
    values$validation <- TRUE
    output$clean.live1 <- renderUI({
      selectizeInput("clean.param.live1", "Live Parameter X", choices = 'FSC-A',
                     selected = "Live/Dead", multiple = FALSE)
    })
    output$clean.live2 <- renderUI({
      selectizeInput("clean.param.live2", "Live Parameter Y", choices = 'UV-BUV496-A',
                     selected = "Live/Dead", multiple = FALSE)
    })
    #values$list.params2
    output$clean.size1 <- renderUI({
      selectizeInput("clean.param.size1", "Size Parameter X", choices = 'FSC-A',
                     selected = "FSC-A")
    })
    output$clean.size2 <- renderUI({
      selectizeInput("clean.param.size2", "Size Parameter Y", choices = 'SSC-A',
                     selected = "SSC-A")
    })
    output$clean.singletsFCS1 <- renderUI({
      selectizeInput("clean.param.singletFSC1", "X Param SingletFSC", choices = 'FSC-A',
                     selected = "FSC-A")
    })

    output$clean.singletsFSC2 <- renderUI({
      selectizeInput("clean.param.singletFSC2", "Y Param SingletFSC", choices = 'FSC-H',
                     selected = "FSC-H")
    })

    output$clean.singletsSSC1 <- renderUI({
      selectizeInput("clean.param.singletSSC1", "X Param SingletSSC", choices = 'SSC-W',
                     selected = "SSC-W")
    })

    output$clean.singletsSSC2 <- renderUI({
      selectizeInput("clean.param.singletSSC2", "Y Param SingletSSC", choices = 'SSC-H',
                     selected = "SSC-H")
    })

  })
  
  observeEvent(input$gating,{
    progress <- Progress$new()
    # if(length(values$flow.frames.quality) == 0) {
    flow.frames <- values$flow.frames
    # } else {
    # flow.frames <- values$flow.frames.quality
    # }
    progress$set(message = "Gating in progress",value = 0.33)
    
    step <- c(1,1,1,1)
    cleaning.x.param.live <-input$clean.param.live1
    cleaning.x.min.live<- input$clean.gate.live1[1]
    cleaning.x.max.live <-input$clean.gate.live1[2]
    cleaning.y.min.live<- input$clean.gate.live2[1]
    cleaning.y.max.live <-input$clean.gate.live2[2]
    cleaning.y.param.live <-input$clean.param.live2
    #cleaning.param.size<- input$clean.param.size
    cleaning.x.param.size<- input$clean.param.size1
    cleaning.y.param.size<- input$clean.param.size2
    cleaning.x.min.size <-input$clean.gate.size1[1]
    cleaning.x.max.size <-input$clean.gate.size1[2]
    cleaning.y.min.size <-input$clean.gate.size2[1]
    cleaning.y.max.size <-input$clean.gate.size2[2]
    cleaning.x.param.singletFSC <- input$clean.param.singletFSC1
    cleaning.y.param.singletFSC <- input$clean.param.singletFSC2
    cleaning.x.param.singletSSC <- input$clean.param.singletSSC1
    cleaning.y.param.singletSSC <- input$clean.param.singletSSC2
    cleaning.singletSSC.bot.left.corner <- c(input$x.border.left,input$y.border.left)
    cleaning.singletSSC.top.right.corner <- c(input$x.border.right, input$y.border.right)
    print('t2akad men l flowframes')
    print(values$flow.frames)
    
     values$cleaning.data <- lapply(flow.frames, function(x) {
   # values$cleaning.data <- mclapply(flow.frames,  mc.cores = 1, FUN=function(x) {  
      # source("cleaning.func.R")
      return(cleaning.CIPHE(flow.frame = x, cleaning.step = step,
                            cleaning.x.param.live = cleaning.x.param.live,
                            cleaning.y.param.live = cleaning.y.param.live,
                            cleaning.x.min.live = cleaning.x.min.live, cleaning.x.max.live = cleaning.x.max.live,
                            cleaning.y.min.live = cleaning.y.min.live, cleaning.y.max.live = cleaning.y.max.live,
                            cleaning.x.min.size = cleaning.x.min.size, cleaning.x.max.size = cleaning.x.max.size,
                            cleaning.y.min.size = cleaning.y.min.size, cleaning.y.max.size = cleaning.y.max.size,
                            cleaning.x.param.size = cleaning.x.param.size,
                            cleaning.y.param.size = cleaning.y.param.size,
                            cleaning.x.param.singletFSC = cleaning.x.param.singletFSC,
                            cleaning.y.param.singletFSC = cleaning.y.param.singletFSC,
                            cleaning.x.param.singletSSC = cleaning.x.param.singletSSC,
                            cleaning.y.param.singletSSC = cleaning.y.param.singletSSC,
                            cleaning.singletSSC.bot.left.corner = cleaning.singletSSC.bot.left.corner,
                            cleaning.singletSSC.top.right.corner = cleaning.singletSSC.top.right.corner)
      )
    })
    print(values$cleaning.data)
    progress$set(message = length(values$cleaning.data),value = 0.66)
    
    names(values$cleaning.data) <- names(values$flow.frames)
    
    values$clean.flow.frames <- lapply(values$cleaning.data, function(x){
      return(x$singlets.SSC$flowFrame)
    })
    progress$set(message = length(values$clean.flow.frames),value = 0.99)
    progress$set(message=names(values$flow.frames))
    names(values$clean.flow.frames) <- names(values$flow.frames)
    progress$close()
  })
  
  # observe({
  #   if(is.null(values$flow.frames)) return(NULL)
  #   output$selectViewPlot <- renderUI({
  #     selectInput("viewPlot","Select View",choices = names(values$flow.frames))
  #   })
  # })
  
  observeEvent(input$viewPlot,{
    flow.frames <- values$flow.frames
    i <- input$viewPlot
    
    output$cleaningPlot <-  renderPlot({
      plot.cleaning.CIPHE(flow.frames[[i]],
                          values$cleaning.data[[i]],
                          #param.live = input$clean.param.live,
                          param.live = c(input$clean.param.live1,input$clean.param.live2),
                          param.size = c(input$clean.param.size1,input$clean.param.size2),
                          param.singletSSC = c(input$clean.param.singletSSC1,input$clean.param.singletSSC2),
                          param.singletFSC = c(input$clean.param.singletFSC1,input$clean.param.singletFSC2))
    }, outputArgs = list(width = 1100,height = 300))
  })
  
  
  
  
  observe({
    if(is.null(values$flow.frames)) return(NULL)
    output$selectViewPlot <- renderUI({
      selectInput("viewPlot","Select View",choices = names(values$flow.frames))
    })
  })

  # observeEvent(input$viewPlot,{
  #   flow.frames <- values$flow.frames
  #   i <- input$viewPlot
  #   
  #   output$cleaningPlot <-  renderPlot({
  #     plot.cleaning.CIPHE(flow.frames[[i]],
  #                         values$cleaning.data[[i]],
  #                         param.live = input$clean.param.live,
  #                         param.size = input$clean.param.size,
  #                         param.singletSSC = c(input$clean.param.singletSSC1,input$clean.param.singletSSC2),
  #                         param.singletFSC = c(input$clean.param.singletFSC1,input$clean.param.singletFSC2))
  #   }, outputArgs = list(width = 1100,height = 300))
  # })
##### UPLOAD FCS PANEL ###########################
  observe({
    print("a")
    if(!is.null(values$flow.frames)){
      print("d")
      shinyjs::disable(id = "fcs_input")
    } else {
      print("e")
      shinyjs::enable(id="fcs_input")
    }
  })
  volumes = getVolumes()

#####################################
  ### For each value selected by the slider, observe a certain event
  #### Returns the destination  but should be fixed.
  ### The value should be returned, but also showed as print as a good format
  # that can also be taken by infinity flow
  ####Function used to reset all the inputs when moving through slider
  reset <- function(){
    values$file_selected <- NULL
  }
#####################################SUBMIT
  observeEvent(input$input_fcs,{
    values$files_selected<- input$files_to_clean
    output$filechosen <- renderText({
      as.character(values$files_selected$name)
    })})
  observeEvent(input$submit,
               {
                   values$path_to_raw_files <- (as.vector(values$files_selected$datapath))
                   values$number_files <- length(values$path_to_raw_files)
                   print(values$number_files)
    

    }


          )

########################CONFIRM EXPLORATORY AND BACKGROUND SELECTION
  ## Checks if there are NA values --> throws an error if value is NA
  ## Then if all the values are not NA --> checks if at least one of the data
  ## acquisition is 'exploratory'


 #Initialized customized table button behavior

####VIEWING PDF
  #UMAP not background corrected
  observeEvent(input$generate_umap, {
    print('loading done')
  output$plot_umap <-
    renderUI({
    tags$iframe(style="height:1000px; width:260%", src="umap_plot_annotated.pdf")
  })
  })
  #UMAP background corrected
  observeEvent(input$generate_umap_corrected, {
    print('loading done')
    output$plot_umap_corrected <-
      renderUI({
        tags$iframe(style="height:1000px; width:260%", src="umap_plot_annotated_backgroundcorrected.pdf")
      })
  })
#########################
  observeEvent(input$refresh_input,{
    print("refresh")
    values$flow.frames = NULL
    values$flowAI = NULL
    values$flowClean = NULL
    values$names.files = c()
    values$means.list = c()
    values$params = list()
    values$params.table = c()
    values$template = NULL
    values$clustering.groups = c()
    values$check.names = FALSE
    values$clustering = NULL
    values$centres = list()
    values$csv = list()
    values$mfi.table = NULL
    values$create.table = NULL
    values$pca.reduc = list()
    values$umap.reduc = NULL
    values$sampling = list()
    values$density = NULL
    values$names.csv = NULL
    values$res1.onesense = NULL
    values$list.plot.output = NULL
    values$manu.annot = NULL
    values$path.big.file = NULL
    values$ref.flow.frames = NULL
    values$ref.table  = NULL
    values$step.ff = list()
    values$step.gates = list()
    values$polygon.temp = NULL
    values$gate.strat = NULL
    values$log = matrix(nrow = 1,ncol = 1,data = "Upload", dimnames = list(NULL,"Log"))
    values$package.table = NULL
    showNotification("all_refresh input done", type = "message")
  })
}
