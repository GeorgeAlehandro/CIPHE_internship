##### shiny package #####

options(shiny.maxRequestSize = 100000*1024^2)

server <- function(session, input, output) {

##### GLOBAL #####################################
  
  values <- reactiveValues(
    flow.frames = NULL,
    flowAI = NULL,
    flowClean = NULL,
    names.files = c(),
    means.list = c(),
    params = list(),
    params.table = c(),
    template = NULL,
    clustering.groups = c(),
    check.names = FALSE,
    clustering = NULL,
    centres = list(),
    csv = list(),
    mfi.table = NULL,
    create.table = NULL,
    pca.reduc = list(),
    umap.reduc = NULL,
    sampling = list(),
    density = NULL,
    names.csv = NULL,
    res1.onesense = NULL,
    list.plot.output = NULL,
    manu.annot = NULL,
    path.big.file = NULL,
    ref.flow.frames = NULL,
    ref.table  = NULL,
    step.ff = list(),
    step.gates = list(),
    polygon.temp = NULL,
    gate.strat = NULL,
    log = matrix(nrow = 1,ncol = 1,data = "Upload", dimnames = list(NULL,"Log")),
    package.table = NULL
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
   
##### LIBRARY ####################################
  
  lib.shiny <- c("shinyjs","DT","shinyDND","shinyjqui","shinyBS","BiocManager","devtools","shinyHeatmaply")
  lib.other <- c("Rcpp","cluster","ggplot2","igraph","plyr","reshape","ggridges","locfit","changepoint","plotly","stringr","d3r",
                 "Rtsne","umap","destiny","FactoMineR","factoextra","hypergate","gplots","jsonlite","scales","grDevices","heatmaply")
  lib.bio <- c("flowCore","flowDensity","FlowSOM","flowViz","FlowRepositoryR","flowAI","flowClean","Biobase","GateFinder","scone","flowPeaks")
  lib.git <- c("uwot","EmbedSOM","flowCut","Rphenograph","oneSENSE")
  rep.git <- c("jlmelville","exaexa","jmeskas","JinmiaoChenLab")
  all <- list("shiny"=lib.shiny,"other"=lib.other,"bio"=lib.bio,"git"=lib.git)
  m <- length(all);  i <- 0
  
  progress <- Progress$new()
  taks.lib <- c()
  for(p in c(1:length(all))){
    taks.lib.lib <- c()
    for(l in all[[p]]){
      i <- i+1
      progress$set(message=paste0("Load ",names(all)[[p]]," library"),detail=l,value=i/m)
      c <- library(l, character.only = T,logical.return = T)
      taks.lib.lib <- c(taks.lib.lib,c)
    }
    taks.lib <- c(taks.lib, list(taks.lib.lib))
  }

  output$messageMenu <- renderMenu({
    val.res <- lapply(taks.lib,function(i){
      return(round((length(which(i==TRUE))/length(i))*100,2))
    })
    status <- length(which(unlist(val.res)==100))==length(unlist(val.res))
    if(status){
      badgeStatus = "success"
    } else {
      badgeStatus = "warning"
    }
    dropdownMenu(
      type = "tasks", badgeStatus = badgeStatus,headerText = "Check Install Packages",
      taskItem(value=val.res[[1]], color="green","Shiny Library"),
      taskItem(value=val.res[[2]], color="aqua","General Library"),
      taskItem(value=val.res[[3]], color= "red","Bio Library"),
      taskItem(value=val.res[[4]], color="yellow","Github Library")
    )
  })
  
  output$lib.bio.table <- renderDataTable({
    list <- unlist(taks.lib)
    which(list==TRUE)
    library <- unlist(all)
    c1 <- list
    c2 <- unlist(lapply(c(1:length(library)),function(j){
      if(c1[[j]]){return(packageDescription(library[j],fields = "Version"))}
      return("FALSE")
    }))
    type <- unlist(lapply(c(1:4),function(j){rep(names(all)[[j]],length(all[[j]]))}))
    package.table <- data.frame("Library"=library,"Status"=c1,"Version"=c2,"Type"=type, stringsAsFactors = FALSE)
    values$package.table <- package.table
    return(package.table)
  },options = list(pageLength = 20, info = FALSE, searching = FALSE),rownames=FALSE)
  
  output$selectUpdate <- renderUI({
    print(values$package.table)
    l <- values$package.table[,"Library"]
    print(l)
    selectInput("package","Select package",choices=l, multiple = F)
  })
  
  observeEvent(input$updatePackage,{
    progress <- Progress$new()
    progress$set(message=input$package, value=0.9)
    source <- values$package.table[which(values$package.table[,"Library"]==input$package),"Type"]

    if(source == "other" || source == "shiny"){
      install.packages(l,dependencies = "Depends")
    }
    if(source =="bio"){
      BiocManager::install(l,ask = FALSE,update = FALSE)
    }
    # if(c==FALSE && names(all)[p]=="git"){
    #   devtools::install_github(paste0(rep.git[[p]],"/",lib.git[[p]]),force = TRUE)
    # }
    progress$close()
  })
  

  output$test <- renderMenu({
      sidebarMenu(
        menuItem("Library", tabName = "library", icon = icon("dashboard"))
      )
    })

  progress$set(message="Read source", value=1)
  source("flowAI/auto-qc.R")
  source("function.R")
  source("oneSENSE.R")
  source("isomap.R")
  progress$close()
  
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
        fcs <- updateKeywords(fcs)
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
  
  observeEvent(input$bigDDl,{

    if(Sys.info()[1] == "Linux"){
      progress <- Progress$new()
      flow.frames <- values$flow.frames
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
        name <- names(values$flow.frames)[i]
        fcs <- flow.frames[[i]]
        path <- gsub(".fcs$",".fcs",name)
        fs <- c(fs, path)
        fcs <- updateKeywords(fcs)
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
      dir <- chooseDir()
      dir <- gsub("\\\\","/",dir)
      progress <- Progress$new()
      flow.frames <- values$flow.frames
      fs <- c()
      for(i in c(1:length(flow.frames))){
        name <- paste0(dir,"/",names(values$flow.frames)[i])
        fcs <- flow.frames[[i]]
        path <- gsub(".fcs$",".fcs",name)
        fs <- c(fs, path)
        fcs <- updateKeywords(fcs)
        write.FCS(fcs, path)
      }
      progress$close()
      showNotification(ui=paste0("Your FCS is write at : ",dir),type="message")
    }
  })
  
##### UPLOAD FCS PANEL ###########################

  getParams <- reactive({
    if(is.null(values$flow.frames)) return(NULL)
    data <- as.matrix(pData(values$flow.frames[[1]]@parameters),stringsAsFactors = F)
    labels <- data[,2]
    params <- data[,1]
    labels[which(is.na(labels))] <- colnames(values$flow.frames[[1]])[c(which(is.na(labels)))]
    labels[which(labels=="<NA>")] <- colnames(values$flow.frames[[1]])[c(which(labels=="<NA>"))]
    names(params) <- labels
    values$params <- params
  })
  
  checkAnnot <- function(fcs,params,th){
    if(length(unique(fcs@exprs[,params]))>th){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
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
  
  observeEvent(input$fcs_input,{
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 0)
    paths <- input$fcs_input$datapath
    names <- as.vector(input$fcs_input$name)
    flow.frames <- lapply(c(1:length(names)), function(i) {
      progress$set(message = paste0("Reading file ", i, "/", length(paths), "."), value=i/length(paths))
      if(length(grep(".csv$",basename(paths[i])))>0){
        csv <- read.csv(paths[i], sep=input$sep, check.names = FALSE,header=TRUE,stringsAsFactors = FALSE)#header = TRUE, stringsAsFactors=FALSE)
        return(catch.createFCSfromCSV(csv))
      } else {
        fcs <- read.FCS.CIPHE(paths[i])
        # fcs <- read.FCS(paths[i])
        if(is.null(fcs)){
          showNotification(ui="FCS Corrupted can't be open !!", type = "error")
        } else {
          fcs@exprs[which(is.na(fcs@exprs))] <- 0
        }
        return(fcs)
      }
    })
    if(is.null(flow.frames[[1]])){
      values$flow.frames <- NULL
      showNotification(ui="FCS/CSV Corrupted or bad format and can't be open !!", type = "error")
      showNotification(ui="Error during upload all values reset", type="error")
    } else {
      names <- str_replace(names,".csv",".fcs")
      values$flow.frames <- flow.frames
      names(values$flow.frames) <- names
      values$names.files <- names
      # values$check.params <- (length(unique(unlist(lapply(flow.frames, function(x){dim(x)[2]})))) == 1)
      getParams()
    }
    progress$close()	
  })
  
  output$listFR <- renderUI({
    selectInput("idFR","Select FlowRepository ID",choices=flowRep.ls(),multiple = F)
  })
  
  output$uiNbrFilesFR <- renderUI({
    if(is.null(input$idFR)) return(NULL)
    print(input$idFR)
    print(paste0(length(flowRep.get(input$idFR)@fcs.files)," Files"))
  })
  
  observeEvent(input$loadFR,{
    progress <- Progress$new()
    progress$set(message="Load FlowRepository..",value=0.9)
    values$flow.frames <- FlowRepositoryReadFCS(flowRep.get(input$idFR))
    values$names.files <- names(values$flow.frames)
    getParams()
    progress$close()
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(values$flow.frames))
  })
  
  output$log <- renderText({
    return(values$log)
  })
  
  observe({
    if(is.null(values$flow.frames)) return(NULL)
    output$myTable <- DT::renderDataTable({
      if(is.null(values$flow.frames)) return(NULL)
      id <- names(values$flow.frames)
      temp  <- unlist(lapply(id, function(i){
        # temp <- paste0(temp, "")
        return(paste0("<input type=\"text\" id=\"",i,"\" value=\"",i,"\" size=\"50\">"))
      }))
      Events <- unlist(lapply(values$flow.frames,function(i){return(dim(i)[1])}))
      Parameters <- unlist(lapply(values$flow.frames,function(i){return(dim(i)[2])}))
      # Tube <- unlist(lapply(values$flow.frames,function(i){return(i@description[["WELL ID"]])}))
      table <- data.frame(
        Names = names(values$flow.frames),
        Rename = temp,
        # Tube = Tube,
        Events = Events,
        Parameters = Parameters
      )
      row.names(table) <- c(1:length(values$flow.frames))
      DT::datatable(table,options = list(orderClasses = TRUE,
        lengthMenu = FALSE,pageLength = 100,searching = FALSE,
        drawCallback= JS('function(settings) {Shiny.bindAll(this.api().table().node());}')
      ),selection='none',escape=F)
    })
  })
  
  observeEvent(input$renameFCS,{
    if(is.null(values$flow.frames)) return(NULL)
    id <- names(values$flow.frames)
    new.names <- unlist(lapply(id, function(i){
      return(paste0(input[[i]],".fcs"))
    }))
    names(values$flow.frames) <- new.names
    values$names.files <- new.names
    # editNames()
  })
  
  observeEvent(input$fcsToCSV,{
    filenames <- paste("data_", Sys.Date(), "_fcsToCSV.zip", sep="")
    progress <- Progress$new()
    fcs <- values$flow.frames
    names.fcs <- names(fcs)
    base <- getwd()
    tmpdir <- tempdir()
    setwd(tempdir())
    csv <- c()
    progress$set(message="Zip files", value=0.5)
    for(i in c(1:length(fcs))){
      cs <- fcs[[i]]@exprs
      path <- paste0(names.fcs[[i]],".csv")
      csv <- c(csv, path)
      write.csv(cs,path)
    }

    zip(zipfile="output.zip", files=csv)
    setwd(base)
    
    path.big.file <- paste0(tmpdir,"/","output.zip")
    path.big.file <- str_replace(path.big.file,"/media/data/html/","http://10.71.1.22/")
    progress$close()
    
    output$lnkBigDDlCSV <- renderText({
      return(paste0("<a href='",path.big.file,"'>Download</a>"))
    })
    
  })
  
##### PREPROCESS DATA ############################
  
  output$transMarkerOutput <- renderUI({
    selectInput("trans_marker","Transform marker",choice=values$params,multiple=TRUE)
  })
  
  observeEvent(input$addAllTransParams,{
    updateSelectInput(session, "trans_marker", selected = values$params)
  })
  
  observeEvent(input$divide,{
    values$flow.frames <- lapply(values$flow.frames, function(i){
      raw <- i@exprs
      if(!is.null(input$trans_marker) || length(input$trans_marker)>1){
        mat <- i@exprs[,input$trans_marker]/input$mult_value 
        raw[,input$trans_marker] <- mat
      } else {
        raw <- raw/input$mult_value
      }
      i@exprs <- raw
      values$log <- rbind(values$log,paste0("Divide by:",input$mult_value))
      return(i)
    })
  })
  
  observeEvent(input$multiple,{
    values$flow.frames <- lapply(values$flow.frames, function(i){
      raw <- i@exprs
      if(!is.null(input$trans_marker) || length(input$trans_marker)>1){
        mat <- i@exprs[,input$trans_marker]*input$mult_value 
        raw[,input$trans_marker] <- mat
      } else {
        raw <- raw*input$mult_value
      }
      i@exprs <- raw
      values$log <- rbind(values$log,paste0("Multiple by:",input$mult_value))
      return(i)
    })
  })
  
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  observeEvent(input$transform,{
    progress <- Progress$new()
    progress$set(message = "Transform...", value = 0)
    if(length(values$flow.frames)<1){return(NULL)}
    if(is.null(input$trans_meth)){return(NULL)}
    if(input$trans_meth == "flowVS"){
      res <- flowVS.CIPHE(values$flow.frames, markers=input$trans_marker)
      values$flow.frames <- res[[1]]
    } else {
      values$flow.frames <- lapply(c(1:length(values$flow.frames)), function(i){
        progress$set(message = "Transform...", value = i/length(values$flow.frames))
        if(input$trans_meth == "logicle"){
          fcs <- logiclTransformCIPHE(values$flow.frames[[i]], value=input$args_trans, marker=input$trans_marker)
        } 
        if(input$trans_meth == "arcsinh") {
          fcs <- arcsinhTransCIPHE(values$flow.frames[[i]], arg=input$args_trans, marker=input$trans_marker)
        }
        return(fcs)
      })
    }
    names(values$flow.frames) <- values$names.files
    values$log <- rbind(values$log,paste0("Transform:",input$trans_meth))
    progress$close()
  })
  
  observeEvent(input$comp,{
    progress <- Progress$new()
    progress$set(message = "Compensate...", value = 1)
    res <- lapply(c(1:length(values$flow.frames)), function(i){ if(i){
      progress$set(message = "Compensate...", value = i/length(values$flow.frames))
      if(is.null(values$flow.frames[[i]]@description[["SPILL"]])){
        showNotification("nocomp",ui = "No comp Matrice",type = "error",duration = 10)
        return(NULL)
      }
      return(compensate(values$flow.frames[[i]],values$flow.frames[[i]]@description[["SPILL"]]))
    }})
    if(!is.null(res[[1]])){
      names(values$flow.frames) <- values$names.files
      values$log <- rbind(values$log,paste0("Compensate"))
    }
    values$flow.frames <- res
    names(values$flow.frames) <- values$names.files
    progress$close()
  })
  
  observeEvent(input$decomp,{
    progress <- Progress$new()
    progress$set(message = "Uncomepsante...", value = 1)
    values$flow.frames <- lapply(c(1:length(values$flow.frames)), function(i){
      return(deCompensateFlowFrame(values$flow.frames[[i]],values$flow.frames[[i]]@description[["SPILL"]]))
    })
    names(values$flow.frames) <- values$names.files
    values$log <- rbind(values$log,paste0("DeCompensate"))
    progress$close()
  })
  
  observeEvent(input$detransform,{
    progress <- Progress$new()
    progress$set(message = "Detransform...", value = 0)
    if(length(values$flow.frames)<1){return(NULL)}
    if(is.null(input$transform)){return(NULL)}
    values$flow.frames <- lapply(c(1:length(values$flow.frames)), function(i){
      if(input$trans_meth == "logicle"){
        values$flow.frames[[i]] <- inversLogiclTransformCIPHE(values$flow.frames[[i]], value=input$args_trans, marker=input$trans_marker)
      } else {
        values$flow.frames[[i]] <- inversArcsinhTransCIPHE(values$flow.frames[[i]], arg=input$args_trans, marker=input$trans_marker)
      }
    })
    names(values$flow.frames) <- values$names.files
    values$log <- rbind(values$log,paste0("DeTransform:",input$trans_meth))
    progress$close()
  })
  
  output$selectFlowViz <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("flowVizParams","Select Markers",choice=values$params)
  })
  
  observeEvent(input$flowVizParams,{
    if(is.null(values$flow.frames))return(NULL)
    updateNumericInput(session,"jpxmin",value=min(values$flow.frames[[1]]@exprs[,input$flowVizParams]))
    updateNumericInput(session,"jpxmax",value=max(values$flow.frames[[1]]@exprs[,input$flowVizParams]))
  })
  
  observe({
    output$textPreview <- renderDataTable({
      if(input$tabs != "preprocess"){return(NULL)}
      if(is.null(values$flow.frames)){return(NULL)}
      return(data.frame(
            round(head(values$flow.frames[[sample(1:length(values$flow.frames),1)]]@exprs,10),2), 
          check.names = F,stringsAsFactors = F))
    },options = list(searching = FALSE,paging = FALSE))
      
    output$overView <- renderPlot({
      if(input$tabs != "preprocess"){return(NULL)}
      if(is.null(values$flow.frames)){return(NULL)}
      if(input$flowVizMod == "FlowViz"){
        readFlowSet <- function(flow.frames){
          out <- tryCatch({
            flowSet(flow.frames)
          },
          error = function(cond){
            return(NULL)
          })
        }
        fs <- readFlowSet(values$flow.frames)
        if(is.null(fs)){
          showNotification("Not Same parameters name",type="error")
          return(NULL)
        }
        if(is.null(input$flowVizParams)) return(NULL)
        return(densityplot(as.formula(paste0("~`",input$flowVizParams,"`")),fs,
                           xlim=c(input$jpxmin,input$jpxmax)))
        } else if(input$flowVizMod=="ggridge") {
          res <- lapply(values$flow.frames, function(i){return(i@exprs[,input$flowVizParams])})
          data <- data.frame(melt(res))
          colnames(data) <- c("value","Files")
          p <- ggplot(data,aes(x=value, y=Files, fill=Files))+ geom_density_ridges2() + 
            xlim(input$jpxmin,input$jpxmax) + ggtitle(input$flowVizParams) + xlab("Intensity") +
            theme(axis.text.y=element_blank(),axis.title.y = element_blank())
          return(p)
        }
    },height=input$height,width=input$width)
  })
  
  observeEvent(input$center,{
    if(is.null(input$trans_marker)){
      showNotification(ui="Select Markers transform",type="error")
      return(NULL)
    }
    if(is.null(values$flow.frames))return(NULL)
    progress <- Progress$new()
    progress$set(message="centering...",value=1)
    flow.frames <- values$flow.frames
    list.means <- c()
    flow.frames <- lapply(flow.frames, function(i){
      fcs <- i
      l.m <- c()
      mat <- apply(fcs@exprs[,input$trans_marker],2,function(j){
        l.m <<- c(l.m, mean(j))
        return(j-mean(j))
        # j <- scale(j,center = T,scale = F)
      })
      fcs@exprs[,input$trans_marker] <- mat
      names(l.m) <- colnames(fcs)
      list.means <<- c(list.means, list(l.m))
      return(fcs)
    })
    names(list.means) <- values$names.files
    names(flow.frames) <- values$names.files
    values$flow.frames <- flow.frames
    values$means.list <- list.means
    progress$close()
  })
  
  observeEvent(input$scale,{
    if(is.null(input$trans_marker)){
      showNotification(ui="Select Markers transform",type="error")
      return(NULL)
    }
    progress <- Progress$new()
    progress$set(message="scaling...",value=1)
    flow.frames <- values$flow.frames
    list.sd <- c()
    flow.frames <- lapply(flow.frames, function(i){
      fcs <- i
      l.s <- c()
      mat <- apply(fcs@exprs[,input$trans_marker],2,function(j){
        l.s <<- c(l.s, sd(j))
        return(j/sd(j))
      })
      fcs@exprs[,input$trans_marker] <- mat
      names(l.s) <- colnames(fcs)
      list.sd <<- c(list.sd, list(l.s))
      return(fcs)
    })
    names(list.sd) <- values$names.files
    names(flow.frames) <- values$names.files
    values$flow.frames <- flow.frames
    values$sd.list <- list.sd
    progress$close()
  })
  
  observeEvent(input$norm,{
    if(is.null(input$trans_marker)){
      showNotification(ui="Select Markers transform",type="error")
      return(NULL)
    }
    progress <- Progress$new()
    progress$set(message="normalize...",value=0.9)
    flow.frames <- values$flow.frames
    flow.frames <- lapply(flow.frames, function(fcs){
      fcs <-  normalize.FCS.CIPHE(fcs, markers=input$trans_marker, quantile=(input$quantConst/100))
      return(fcs)
    })
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    values$log <- rbind(values$log,paste0("Normalize"))
    progress$close()
  })

  observeEvent(input$clr_fn,{
    if(is.null(input$trans_marker)){
      showNotification(ui="Select Markers transform",type="error")
      return(NULL)
    }
    progress <- Progress$new()
    progress$set(message="normalize...",value=0.9)
    flow.frames <- values$flow.frames
    flow.frames <- lapply(flow.frames, function(fcs){
      fcs@exprs[,input$trans_marker] <- CLR_FN(data.matrix(fcs@exprs[,input$trans_marker]))
      return(fcs)
    })
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    values$log <- rbind(values$log,paste0("Normalize CLR"))
    progress$close()
  })
    
##### VIEW DATA ##################################
  
  output$selectFilePreview <- renderUI({
    if(length(values$names.files)<1){return(NULL)}
    selectInput("idpreview","Select Files",choices = values$names.files)
  })
  
  output$selectXPreview <- renderUI({
    if(length(values$params)<1){return(NULL)}
    selectInput("xpreview","X Parameters",choices=values$params,selected=values$params[1])
  })
  
  output$selectYPreview <- renderUI({
    if(length(values$params)<1){return(NULL)}
    selectInput("ypreview","Y Parameters",choices=values$params,selected=values$params[2])
  })
  
  output$selectParPreview <- renderUI({
    if(length(values$params)<1){return(NULL)}
    selectInput("ppreview","Parameters", choices=values$params, multiple = T,selected=values$params[1])
  })
  
  observeEvent(input$xpreview, {
    if(length(values$flow.frames)<1) return(NULL)
    output$xlimOutput <- renderUI({
      min <- round(min(values$flow.frames[[input$idpreview]]@exprs[,input$xpreview]),2)-1
      max <- round(max(values$flow.frames[[input$idpreview]]@exprs[,input$xpreview]),2)+1
      sliderInput("xlim","X limit",min=min, max=max, step=0.01,value=c(min,max))
    })
  })
  
  observeEvent(input$ypreview, {
    if(length(values$flow.frames)<1) return(NULL)
    output$ylimOutput <- renderUI({
      min <- round(min(values$flow.frames[[input$idpreview]]@exprs[,input$ypreview]),1)-1
      max <- round(max(values$flow.frames[[input$idpreview]]@exprs[,input$ypreview]),1)+1
      sliderInput("ylim","Y limit",min=min, max=max, step=0.01,value=c(min,max))
    })
  })

  output$selectHeatmapDim <- renderUI({
    selectInput("selectHeatmapDim","Select Heatmap Dim",choices=values$params, multiple=TRUE)
  })
  
  observe({
    if(is.null(values$flow.frames)){return(NULL)}
    height <- 700 ; width <- 700;

    # if(input$view_mod=="Densityplot"){
    #   if(!is.null(input$ppreview)){
    #     height <- 80*length(input$ppreview)
    #   }
    # }
    # if(input$number=="Ones"){height <- 800 ; width <- 800;}
    
    output$plotPreview <- renderPlot({
      if(length(values$flow.frames)<1) return(NULL)
      if(length(values$params)<1){return(NULL)}
      # if(is.null(input$ypreview) || is.null(input$xpreview)) return(NULL)

      if(!is.null(input$annotID)){ #select sub flowframe filter by annot
        id <- c()
        for(i in input$annotID){
          id.t <- which(values$flow.frames[[input$idpreview]]@exprs[,input$annotParams]==i)
          id <- c(id,id.t)
        }
        fcs <- values$flow.frames[[input$idpreview]][id,]
      } else {
        fcs <- values$flow.frames[[input$idpreview]]
      }

      fcs <- fcs[sample(1:dim(fcs)[1],round(input$samplingView*dim(fcs)[1])),]

      if(input$view_mod=="Scatterplot"){
        if(input$colorAxes == "Density"){
          plotDens(fcs,c(input$xpreview,input$ypreview),main=input$idpreview,xlim=input$xlim,ylim=input$ylim)
        } else {
          if(checkAnnot(fcs,input$colorAxes,100)==TRUE){
            palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
            colors <- palette(20)[as.numeric(cut(fcs@exprs[,c(input$colorAxes)],breaks=20))]
          } else {
            colors <- rainbow(length(unique(fcs@exprs[,c(input$colorAxes)]))+1)[fcs@exprs[,c(input$colorAxes)]+1]
          }
          plot(fcs@exprs[,c(input$xpreview,input$ypreview)],
            main=input$idpreview,xlim=input$xlim,cex=input$cex,
            ylim=input$ylim,pch=20,col=colors)
        }
      }

      if(input$view_mod=="Heatmap"){
        if(length(input$selectHeatmapDim)<2)return(NULL)
        my_palette <- colorRampPalette(c("green", "black", "red"))(n = 299)
        width <- 1200
        if(input$colorAxes == "Density" || checkAnnot(fcs,input$colorAxes,100)==TRUE){
          rowsidecolors <- as.vector(as.character(rep(1,dim(fcs)[1])))
        } else {
          rowsidecolors <- as.vector(as.character(fcs@exprs[,input$colorAxes]))
        }
        fcs <- fcs[,input$selectHeatmapDim]
        heatmap.2(as.matrix(t(fcs@exprs)),
          symkey=FALSE,
          symbreaks=FALSE,
          Rowv = NULL,
          density.info="none", trace="none", scale="row",dendrogram="none",
          col=my_palette,
          cexRow=1.5, cexCol=1.5,
          margins=c(10,25),
          keysize = 0.5,
          key.par=list(mar=c(1.5,1.5,1,1),cex.axis=0.5),
          key.title="",
          key.xlab = "none",
          breaks = seq(-1, 1, length.out = 300),
          ColSideColors = rowsidecolors,
          labCol = FALSE
        )
      }

    },height = height, width = width)
    
  })
  
  output$selectFilePreview2 <- renderUI({
    if(length(values$names.files)<1){return(NULL)}
    selectInput("idpreview2","Select Files",choices = values$names.files)
  })
  
  output$selectFileMultiple <- renderUI({
    if(length(values$names.files)<1){return(NULL)}
    selectInput("idpreviewM","Select Files",choices = values$names.files,multiple = TRUE)
  })
  
  output$colorMods <- renderUI({
    if(is.null(values$params)){return(NULL)}
    colors <- c("Density",values$params)
    selectInput("colorAxes","Select Cluster", choices = colors, multiple = F)
  })
  
  output$uiAnnotParams <- renderUI({
    if(is.null(values$params)){return(NULL)}
    selectInput("annotParams","Select subframe",choices=values$params)
  })
  
  output$uiAnnotID <- renderUI({
    if(is.null(input$annotParams)) return(NULL)
    if(checkAnnot(values$flow.frames[[input$idpreview]],params = input$annotParams,th=100)){
      choice <- NULL 
    } else {
      choice <- unique(values$flow.frames[[input$idpreview]]@exprs[,input$annotParams])
    }
    selectInput("annotID"," = ",choices=choice,multiple=TRUE)
  })
  
##### CONFIGURATION ##############################

  observe({
    if(is.null(values$flow.frames)){return(NULL)}
    output$keywordsTable <- renderRHandsontable({
      data <- lapply(values$flow.frames, function(i){return(i@description)})
      
      data <- lapply(data, function(j){
        if(length(grep("SPILL",names(j)))>0){
          return(j[-grep("SPILL",names(j))])
        } else {
          return(j)
        }
      })
      
      mat <- do.call(cbind,data)
      mat <- as.data.frame(apply(data.frame(mat), 2, as.character),check.names=F,stringsAsFactors = FALSE)

      if(dim(mat)[1] != length(names(values$flow.frames[[1]]@description))){
        row.names(mat) <- names(values$flow.frames[[1]]@description)[-grep("SPILL",names(values$flow.frames[[1]]@description))]
      } else {
        row.names(mat) <- names(values$flow.frames[[1]]@description)
      }

      # values$keywords <- row.names(mat)
      rhandsontable(mat, readOnly = FALSE,height=800,rowHeaderWidth = 200) %>%
        hot_context_menu(
          customOpts = list(
            csv = list(name = "Download to CSV",
                       callback = htmlwidgets::JS(
                         "function (key, options) {
                           var csv = csvString(this, sep=',', dec='.');

                           var link = document.createElement('a');
                           link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                             encodeURIComponent(csv));
                           link.setAttribute('download', 'data.csv');

                           document.body.appendChild(link);
                           link.click();
                           document.body.removeChild(link);
         }"))))
    })
  })
  
  observeEvent(input$applyKeywords,{
    progress <- Progress$new()
    progress$set("Change keywords...",value=0.9)
    mat <- do.call(rbind,input$keywordsTable$data)
    flow.frames <- values$flow.frames
    row.names(mat) <- values$keywords
    for(i in c(1:length(flow.frames))){
      for(j in names(flow.frames[[i]]@description)){
        if(length(grep(j,values$keywords))>0){
          flow.frames[[i]]@description[[j]] <- unlist(mat[,i][[j]])
        }
      }
    }
    progress$close()
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    showNotification("Keywords Changed",type="message")
    values$log <- rbind(values$log,paste0("Keywords Edit"))
  })
  
  observeEvent(input$template,{
    paths <- input$template$datapath
    names <- as.vector(input$template$name)
    if(length(grep(".fcs$",names))>0){print("fcs")}
    if(length(grep(".csv$",names))>0){print("csv")}
  })

  observe({
    if(length(unique(unlist(lapply(values$flow.frames, function(x){dim(x)[2]})))) == 1) return(NULL)
    if(length(values$flow.frames)<1){return(NULL)}
    showNotification("Vos fichiers n'ont pas le m??me nombre de param??tres !")
  })

  observe({ ## View and edit table
    if(length(unique(unlist(lapply(values$flow.frames, function(x){dim(x)[2]})))) != 1) return(NULL)
    if(is.null(input$selectViewLabels)) return(NULL)
    if(input$tabs != "conf" || is.null(input$tabs)) {return(NULL)}
    if(!is.null(values$flow.frames)){
      if(length(values$flow.frames) == 1) {
        flow.frame <- values$flow.frames[[1]]
        mat <- as.data.frame(pData(flow.frame@parameters[,c(1,2)]),stringsAsFactors = F)
        colnames(mat) <- c("names","labels")
        output$paramsTable <- renderRHandsontable({rhandsontable(mat, selectCallback = TRUE)})
    }
    else {
        res <- lapply(values$flow.frames, function(x){return(pData(x@parameters)[,c(1,2)])})
        mat <- do.call(cbind, res)
        cols1 <- paste0(names(values$flow.frames),"_names")
        cols2 <- paste0(names(values$flow.frames),"_labels")
        colnames(mat) <- as.vector(rbind(cols1,cols2))
        values$params.table <- mat
        output$paramsTable <- renderRHandsontable({
          colnames <- colnames(values$params.table)[grep(input$selectViewLabels,colnames(values$params.table))]
          mat <- as.data.frame(values$params.table[,grep(input$selectViewLabels,colnames(values$params.table))], col.names=colnames,stringsAsFactors = F)
          colnames(mat) <- colnames
          rhandsontable(mat, readOnly = F)
        })
      }
    } else {
      output$paramsTable <- renderRHandsontable({
        rhandsontable(data.frame(NULL,stringsAsFactors = F,check.names = FALSE,readOnly = F)) %>%  
          hot_cell(
            
          )
      })
    }
  })
  
  observeEvent(input$addPattern,{
    x <- input$paramsTable_select$select$r
    y <- input$paramsTable_select$select$c
    print(input$paramsTable_select$data[[x]][[y]])
  })

  observeEvent(input$saveParams,{
    if(length(unique(unlist(lapply(values$flow.frames, function(x){dim(x)[2]})))) != 1) return(NULL)
    if(length(values$flow.frames)<1){return(NULL)}
    mat <- do.call(rbind,input$paramsTable$data)
    
    if(!is.null(values$flow.frames)){
      if(input$selectViewLabels == "labels"){
        flow.frames <- values$flow.frames
        progress <- Progress$new()
        progress$set(message = "rename labels", value = 0.5)
        flow.frames <- lapply(c(1:length(flow.frames)),function(x){
          fcs <- flow.frames[[x]]
          temp <- mat[,x]
          temp[unlist(lapply(temp , is.null))] <- NA 
          pData(fcs@parameters)[,"desc"] <- unlist(temp)
          for(i in c(1:length(colnames(fcs)))){
            fcs@description[[paste0("$P",i,"S")]] <- unlist(mat[i,x])
          }
          return(fcs)
        })
        names(flow.frames) <- names(values$flow.frames)
        values$flow.frames <- flow.frames
        progress$close()
      } else if (input$selectViewLabels == "names"){
        flow.frames <- values$flow.frames
        progress <- Progress$new()
        progress$set(message = "rename names", value = 0.5)
        flow.frames <- lapply(c(1:length(flow.frames)), function(x){
          fcs <- flow.frames[[x]]
          temp <- mat[,x]
          temp[unlist(lapply(temp , is.null))] <- NA
          pData(fcs@parameters)[,"name"] <- unlist(temp)
          for(i in c(1:length(colnames(fcs)))){
            fcs@description[[paste0("$P",i,"N")]] <- unlist(mat[i,x])
            colnames(fcs)[i] <- unlist(mat[i,x])
          }
          return(fcs)
        })
        names(flow.frames) <- names(values$flow.frames)
        values$flow.frames <- flow.frames
        progress$close()
      }
    } else if (length(values$flow.frames)==1){
      progress <- Progress$new()
      progress$set(message = "rename labels & names", value = 0.5)
      fcs <- values$flow.frames[[1]]
      temp <- mat[,2]
      temp[unlist(lapply(temp , is.null))] <- NA 
      pData(fcs@parameters)[,"desc"] <- unlist(temp)
      for(i in c(1:length(colnames(fcs)))){
        fcs@description[[paste0("$P",i,"S")]] <- unlist(mat[i,2])
      }
      temp <- mat[,1]
      temp[unlist(lapply(temp , is.null))] <- NA 
      pData(fcs@parameters)[,"name"] <- unlist(temp)
      for(i in c(1:length(colnames(fcs)))){
        fcs@description[[paste0("$P",i,"N")]] <-unlist(mat[i,1])
        colnames(fcs)[i] <- unlist(mat[i,1])
      }
      flow.frames <- list(fcs)
      names(flow.frames) <- names(values$flow.frames)
      values$flow.frames <- flow.frames
      progress$close()
    }
    getParams()
    showNotification(ui="Label Edit Done", type="message")
    values$log <- rbind(values$log,paste0("Labelled"))
  })
  
  output$uiSelectEditAnnot <- renderUI({
    if(is.null(values$params))return(NULL)
    selectInput("selectEditAnnot","Select Params",choices = values$params,multiple=F)
  })
  
  output$deleteRowExprs <- renderRHandsontable({
    if(is.null(input$selectEditAnnot))return(NULL)
    if(checkAnnot(values$flow.frames[[1]],input$selectEditAnnot,201)==TRUE)return(NULL)
    cols2 <- c(1:max(unique(values$flow.frames[[1]]@exprs[,input$selectEditAnnot])))
    colsp <- c()
    lapply(values$flow.frames,function(fcs){
      temp <- unlist(lapply(cols2,function(j){
        return(length(which(fcs@exprs[,input$selectEditAnnot]==j)))
      }))
      colsp <<- c(colsp,list(temp))
    })
    cols1 <- rep(TRUE,length(cols2))
    table <- data.frame(do.call(cbind,colsp))
    table <- cbind(rep(FALSE,length(cols2)),cols2,table)
    colnames(table) <- c("delete","value",values$names.files)
    hot <- rhandsontable(table,height = 800)
    hot <- hot_cols(hot, fixedColumnsLeft = 3, renderer = "
            function(instance, td, row, col, prop, value, cellProperties) {
                if(col == 0)
                    Handsontable.renderers.CheckboxRenderer.apply(this, arguments)
                else {
                    Handsontable.renderers.TextRenderer.apply(this, arguments)
                    if(instance.params != null) {
                        if(instance.params.data[row][0]){
                            td.style.background = 'lightgrey'
                        }
                    }
                }
                return(td)
            }"
      )
    return(hot)
  })
  
  observeEvent(input$deleteRow,{
    if(is.null(input$deleteRowExprs))return(NULL)
    DF <- data.frame(hot_to_r(input$deleteRowExprs),stringsAsFactors = FALSE)
    id <- which(DF[,1]==TRUE)
    if(length(id)==0)return(NULL)
    flow.frames <- values$flow.frames
    flow.frames <- lapply(flow.frames, function(i){
      fcs <- i
      id.delete <- as.vector(unlist(lapply(id,function(j){
        return(which(fcs@exprs[,input$selectEditAnnot]==j))
      })))
      if(length(id.delete)>0){
        fcs <- fcs[-id.delete,]
      }
      return(fcs)
    })
    showNotification(ui="Your rows has been delete", type = "message")
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    getParams()
  })
  
##### DESCRIPTIVE STATS ##########################
  
  getMode <- function(v){
    uniqv <- unique(v)
    return(uniqv[which.max(tabulate(match(v, uniqv)))])
  }
   
  output$mfiMarkerUI <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("mfiMark","Select MFi Marker's",choices=c(values$params),multiple=TRUE)	
  })
  
  observeEvent(input$computeMFI,{
    if(length(input$mfiMark)<1)return(NULL)
    stats <- c(input$stats,c(input$quantiles))
    print(stats)
    res <- list()
    for(i in input$mfiMark){
      l1 <- list()
      for(j in stats){
        l <- unlist(lapply(values$flow.frames,function(k){
          if(j == "Mean"){return(round(mean(k@exprs[,i]),2))}
          if(j == "Median"){return(round(median(k@exprs[,i]),2))}
          if(j == "Mode"){return(round(getMode(k@exprs[,i]),2))}
          if(j == "SD"){return(round(sd(k@exprs[,i]),2))}
          if(j == "Min"){return(round(min(k@exprs[,i]),2))}
          if(j == "Max"){return(round(max(k@exprs[,i]),2))}
          if(is.numeric(as.numeric(j))){
            return(round(quantile(k@exprs[,i],probs = as.numeric(j)/100)))
          }
        }))
        l1 <- c(l1,list(l))
      }
      res <- c(res, list(l1))
    }
    
    res2 <- lapply(res, function(i){t(do.call(rbind,i))})
    table <- do.call(cbind, res2)
    row.names(table) <- values$names.files
    colnames(table) <-as.vector(outer(stats, input$mfiMark, paste, sep="#"))
    values$mfi.table <- table
  })
  
  output$mfiTable <- renderDataTable({
    if(is.null(values$mfi.table)) return(NULL)
    return(values$mfi.table)
  },rownames = TRUE)
  
  output$ddlMFItable <- downloadHandler(
    filename = function(){
      return(paste0("Stats.csv"))
    },
    content = function(file){
      write.csv(values$mfi.table, file)
    }
  )
  
  output$uiParamsSelect <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    if(is.null(values$params))return(NULL)
    selectInput("ParamsSelect","Select Cluster ID",choices=values$params,multiple=FALSE)
  })
  
  output$uiMarkerSelect <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    if(is.null(values$params))return(NULL)
    if(checkAnnot(fcs=values$flow.frames[[1]],input$ParamsSelect,th = 300)==TRUE)return(NULL)
    selectInput("markerSelect","Params Export",choices=values$params,multiple=TRUE)
  })
  
  observeEvent(input$addAllCreateFCS,{
    updateSelectInput(session,"markerSelect",selected=values$params)
  })
  
  observeEvent(input$clearAllCreateFCS,{
    updateSelectInput(session,"markerSelect",selected=1)
  })
  
  observeEvent(input$computeCreate,{
    if(is.null(input$markerSelect))return(NULL)
    stats <- c(input$statsCreate)
    print(stats)
    progress <- Progress$new()
    progress$set(message="Compute progress", value=0.9)
    fcs <- c()
    for(x in values$flow.frames){
      cl <- c()
      for(c in c(1:length(unique(x@exprs[,input$ParamsSelect])))){
        k <- x[which(x@exprs[,input$ParamsSelect]==c),]
        m <- c()
        for(y in input$markerSelect){
          v <- c()
          for(z in input$statsCreate){
            a <- NULL
            if(z == "Mean"){a <- round(mean(k@exprs[,y]),2)}
            if(z == "Median"){a <- round(median(k@exprs[,y]),2)}
            if(z == "Mode"){a <- round(getMode(k@exprs[,y]),2)}
            if(z == "SD"){a <- round(sd(k@exprs[,y]),2)}
            if(z == "Min"){a <- round(min(k@exprs[,y]),2)}
            if(z == "Max"){a <- round(max(k@exprs[,y]),2)}
            v <- c(v, a)
          }
          m <- c(m, v)
        }
        if(length(grep("size",input$statsCreate))>0){
          m <- c(m, dim(k)[1])
        }
        cl <- c(cl,list(m))
      }
      mat <- do.call(rbind,cl)
      fcs <- c(fcs, list(mat))
    }
    fcs <- lapply(fcs, function(i){
      row <- unique(x@exprs[,input$ParamsSelect])
      d <- input$statsCreate[-grep("size",input$statsCreate)]
      p <- names(values$params)[which(values$params%in%input$markerSelect)]
      colnames(i)[c(1:((dim(i)[2])-1))] <- as.vector(outer(d,p,paste,sep="_"))
      i <- cbind(row,i)
      colnames(i)[1] <- input$ParamsSelect
      return(i)
    })
    values$create.table <- fcs
    names(values$create.table) <- values$names.files
    progress$close()
  })
  
  output$previewTableCreateFCS <- renderRHandsontable({
    if(is.null(values$create.table))return(NULL)
    return(rhandsontable(values$create.table[[1]]))
  })
  
  output$ddlCSVCreateFromClsuter <- downloadHandler(
    filename = function(){
      paste0("output.zip")
    },
    content = function(file){
      csv <- values$create.table
      root <- getwd()
      tmpdir <- tempdir()
      setwd(tempdir())
      fs <- c()
      progress$set(message="Zip files", value=0.5)
      for(i in c(1:length(csv))){
        name <- names(values$flow.frames)[i]
        fcs <- csv[[i]]
        path <- gsub(".fcs$",".csv",name)
        fs <- c(fs, path)
        write.csv(fcs, path)
      }
      zip(zipfile=file, files=fs)
      setwd(root)
      progress$close()
    },
    contentType = "application/zip"
  )
  
##### DIFERENTILE STATS ##########################
  
  observe({
    if(is.null(values$params)) return(NULL)
    if(input$tabs != "difstats")return(NULL)
    updateSelectInput(session,"annotGroup",choices=values$params)
  })
  
  observeEvent(input$annotGroup,{
    output$uiGroup1 <- renderUI({
      if(checkAnnot(values$flow.frames[[1]],input$annotGroup,50)==TRUE){return(NULL)}
      selectInput("group1","Group One",choice=unique(values$flow.frames[[1]]@exprs[,input$annotGroup]),multiple = F)
    })
    output$uiGroup2 <- renderUI({
      if(checkAnnot(values$flow.frames[[1]],input$annotGroup,50)==TRUE){return(NULL)}
      selectInput("group2","Group One",choice=unique(values$flow.frames[[1]]@exprs[,input$annotGroup]),multiple = F)
    })
    output$uiVariable <- renderUI({
      selectInput("variable","Variable",choice=values$params, multiple=T)
    })
    output$uiAnnotPoint <- renderUI({
      if(checkAnnot(values$flow.frames[[1]],input$annotGroup,50)==TRUE){return(NULL)}
      selectInput("annotPoint","Annot",choice=values$params, multiple=F)
    })
  })
  
  observeEvent(input$addAllVolcano,{
    updateSelectInput(session, "variable", choices=values$params,selected=values$params)
  })
  
  observeEvent(input$clearAllVolcano,{
    updateSelectInput(session, "variable", choices=values$params,selected=NULL)
  })
  
  observeEvent(input$addGrepVolcano,{
    if(is.null(input$grepVolcano))return(NULL)
    id <- grep(input$grepVolcano,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session, "variable", choices=values$params,selected=values$params[id])
  })
  
  # observeEvent(input$annotPoint,{
  #   output$uiSelectAnnotPoint <- renderUI({
  #     if(checkAnnot(values$flow.frames[[1]],input$annotPoint,50)==TRUE){return(NULL)}
  #     selectInput("selectAnnotPoint","=",choice=unique(values$flow.frames[[1]]@exprs[,input$annotPoint]),multiple = F)
  #   })
  # })
  
  observeEvent(input$plotVolcano,{
    if(is.null(input$variable)) return(NULL)
    if(is.null(input$annotPoint)) return(NULL)

    data <- values$flow.frames[[1]]@exprs
    g1 <- data[which(data[,input$annotGroup]==input$group1),]
    g2 <- data[which(data[,input$annotGroup]==input$group2),]
    if(checkAnnot(values$flow.frames[[1]],input$annotPoint,50)==TRUE){return(NULL)}
    
    markers <- c()
    logFC <- c()
    logPvalue <- c()
    groupe <- c()
    progress <- Progress$new()
    progress$set(message="Compute VolvanoPlot",value=0.9)
    for(i in unique(values$flow.frames[[1]]@exprs[,input$annotPoint])){
      for(j in input$variable){
        if(length(which(g1[,input$annotPoint]==i))>1 && length(which(g2[,input$annotPoint]==i))>1){
          a <- g1[which(g1[,input$annotPoint]==i),j]
          b <- g2[which(g2[,input$annotPoint]==i),j]
          markers <- c(markers, list(j))
          logFC <- c(logFC, list(log2(mean(a)/mean(b))))
          logPvalue <- c(logPvalue, list(-log10(t.test(x=a,y=b)$p.value)))
          groupe <- c(groupe, list(i))
        } else {
          logFC <- c(logFC, list(0))
          logPvalue <- c(logPvalue, list(0))
          markers <- c(markers, list(j))
          groupe <- c(groupe, list(i))
        }
      }
    }
    
    data <- data.frame(
      m = unlist(markers),
      log2FC = unlist(logFC),
      logPvalue = unlist(logPvalue),
      groupe = unlist(groupe)
    )
  
    progress$close()
    
    output$volcanoPlot <- renderPlotly({
      plot_ly(data, x = ~log2FC,y=~logPvalue,color=~as.factor(groupe),text=~m,
        mode='markers',marker=list(size=5))%>%layout(dragmode="select2")
    })
    
  })
  
##### MANAGE CONCATENATE #########################

  output$clusteringui3 <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    dd <- values$clustering.groups
    if(length(values$clustering.groups)==0){values$clustering.groups <- NULL}
    return(tagList(mapply(get_cluster_groups_table, dd, names(dd), SIMPLIFY = F)))
  })
  
  output$clusteringui2 <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("clusteringui_files_list",label = "File list",
                choices = names(values$flow.frames),selectize = F,
                multiple = T,width = "100%"
    )
  })
  
  observeEvent(input$clusteringui_add_all_groups, {
    if(is.null(input$clusteringui_add_all_groups)) return(NULL)
    for(name in names(values$flow.frames)){
      values$clustering.groups <- c(values$clustering.groups, setNames(list(name), name[1]))
    }
  })
  
  observe({
    if(input$tabs != "concat"){return(NULL)}
    key <- input$clusteringui_remove_clustering_group$key
    if(!is.null(key) && key != ""){ isolate({values$clustering.groups[key] <- NULL})}
  })
  
  observeEvent(input$clusteringui_add_clustering_group, {
    if(is.null(input$clusteringui_files_list)) {return(NULL)}
    files_list <- isolate({input$clusteringui_files_list})
    values$clustering.groups <- c(values$clustering.groups, setNames(list(files_list), files_list[1]))
  })
  
  output$fileGrouped <- reactive({
    if(is.null(values$clustering.groups)) return(FALSE)
    return(length(values$clustering.groups)>0)
  })
  
  outputOptions(output, "fileGrouped", suspendWhenHidden = FALSE)
  
  output$fileEnriched <- reactive({
    if(is.null(values$centres)) return(FALSE)
    return(length(values$centres)>0)
  })
  
  outputOptions(output, "fileEnriched", suspendWhenHidden = FALSE)
  
  observeEvent(input$addConcat,{
    if(is.null(input$clusteringui_files_list)) {return(NULL)}
    if(is.null(values$clustering.groups)) {return(NULL)}
    flow.frames <- c()
    for(i in values$clustering.groups){
      flow.frames <- c(flow.frames,concatenateCIPHE(values$flow.frames[c(i)],input$concatParams))
    }
    names(flow.frames) <- names(values$clustering.groups)
    values$names.files <- names(flow.frames)
    values$flow.frames <- flow.frames
    values$log <- rbind(values$log,paste0("Concatenate"))
    getParams()
  })
  
  output$selectSepMarkers <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    selectInput("sepMarkers","Select Split Markers",choices=values$params,multiple=FALSE)
  })
  
  output$theoricNbrEvents <- renderTable({
    if(is.null(input$sepMarkers)){return(NULL)}
    nbr.files <- unlist(lapply(values$flow.frames, function(i){length(unique(i@exprs[,input$sepMarkers]))}))
    size <- unlist(lapply(values$flow.frames,function(i){return(dim(i)[1])}))
    return(data.frame(files=names(values$flow.frames),size = size,nbr.files=nbr.files))
  })
  
  observeEvent(input$selectSep,{
    if(is.null(input$sepMarkers)){return(NULL)}
    if(length(unique(values$flow.frames[[1]]@exprs[,input$sepMarkers]))>100){return(NULL)}
    cols <- c(input$sepMarkers,names(values$flow.frames))
    nbr.files <-unlist(lapply(values$flow.frames, function(i){length(unique(i@exprs[,input$sepMarkers]))}))
    id.max <- max(nbr.files)
    print(id.max)
    flow.frames <- c()
    for(i in values$flow.frames){
      flow.frames.temp <- lapply(c(1:(id.max)),function(j){
        id <- which(i@exprs[,input$sepMarkers]==j)
        return(i[id,])
      })
      flow.frames <- c(flow.frames, flow.frames.temp)
    }
    values$flow.frames <- flow.frames
    values$names.files <- paste0(values$names.files,c(1:id.max),".fcs")
    names(values$flow.frames) <- values$names.files
    getParams()
    # editNames()
    values$log <- rbind(values$log,paste0("Unconcatenate"))
    # output$textNewfiles <- renderHandsonTable({
    #   return()
    # })
  })
  
##### SAMPLING RANDOM ############################
  
  output$sizeSampling <- renderDataTable({
    if(is.null(values$flow.frames)){return(NULL)}
    return(data.frame(filesname=names(values$flow.frames),size=unlist(lapply(values$flow.frames,function(j){return(dim(j)[1])}))))
  },options=list(paging=FALSE,searching=FALSE,pageLength=20))
  
  output$sliderInputByFiles <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    slider_output_list <- lapply(names(values$flow.frames), function(x){
        id_name <- paste0(x,"_freq")
        object <- sliderInput(id_name,paste0("Percentiles of ",x),min=1,max=100,value=100,step=1)
    })
    do.call(tagList, slider_output_list)
  })
  
  output$numberInputByFiles <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    number_output_list <- lapply(names(values$flow.frames), function(x){
      id_name <- paste0(x,"_evt")
      object <- numericInput(id_name,paste0("Number events of ",x),min=1,step=1,value=dim(values$flow.frames[[x]])[1])
    })
    do.call(tagList, number_output_list)
  })
  
  observeEvent(input$applyAllSampling,{
    progress <- Progress$new()
    progress$set(message="Random Downsampling All...",value=0.9)
    flow.frames <- values$flow.frames
    flow.frames <- lapply(flow.frames, function(j){
      if(input$sampleMode == "percent"){
        return(j[sample(c(1:dim(j)[1]),round((input$percentile/100)*dim(j)[1])),])
      } else {
        return(j[sample(c(1:dim(j)[1]),input$events),])
      }
    })
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    showNotification("sample",type = "message",ui = "Sampling Done")
    values$log <- rbind(values$log, paste0("Random Dowsampl"))
    progress$close()
  })
  
  observeEvent(input$applyIndSampling,{
    progress <- Progress$new()
    progress$set(message="Random Downsampling Individual...",value=0.9)
    flow.frames <- values$flow.frames
    flow.frames <- lapply(c(1:length(flow.frames)), function(j){
      f <- input[[paste0(names(values$flow.frames)[j],"_freq")]]
      e <- input[[paste0(names(values$flow.frames)[j],"_evt")]]

      fcs <- flow.frames[[j]]
      if(input$sampleMode == "percent"){
        return(fcs[sample(c(1:dim(fcs)[1]),round((f/100)*dim(fcs)[1])),])
      } else {
        return(fcs[sample(c(1:dim(fcs)[1]),e),])
      }
    })
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    progress$close()
    showNotification("sample",type = "message",ui = "Sampling Done")
    values$log <- rbind(values$log, paste0("Random Dowsampl"))
  })
  
##### SAMPLING PERCENTILE ########################
  
  output$selectXdens <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    selectInput("Xdens","X",choices = values$params,selected = values$params[[1]])
  })
  
  output$selectYdens <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    selectInput("Ydens","Y",choices = values$params,selected = values$params[[2]])
  })
  
  output$selectFilesDens <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    selectInput("IdDens","Select Files", choices=values$names.files)
  })
  
  observeEvent(input$activeDens,{
    if(is.null(input$Xdens) || is.null(input$Ydens)) return(NULL)
    if(is.null(input$probDensity1)) return(NULL)
    progress <- Progress$new()
    progress$set(message="Prepare density plot...", value=0.9)
    xydata <- values$flow.frames[[input$IdDens]]@exprs[,c(input$Xdens,input$Ydens)]
    fit.trim<- locfit( ~ lp(xydata[,1], xydata[,2], nn=input$probDensity1/100,scale=FALSE))
    dens <-(fitted(fit.trim))
    med.dens <- median(dens)
    values$density <- list(xydata,dens)
    progress$close()
  })
  
  output$uiSliderProdDensity2 <- renderUI({
    if(is.null(values$density)){return(NULL)}
    sliderInput("probDensity2","Sensibility of Density",min=0.1,max=50, sep=0.1,value=10)
  })
  
  observeEvent(input$probDensity2,{
    if(is.null(values$density))return(NULL)
    dens <- values$density[[2]]
    xydata <- values$density[[1]]
    high <- dens>quantile(dens,probs=input$probDensity2/100)
    output$densityPlotSample <- renderPlot({
      plot(xydata,pch=".",col="lightgrey",main=paste(input$probDensity2,"% excluded"))
      points(xydata[high,],pch=".",col="black")
    })
    down.sample <- values$flow.frames[[input$IdDens]][high,]
    output$densityTableSample <- renderTable({
      nbr <- c(dim(xydata)[1],dim(xydata[high,])[1],(dim(xydata)[1]-dim(xydata[high,])[1]))
      mat <- t(data.frame("All Events"=nbr[1],"Selected Events"=nbr[2],"Deleted Events"=nbr[3],check.names = F))
      return(mat)
    },rownames = T,colnames = F)
    output$histPlotX <- renderPlot({
      plot(hist(xydata[,1],breaks=dim(xydata)[1]/input$probDensity1,col="blue",border="blue")$density,type="l",col="red")
    })
    output$histPlotY <- renderPlot({
      plot(hist(xydata[,2],breaks=dim(xydata)[1]/input$probDensity1,col="blue",border="blue")$density,type="l",col="red")
    })
  })
  
  observeEvent(input$applyDensSample,{
    if(is.null(values$density)) return(NULL)
    progress <- Progress$new()
    progress$set(message = "Apply to your file...",value=0.9)
    dens <- values$density[[2]]
    xydata <- values$density[[1]]
    high <- dens>quantile(dens,probs=input$probDensity2/100)
    fcs <- values$flow.frames[[input$IdDens]][high,]
    values$flow.frames[[input$IdDens]] <- fcs
    names(values$flow.frames) <- values$names.files
    showNotification("Apply Complete",type="message")
    progress$close()
    values$density <- NULL
  })
  
  observeEvent(input$enrichDensSample,{
    if(is.null(values$density)) return(NULL)
    progress <- Progress$new()
    progress$set(message = "Apply to your file...",value=0.9)
    dens <- values$density[[2]]
    xydata <- values$density[[1]]
    high <- dens>quantile(dens,probs=input$probDensity2/100)
    new.column <- rep(0,dim(values$flow.frames[[input$IdDens]])[1])
    new.column[high] <- 1
    fcs <- enrich.FCS.CIPHE(values$flow.frames[[input$IdDens]],new.column,"Density.sample")
    values$flow.frames[[input$IdDens]] <- fcs
    names(values$flow.frames) <- values$names.files
    getParams()
    progress$close()
    showNotification("Enrich density done",type="message")
    values$log <- rbind(values$log, paste0("Density Dowsampl"))
  })
  
##### CLUSTERING #################################

  output$clusterParamsOutput <- renderUI({
    selectInput("clusterParams","Params for Clustering",choice=values$params, multiple=TRUE)
  })
  
  observe({
    if(input$tabs != "clusters") return(NULL)
    if(is.null(values$flow.frames)) return(NULL)
    updateSelectInput(session, "clusterParams", choices=values$params)
  })
  
  observeEvent(input$allParamsClust,{
    updateSelectInput(session, "clusterParams", choices=values$params,selected=values$params)
  })
  
  observeEvent(input$clearParamsClust,{
    updateSelectInput(session, "clusterParams", choices=values$params,selected=NULL)
  })
  
  observeEvent(input$addGrepParamsClust,{
    if(is.null(input$grepParamsClust))return(NULL)
    id <- grep(input$grepParamsClust,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session, "clusterParams", choices=values$params,selected=values$params[id])
  })
  
  observeEvent(input$runClustering,{
    if(is.null(input$clusterParams)){return(NULL)}
    progress <- Progress$new()
    flow.frames <- values$flow.frames
    progress$set(message="Clustering",value=0)
    if(input$clustering_meth == "HClust"){
      if(dim(flow.frames[[1]])[1]>15000){
        progress$close()
        showNotification("Use Hierarchical Clustering with less than 10 000 events please",type = "error")
        return(NULL)
      }
    }
    clustering <- list()
    # centres <- list()
    v <- 0
    for(i in flow.frames){
      v <- v +0.1
      progress$set(message="Clustering",value=v)
      data <- i@exprs[,input$clusterParams]
      if(input$clustering_meth == "k-means"){
        km <- kmeans(data, centers=input$ncluster, iter.max=200)
        clustering <- c(clustering,list(km$cluster))
        # centres <- c(centres,list(as.matrix(km$centers)))
      }
      if(input$clustering_meth == "RPhenograph"){
        res <- Rphenograph(data,k=input$knn)
        clusters <- membership(res[[2]])
        print("test")
        # print(clusters)
        clustering <- c(clustering, list(clusters))
        # centres <- c(centres,)
      }
      if(input$clustering_meth == "HClust"){
        hcl <- hclust(dist(data))
        id <- cutree(hcl, k=input$ncluster)
        clustering <- c(clustering, list(id))
      }
      if(input$clustering_meth == "CLARA"){
        clara <- clara(data, samples=50,k=input$ncluster)
        clustering <- c(clustering, list(clara[[4]]))
      }
      if(input$clustering_meth == "FlowSOM"){
        fcs <- flowFrame(data)
        fs <- FlowSOM::BuildSOM(ReadInput(fcs,compensate = F,transform = F), 
                  xdim=input$Xgrid, ydim=input$Ygrid, colsToUse=c(1:dim(fcs)[2]))
        mst <- BuildMST(fs)
        if(input$methodsMeta != "None"){
          cluster<-MetaClustering(data=mst$map$codes,method=input$methodsMeta, nClus=input$maxCluster)[mst$map$mapping[,1]]
        } else {
          cluster <- mst$map$mapping[,1]
        }
         clustering <- c(clustering, list(cluster))
      }
      if(input$clustering_meth == "flowPeaks"){
        fp <- flowPeaks(data, tol=input$tol)
        cluster <- fp$peaks.cluster
        clustering <- c(clustering, list(cluster))
      }
    }
    values$clustering <- clustering
    # values$centres <- centres
    # if(length(centres)==length(flow.frames))names(values$centres) <- names(flow.frames)
    names(values$clustering) <- names(flow.frames)
    progress$close()
  })
  
  output$fileClustered <- reactive({
    return(length(values$clustering)>0)
  })
  
  outputOptions(output, "fileClustered", suspendWhenHidden = FALSE)
  
  output$selectXCluster <- renderUI({
    selectInput("XClusterPreviw","Select X", choices=values$params,selected=values$params[1])
  })
  output$selectYCluster <- renderUI({
    selectInput("YClusterPreviw","Select Y", choices=values$params,selected=values$params[2])
  })
  
  output$IDClusteredPreviewxOutput <- renderUI({
    selectInput("idpreviewClustered","Select Files 1",choice=names(values$flow.frames))
  })
  
  output$IDClusteredPreviewxOutput2 <- renderUI({
    selectInput("idpreviewClustered2","Select Files 2",choice=names(values$flow.frames))
  })
  
  output$plotPreviewCluster <- renderPlot({
    if(length(values$clustering)<1){return(NULL)}
    if(is.null(input$XClusterPreviw) || is.null(input$YClusterPreviw)) return(NULL)
    col <- rainbow(length(unique(values$clustering[[input$idpreviewClustered]])))[values$clustering[[input$idpreviewClustered]]]
    graphics::plot(values$flow.frames[[input$idpreviewClustered]]@exprs[,c(input$XClusterPreviw,input$YClusterPreviw)],main=input$idpreviewClustered, pch=".",col=col)
  },width=600, height=600)
  
  output$plotPreviewCluster2 <- renderPlot({
    if(length(values$clustering)<1){return(NULL)}
    if(is.null(input$XClusterPreviw) || is.null(input$YClusterPreviw)) return(NULL)
    graphics::plot(values$flow.frames[[input$idpreviewClustered2]]@exprs[,c(input$XClusterPreviw,input$YClusterPreviw)],main=input$idpreviewClustered2, pch=".",col=colors()[values$clustering[[input$idpreviewClustered2]]])
  },width=600, height=600)
  
  observeEvent(input[["enrich_cluster"]],{
    c <- TRUE
    for(i in input$clustering_meth){
      nw <- input[[paste0("dim_",i)]]
      if(length(grep(nw,colnames(values$flow.frames[[1]])))>0){
        c <- FALSE
      }
    }
    if(c){
    flow.frames <- c()
    for(i in names(values$flow.frames)){
      clusters <- values$clustering[[i]]
      names(clusters) <- input[[paste0("dim_",input$clustering_meth)]]
      fcs <- enrich.FCS.CIPHE(values$flow.frames[[i]], clusters, nw.names=input[[paste0("dim_",input$clustering_meth)]])
      flow.frames <- c(flow.frames, fcs)
    }
    names(flow.frames) <- names(values$flow.frames)
    values$flow.frames <- flow.frames
    values$log <- rbind(values$log, paste0("Clustering :",input$clustering_meth))
    showNotification("Add new dimensions to FCS",type = "message")
    removeModal()
    getParams()
    } else {
      showNotification("One Dimension Already Exists",type = "error")
    }
  })
  
  observeEvent(input$addClusterID,{
    if(is.null(values$clustering))return(NULL)
    showModal(checkDimName(input$clustering_meth, values$flow.frames[[1]], button_id = "enrich_cluster"))
  })
 
  output$clusterCenters <- downloadHandler(
    filename = function() {
      paste("data_", Sys.Date(), "_centers.zip", sep="")
    },
    content = function(file) {
      progress <- Progress$new()
      fcs <- values$centres
      names.fcs <- names(fcs)
      tmpdir <- tempdir()
      setwd(tempdir())
      csv <- c()
      progress$set(message="Zip files", value=0.5)
      for(i in c(1:length(fcs))){
        cs <- fcs[[i]]
        path <- paste0(names.fcs[[i]],".csv")
        csv <- c(csv, path)
        write.csv(cs,path)
      }
      zip(zipfile=file, files=csv)
      setwd(getwd())
      progress$close()
    }
  )
  
##### CLUSTERING ANALYSIS ########################
  
  output$uiSelectClusterId <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    selectInput("selectClusterID","Cluster Params",choices=values$params, multiple=FALSE)
  })
  
  output$uiClAnalXparam <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    selectInput("clAnalXParam","X",choices=values$params, multiple=FALSE,selected=1)
  })
  
  output$uiClAnalYparam <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    selectInput("clAnalYParam","Y",choices=values$params, multiple=FALSE,selected=2)
  })
  
  output$useViewPlot <- renderUI({
    outputPlotObject <- lapply(values$flow.frames, function(i){
      plotObject <- renderPlot({
        par(mar=c(5,5,2,2))
        if(is.null(input$clAnalXParam) || is.null(input$clAnalYParam)) return(NULL)
        if(!is.null(input$selectClusterID) && checkAnnot(i,input$selectClusterID,50)==FALSE){
          dat <- as.vector(unlist(i@exprs[,input$selectClusterID]))
          color <- rainbow(length(unique(dat)))
          plot(i@exprs[,c(input$clAnalXParam,input$clAnalYParam)],pch=".",col=color[dat])
        } else {
          plotDens(i,c(input$clAnalXParam,input$clAnalYParam))
        }
      },width = 300,height = 300)
      return(plotObject)
    })
    do.call(tagList, outputPlotObject)
  })
  
  # output$uiMultiPlotXParam <- renderUI({
  #   if(is.null(values$flow.frames))return(NULL)
  #   selectInput("multiPlotXParam","X",choices=values$params, multiple=FALSE,selected=1)
  # })
  
  # output$uiMultiPlotYParam <- renderUI({
  #   if(is.null(values$flow.frames))return(NULL)
  #   selectInput("multiPlotYParam","Y",choices=values$params, multiple=FALSE,selected=2)
  # })
  
  output$uifileViewMDP <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    selectInput("fileViewMDP","Files",choices=names(values$flow.frames),multiple=FALSE)
  })
  
  output$markerDensityPlot <- renderUI({
    if(is.null(input$clAnalXParam) || is.null(input$clAnalYParam))return(NULL)
    outputMultiPlot <- lapply(input$selectParamsHeatmap,function(i){
      data <- values$flow.frames[[input$fileViewMDP]]@exprs[,c(input$clAnalXParam,input$clAnalYParam)]
      palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
      palette <- palette(20)[as.numeric(cut(values$flow.frames[[input$fileViewMDP]]@exprs[,i],breaks=20))]
      return(renderPlot({
        par(mar=c(5,5,2,2))
        plot(data, pch=".",col=palette,main=i,xlab=input$clAnalXParam,ylab=input$clAnalYParam)
      },width=250,height = 250))
    })
    outputMultiPlot <- do.call(tagList,outputMultiPlot)
    return(outputMultiPlot)
  })
  
  output$uiSelectClusterIdUnique <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    if(is.null(input$selectClusterID)) return(NULL)
    if(checkAnnot(values$flow.frames[[1]],params = input$selectClusterID,th = 226)==FALSE){
      selectInput("selectClusterIdUnique","Select cluster",choices=unique(values$flow.frames[[1]]@exprs[,input$selectClusterID]),multiple=TRUE,selected =unique(values$flow.frames[[1]]@exprs[,input$selectClusterID]) )
    }
  })
  
  output$uiSelectSplitParams <- renderUI({
    if(is.null(values$flow.frames))return(NULL)
    if(input$checkSplitClAnal){
      selectInput("selectSplitParams","Split by",choices=values$params, multiple=F)
    } else {
      return(NULL)
    }
  })
  
  output$uiSelectParamsHeatmap <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("selectParamsHeatmap","Select Params",choices=values$params, multiple=TRUE)
  })
  
  observeEvent(input$allParamsClustAnal,{
    updateSelectInput(session,"selectParamsHeatmap",choices = values$params, selected=values$params)
  })
  
  observeEvent(input$clearParamsClustAnal,{
    updateSelectInput(session,"selectParamsHeatmap",choices = values$params, selected=NULL)
  })
  
  observeEvent(input$addGrepParamsClustAnal,{
    if(is.null(input$grepParamsClustAnal))return(NULL)
    id <- grep(input$grepParamsClustAnal,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session,"selectParamsHeatmap",choices = values$params, selected=values$params[id])
  })
  
  observeEvent(input$runHcl,{
    if(is.null(input$selectClusterIdUnique)) return(NULL)
    progress <- Progress$new()
    progress$set(message="Compute Heatmap...",value=0.9)
    
    output$heatmapCluster <- renderUI({
      
      outputPlotObject <-lapply(values$flow.frames,function(i){
        
        plotObject <- renderPlotly({
          mat <- lapply(input$selectClusterIdUnique, function(j){
            dat <-i@exprs[,input$selectParamsHeatmap]
            colnames(dat) <- names(values$params)[which(values$params%in%input$selectParamsHeatmap)]
            res <- apply(dat[which(i@exprs[,input$selectClusterID]==j),],2,mean)
            return(res)
          })
          mat <- do.call(rbind,mat)
          row.names(mat) <- input$selectClusterIdUnique
          heatmaply(as.matrix(mat),
                    colors = plasma(256, alpha = 1, begin = 0, end = 1, direction = 1),
                    scale = "none", Colv = "row", dendrogram = "row",
                    height=1000,width = 1000)
        })
        return(plotObject)
      })
      outputPlotObject <- do.call(tagList, outputPlotObject)
      return(outputPlotObject)
    })
    progress$close()
  })
  
##### FLOWAI #####################################
  
  output$selectFlowAIStep <- renderUI({
    return(
      checkboxGroupInput("FlowAIStep", label = h3("FlowAI Steps"),
        choices = list("Flow Rate" = "FR",
                        "Signal Acquisition" ="FS",
                        "Dynamic Range" = "FM"),
        selected = c("FR","FM","FS")
      )
    )
  })
  
  observeEvent(input$FlowAIStep,{
    output$selectChFM <- renderUI({
      if(length(grep("FS",input$FlowAIStep))>0){
        return(selectInput("ChFS","Remove Channels S.A.", choices = values$params,multiple=TRUE))
      } else {return(NULL)}
    })
    output$selectChFS <- renderUI({
      if(length(grep("FM",input$FlowAIStep))>0){
        return(selectInput("ChFM","Remove Channels D.R.", choices = values$params,multiple=TRUE))
      } else {return(NULL)}
    })
  })
  
  observeEvent(input$clean,{
    progress <- Progress$new()
    progress$set(message="Cleanning in progress...",value=0.5)
    if(!is.null(input$FlowAIStep)){
      remove_from <- input$FlowAIStep
      if(length(remove_from) == 3){
        remove_from <- "all"
      } else {
        remove_from <- paste0(remove_from,collapse = "_")
      }
      ChRemoveFS <- input$ChFS
      ChFM <- input$ChFM
      
      flow.frames <- values$flow.frames
      flowAI <- lapply(flow.frames, function(x){
        res <- flow_auto_qc.CIPHE(x,output = 2, folder_results = FALSE,
                            html_report = FALSE, mini_report = FALSE, fcs_QC = FALSE,
                            remove_from = remove_from,ChExcludeFS = ChRemoveFS,
                            ChExcludeFM = ChFM
        )
        return(res)
      })
      progress$set(message = "Cleanning in progress...",value = 1)
      progress$close()
      values$flowAI <- lapply(flowAI, function(i){return(i$flow.frame)})
      values$data_flowAI <- flowAI
      names(values$flowAI) <- names(values$flow.frames)
    }
  })
  
  output$selectFlowRate <- renderUI({
    if(is.null(values$data_flowAI)){return(NULL)}
    selectInput("idFlowRare","Select File",choice=names(values$data_flowAI))
  })
  
  output$selectSignal <- renderUI({
    if(is.null(values$data_flowAI)){return(NULL)}
    selectInput("idSignal","Select File",choice=names(values$data_flowAI))
  })
  
  output$selectDynamic <- renderUI({
    if(is.null(values$data_flowAI)){return(NULL)}
    selectInput("idDynamic","Select File",choice=names(values$data_flowAI))
  })
  
  output$plotFlowRate <- renderPlot({
    if(is.null(input$idFlowRare)){return(NULL)}
    return(flow_rate_plot(values$data_flowAI[[input$idFlowRare]]$data[[1]]))
  })
  
  output$plotSignal <- renderPlot({
    if(is.null(input$idSignal)){return(NULL)}
    return(flow_signal_plot(values$data_flowAI[[input$idSignal]]$data[[2]]))
  })
  
  output$plotDynamic <- renderPlot({
    if(is.null(input$idDynamic)){return(NULL)}
    return(flow_margin_plot(values$data_flowAI[[input$idDynamic]]$data[[3]]))
  })
  
  output$fileCleanAI <- reactive({
    return(length(values$flowAI)>0)
  })
  
  outputOptions(output, "fileCleanAI", suspendWhenHidden = FALSE)
  
  observe({
    if(input$tabs != "flowai"){return(NULL)}
    if(is.null(values$flowAI)){return(NULL)}
    progress <- Progress$new()
    progress$set(message = "compute quality", value = 1)
    output$tableQC_ui <- renderTable({
      flow.frames.q <- values$flowAI
      res <- lapply(flow.frames.q, function(i){
        lo <- length(which(i@exprs[,dim(i)[2]]>10000))
        hi <- length(which(i@exprs[,dim(i)[2]]<10000))
        tot <- dim(i)[1]
        freq.lo <- round(100*(lo/tot),2)
        freq.hi <- round(100*(hi/tot),2)
        ret <- c(1,tot,lo,hi,freq.lo,freq.hi)
        return(ret)
      })
      
      mat <- do.call(rbind, res)
      # mat <- iris
      mat[,1] <- names(values$flow.frames)
      colnames(mat) <- c("filesnames","total #","lowQC #","highQC #","lowQC %","highQC %")
      return(data.frame(mat, check.names = FALSE))
    },  rownames = FALSE, colnames = TRUE)
    progress$close()
  })
  
  output$selectXFlowAI <- renderUI({
    selectInput("XFlowAIPreviw","Select X", choices=values$params,selected=values$params[1])
  })
  output$selectYFlowAI <- renderUI({
    selectInput("YFlowAIPreviw","Select Y", choices=values$params,selected=values$params[2])
  })
  
  output$IDFlowAIPreviewxOutput <- renderUI({
    selectInput("idpreviewFlowAI","Select Files 1",choice=names(values$flowAI))
  })
  
  output$plotPreviewFlowAI <- renderPlot({
    if(length(values$flowAI)<1){return(NULL)}
    color <- rep(3, dim(values$flowAI[[input$idpreviewFlowAI]])[1])
    print(length(color))
    color[which(values$flowAI[[input$idpreviewFlowAI]]@exprs[,"QCvector"]>10000)] <- 2
    plot(values$flow.frames[[input$idpreviewFlowAI]]@exprs[,c(input$XFlowAIPreviw,input$YFlowAIPreviw)],
         main=input$idpreviewFlowAI, pch=".", col=color)
  },width=500, height=500)
  
  observeEvent(input$addFlowAI,{
    flow.frames <- c()
    for(i in names(values$flow.frames)){
      clusters <- values$flowAI[[i]]@exprs[,"QCvector"]
      names(clusters) <- input$flagQC
      fcs <- enrich.FCS.CIPHE(values$flow.frames[[i]], clusters, nw.names=input$flagQC)
      flow.frames <- c(flow.frames, fcs)
    }
    names(flow.frames) <- names(values$flow.frames)
    values$flow.frames <- flow.frames
    values$log <- rbind(values$log,paste0("FlowAI Enrich"))
    getParams()
  })

##### FLOWCLEAN ##################################
  
  observe({
    if(input$tabs != "flowclean") return(NULL)
    if(is.null(values$flow.frames)) return(NULL)
    updateSelectInput(session, "vectMarkers", choices=values$params)
  })
  
  observeEvent(input$allVectMarkers,{
    updateSelectInput(session, "vectMarkers", choices=values$params,selected=values$params)
  })
  
  observeEvent(input$clearVectMarkers,{
    updateSelectInput(session, "vectMarkers", choices=values$params,selected=NULL)
  })
  
  observeEvent(input$addGrepVectMarkers,{
    if(is.null(input$grepVectMarkers))return(NULL)
    id <- grep(input$grepVectMarkers,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session, "vectMarkers", choices=values$params,selected=values$params[id])
  })
  
  observeEvent(input$runFlowClean,{
    if(is.null(input$vectMarkers))return(NULL)
    progress <- Progress$new()
    v <- 0
    flow.frames <- values$flow.frames
    res <- lapply(flow.frames, function(i){
      v <<- v + 1
      progress$set(message="flowClean progress...",value=v/length(values$flow.frames))
      id <- (colnames(i)%in%input$vectMarkers)
      return(clean(i, vectMarkers=id, diagnostic=FALSE, 
                   binSize=as.numeric(input$binSize),nCellCutoff=as.numeric(input$nCellCutoff),
              filePrefixWithDir=NULL, ext=NULL)@exprs[,"GoodVsBad"])
    })
    names(res) <- values$names.files
    values$flowClean <- res
    progress$close()
  })

  observe({
    if(input$tabs != "flowclean") return(NULL)
    if(is.null(values$flowClean)) return(NULL)
    output$uiSelectXflowClean <- renderUI({
      selectInput("XflowClean","X Param",choices=values$params, multiple=FALSE,selected=values$params[[1]])
    })
    output$uiSelectYflowClean <- renderUI({
      selectInput("YflowClean","Y Param",choices=values$params, multiple=FALSE,selected=values$params[[2]])
    })
    output$uiSelectIDflowClean <- renderUI({
      selectInput("IDflowClean","Files",choices=values$names.files, multiple=FALSE)
    })
    output$tablResultFlowCLEAN <- renderTable({
      row <- values$names.files
      res <- lapply(values$flowClean,function(j){
        g <- length(which(j<=10000))
        b <- length(which(j>10000))
        return(c(g,b,g/length(j)*100,b/length(j)*100))
      })
      mat <- do.call(rbind,res)
      colnames(mat) <- c("Good#events","Bad#events","Good%Percents","Bad%Percents")
      row.names(mat) <- row
      return(mat)
    },rownames = TRUE,colnames = TRUE)
  })
  
  observe({
    if(input$tabs != "flowclean") return(NULL)
    if(is.null(input$XflowClean) || is.null(input$YflowClean)) return(NULL)
    fcs <- values$flow.frames[[input$IDflowClean]]
    res <- values$flowClean[[input$IDflowClean]]
    output$plotFlowClean <- renderPlot({
      col <- rep(3,length(res))
      col[which(res>10000)] <- 2
      plot(exprs(fcs)[,c(input$XflowClean,input$YflowClean)],pch=".",col=col,xlab=input$XflowClean,ylab=input$YflowClean)
    })
  })
  
  observeEvent(input$enrichFlowClean,{
    if(is.null(values$flowClean)) return(NULL)
    flow.frames <- lapply(c(1:length(values$flow.frames)),function(i){
      nw <- values$flowClean[[i]]
      fcs<- values$flow.frames[[i]]
      return(enrich.FCS.CIPHE(fcs,nw,"flowClean"))
    })
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    showNotification(ui="Enrich FlowClean Done !",type = "message")
    values$log <- rbind(values$log,paste0("FlowClean Enrich"))
    getParams()
  })
    
##### DIMENSION REDUCTION PCA ####################
  
  output$selectPCAParams <- renderUI({
    selectInput("PCAParams","Select Params",choices=values$params,multiple=TRUE)
  })
  
  observeEvent(input$addAllPCADim,{
    updateSelectInput(session,"PCAParams",choices=values$params,selected=values$params)
  })
  
  observeEvent(input$addGrepPCA,{
    if(is.null(input$grepPCA)) return(NULL)
    id <- grep(input$grepPCA,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session,"PCAParams",choices=values$params, selected=values$params[id])
  })
  
  observeEvent(input$clearAllPCADim,{
    updateSelectInput(session,"PCAParams",choices=values$params,selected=NULL)
  })
  
  observeEvent(input$runPCA,{
    if(is.null(input$PCAParams)){return(NULL)}
    progress <- Progress$new()
    progress$set(message="Dim reduction", value=0.5)
    flow.frames <- values$flow.frames
    dim.reduc <- list()
    for(i in names(flow.frames)){
      data <- PCA(flow.frames[[i]]@exprs[,input$PCAParams],graph=F,ncp=length(input$PCAParams))
      data$svd$U <- data.frame(data$svd$U)
      colnames(data$svd$U) <- sprintf("PCA%d", c(1:length(input$PCAParams)))
      dim.reduc <- c(dim.reduc,list(data))
    }
    names(dim.reduc) <- names(flow.frames)
    values$pca.reduc <- dim.reduc
    progress$close()
  })
  
  output$fileReducPCA <- reactive({
    return(length(values$pca.reduc)>0)
  })
  
  outputOptions(output, "fileReducPCA", suspendWhenHidden = FALSE)
  
  output$ZPCAPlot <- renderUI({
    if(is.null(values$params)){return(NULL)}
    selectInput("ZPCA","Select Color",choices=values$params)
  })
  
  output$XPCAPlot <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("XPCA","X Params",choices=colnames(values$pca.reduc[[1]]$svd$U),selected=colnames(values$pca.reduc[[1]]$svd$U)[[1]])
  })
  
  output$YPCAPlot <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("YPCA","Y Params",choices=colnames(values$pca.reduc[[1]]$svd$U),selected=colnames(values$pca.reduc[[1]]$svd$U)[[2]])
  })
  
  output$selectPCA <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("PCAID","Select File",choices=names(values$flow.frames))
  })
  
  output$selectPCA2 <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("PCAID2","Select File",choices=names(values$flow.frames))
  })
  
  output$selectPCAToAdd <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("PCASelect","Add Params",choices=colnames(values$pca.reduc[[1]]$svd$U),multiple = TRUE)
  })
  
  output$PCAPlot <- renderPlot({
    if(is.null(values$pca.reduc)){return(NULL)}
    if(is.null(input$PCAID)){return(NULL)}
    z <- values$flow.frames[[input$PCAID]]@exprs[,input$ZPCA]
    dat <- data.frame(values$pca.reduc[[input$PCAID]]$svd$U[,c(input$XPCA,input$YPCA)])
    colnames(dat) <- c("x","y")
    palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
    col <- palette(20)[as.numeric(cut(z,breaks=20))]
    plot(dat[,c("x","y")],col=col,pch=20)
  },width = 600, height=600)
  
  # output$PCAPlot2 <- renderPlot({
  #   if(is.null(values$pca.reduc)){return(NULL)}
  #   if(is.null(input$PCAID2)){return(NULL)}
  #   col <- sort(values$flow.frames[[input$PCAID2]]@exprs[,input$ZPCA])
  #   dat <- data.frame(values$pca.reduc[[input$PCAID2]]$svd$U[,c(input$XPCA,input$YPCA)])
  #   colnames(dat) <- c("x","y")
  #   qplot(x, y, data=dat, colour=col) + scale_colour_gradient(low="red", high="blue")
  # })
  
  output$uiInfoPCA <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    selectInput("infoPCAID","Select Info",choices = c("Variance","Contribution","Qualit??","Contribution"),multiple=FALSE)
  })
  
  output$uiContribAxes <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    if(input$infoPCAID != "Contribution"){return(NULL)}
    selectInput("contribAxes","Axes",choices=c(1:dim(values$flow.frames[[1]])[2]),selected=1)
  })
  
  output$uiContribTop <- renderUI({
    if(is.null(values$pca.reduc)){return(NULL)}
    if(input$infoPCAID != "Contribution"){return(NULL)}
    selectInput("contribTop","Top",choices=c(1:dim(values$flow.frames[[1]])[2]),selected=1)
  })
  
  output$plotInfo1 <- renderPlot({
    if(is.null(values$pca.reduc)){return(NULL)}
    if(is.null(input$PCAID)){return(NULL)}
    if(input$infoPCAID == "Variance"){
      return(fviz_eig(values$pca.reduc[[input$PCAID]], addlabels = TRUE))
    }
    if(input$infoPCAID == "Contribution"){
      return(fviz_contrib(values$pca.reduc[[input$PCAID]], choice = "var",
                          axes = as.numeric(input$contribAxes),top=as.numeric(input$contribTop)))
    }
  })
  
  output$plotInfo2 <- renderPlot({
    if(is.null(values$pca.reduc)){return(NULL)}
    if(is.null(input$PCAID2)){return(NULL)}
    return(fviz_eig(values$pca.reduc[[input$PCAID2]], addlabels = TRUE))
  })
  
  observeEvent(input[["enrich_pca"]],{
    if(is.null(input$PCASelect))return(NULL)
    c <- TRUE
    for(i in input$PCASelect){
      nw <- input[[paste0("dim_",i)]]
      if(length(grep(nw,colnames(values$flow.frames[[1]])))>0){
        c <- FALSE
      }
    }
    if(c){
      progress <- Progress$new()
      progress$set(message="Add new dimension...",value=0.9)
      flow.frames <- lapply(c(1:length(values$flow.frames)),function(j){
        fcs <- values$flow.frames[[j]]
        pca <- values$pca.reduc[[j]]
        for(i in input$PCASelect){
          n <- input[[paste0("dim_",i)]]
          fcs <- enrich.FCS.CIPHE(fcs, pca$svd$U[,i],n)
        }
        return(fcs)
      })
      values$flow.frames <- flow.frames
      names(values$flow.frames) <- values$names.files
      getParams()
      progress$close()
      values$log <- rbind(values$log,paste0("PCA Enrichment"))
      showNotification("Add new dimensions to FCS",type = "message")
      removeModal()
    } else {
      showNotification("One Dimension Already Exists",type = "error")
    }
  })
  
  observeEvent(input$addPCA,{
    if(is.null(values$pca.reduc))return(NULL)
    if(is.null(input$PCASelect))return(NULL)
    showModal(checkDimName(input$PCASelect, values$flow.frames[[1]],button_id="enrich_pca"))
  })
  
##### DIMENSION REDUCTION UNLINEAR ###############
  
  observeEvent(input$addAllUMAP,{
    updateSelectInput(session,"UMAPParams",choices=values$params, selected=values$params)
  })
  
  observeEvent(input$clearAllUMAP,{
    updateSelectInput(session,"UMAPParams",choices=values$params, selected=NULL)
  })
  
  observeEvent(input$addGrepUMAP,{
    if(is.null(input$grepUMAP)) return(NULL)
    id <- grep(input$grepUMAP,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session,"UMAPParams",choices=values$params, selected=values$params[id])
  })
  
  output$selectUMAPParams <- renderUI({
    selectInput("UMAPParams","Select Params",choices=values$params,multiple=TRUE)
  })
  
  observeEvent(input$runUMAP,{
    if(is.null(input$UMAPParams)){return(NULL)}
    progress <- Progress$new()
    progress$set(message="Dim reduction", value=0)
    flow.frames <- values$flow.frames
    dim.reduc <- list()
    j <- 0
    m <- length(names(flow.frames))
    for(i in names(flow.frames)){
      j <- j+1
      progress$set(message=paste0("Dim reduction ",j,'/',m), value=j/m)
      if(input$unlineartReduc == "UMAP"){
        if(input$advOptUmap == FALSE){
          print("umap")
          data <- data.frame(umap::umap(flow.frames[[i]]@exprs[,input$UMAPParams])$layout)
          colnames(data) <- c("UMAP1","UMAP2")
          dim.reduc <- c(dim.reduc,list(data))
        } else {
          print("uwot")
          data <- data.frame(uwot::umap(flow.frames[[i]]@exprs[,input$UMAPParams],n_neighbors = input$n_neighbors,metric = input$metric))
          colnames(data) <- c("UMAP1","UMAP2")
          dim.reduc <- c(dim.reduc,list(data))
        }
      }
      if(input$unlineartReduc == "R-TSNE"){
        tsne <- Rtsne(unique(flow.frames[[i]]@exprs[,input$UMAPParams]),
          perplexity=input$knn, max_iter=input$iter, theta=input$theta, eta=input$eta)

        unique.id <- which(duplicated(flow.frames[[i]]@exprs[,input$UMAPParams]))

        data <- tsne$Y
        for(i in unique.id){
          data <- rbind(data[1:(i-1),],c(0,0),data[i:dim(data)[1],])
        }


        # data <- data.frame(rbind(tsne$Y,as.matrix(mat.add.unique)))
        colnames(data) <- c("tSNE1","tSNE2")
      }
      if(input$unlineartReduc == "IsoMap"){
        iso <- isomap(flow.frames[[i]]@exprs[,input$UMAPParams])
        data <- data.frame(iso)
        colnames(data) <- c("IsoMap1","IsoMap2")
      }
      if(input$unlineartReduc == "DiffusionMap"){
        dm <- DiffusionMap(unique(flow.frames[[i]]@exprs[,input$UMAPParams]))
        add.unique <- rep(0,length(flow.frames[[i]]@exprs[,input$UMAPParams])-length(unique(flow.frames[[i]]@exprs[,input$UMAPParams])))
        mat.add.unique <- data.frame("1"=add.unique,"2"=add.unique)
        data  <- data.frame(rbind(cbind(dm$DC1,dm$DC2),mat.add.unique))
        colnames(data) <- c("DiffusionMap1","DiffusionMap2")
      }
      if(input$unlineartReduc == "EmbedSOM"){
        fcs <- flowFrame(flow.frames[[i]]@exprs[,input$UMAPParams])
        fs <- FlowSOM::BuildSOM(ReadInput(fcs,compensate = F,transform = F), xdim=20, ydim=20, colsToUse=c(1:dim(fcs)[2]))
        data <- data.frame(EmbedSOM::EmbedSOM(fs))
        colnames(data) <- c("EmbedSOM1","EmbedSOM2")
      }
      dim.reduc <- c(dim.reduc,list(data))
    }

    names(dim.reduc) <- names(flow.frames)
    values$umap.reduc <- dim.reduc
    progress$close()
  })
  
  output$fileReducUMAP <- reactive({
    return(length(values$umap.reduc)>0)
  })
  
  outputOptions(output, "fileReducUMAP", suspendWhenHidden = FALSE)
  
  output$ZumapPlot <- renderUI({
    selectInput("ZUMAP","Select Color",choices=values$params)
  })
  
  output$XumapPlot <- renderUI({
    if(is.null(values$umap.reduc)){return(NULL)}
    selectInput("XUMAP","X Params",choices=colnames(values$umap.reduc[[1]]),selected=colnames(values$umap.reduc[[1]])[1])
  })
  
  output$YumapPlot <- renderUI({
    if(is.null(values$umap.reduc)){return(NULL)}
    selectInput("YUMAP","Y Params",choices=colnames(values$umap.reduc[[1]]),selected=colnames(values$umap.reduc[[1]])[2])
  })
  
  output$selectUMAP <- renderUI({
    if(is.null(values$umap.reduc)){return(NULL)}
    selectInput("UMAPID","Select File",choices=names(values$umap.reduc))
  })
  
  output$selectUMAP2 <- renderUI({
    if(is.null(values$umap.reduc)){return(NULL)}
    selectInput("UMAPID2","Select File",choices=names(values$umap.reduc))
  })
  
  output$selectUMAPToAdd <- renderUI({
    if(is.null(values$umap.reduc)){return(NULL)}
    selectInput("UMAPSelect","Add Params",choices=colnames(values$umap.reduc[[1]]),multiple = TRUE,selected=colnames(values$umap.reduc[[1]]))
  })
  
  output$UMAPPlot <- renderPlot({
    if(is.null(values$umap.reduc)){return(NULL)}
    if(is.null(input$UMAPID)){return(NULL)}
    z <- values$flow.frames[[input$UMAPID]]@exprs[,input$ZUMAP]
    dat <- data.frame(values$umap.reduc[[input$UMAPID]][,c(input$XUMAP,input$YUMAP)])
    colnames(dat) <- c("x","y")
    palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
    col <- palette(20)[as.numeric(cut(z,breaks=20))]
    plot(dat[,c("x","y")],col=col,pch=20)
  },width = 500, height =500)

  output$UMAPPlot2 <- renderPlot({
    if(is.null(values$umap.reduc)){return(NULL)}
    if(is.null(input$UMAPID2)){return(NULL)}
    z <- values$flow.frames[[input$UMAPID2]]@exprs[,input$ZUMAP]
    dat <- data.frame(values$umap.reduc[[input$UMAPID2]][,c(input$XUMAP,input$YUMAP)])
    colnames(dat) <- c("x","y")
    palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
    col <- palette(20)[as.numeric(cut(z,breaks=20))]
    plot(dat[,c("x","y")],col=col,pch=20)
  },width = 500, height=500)
  
  observeEvent(input[["enrich_umap"]],{
    if(is.null(input$UMAPSelect))return(NULL)
    c <- TRUE
    for(i in input$UMAPSelect){
      nw <- input[[paste0("dim_",i)]]
      if(length(grep(nw,colnames(values$flow.frames[[1]])))>0){
        c <- FALSE
      }
    }
    if(c){
      progress <- Progress$new()
      progress$set(message="Add new dimension...",value=0.9)
      flow.frames <- lapply(c(1:length(values$flow.frames)),function(j){
        fcs <- values$flow.frames[[j]]
        pca <- values$umap.reduc[[j]]
        for(i in input$UMAPSelect){
          n <- input[[paste0("dim_",i)]]
          fcs <- enrich.FCS.CIPHE(fcs, pca[,i], n)
        }
        return(fcs)
      })
      values$flow.frames <- flow.frames
      names(values$flow.frames) <- values$names.files
      getParams()
      progress$close()
      showNotification("Add new dimensions to FCS",type = "message")
      values$log <- rbind(values$log,paste0("Unlinear :",input$unlineartReduc))
      removeModal()
    } else {
      showNotification("One Dimension Already Exists",type = "error")
    }
  })
  
  observeEvent(input$addUMAP,{
    if(is.null(values$umap.reduc))return(NULL)
    if(is.null(input$UMAPSelect))return(NULL)
    showModal(checkDimName(input$UMAPSelect, values$flow.frames[[1]],button_id="enrich_umap"))
  })
  
##### BACKGATING #################################
  
  observe({
    if(input$tabs != "backgating"){return(NULL)}
    if(is.null(values$flow.frames) || is.null(input$markerX) || is.null(input$markerY) || is.null(input$annotation)) return(NULL)
    col <- rep(2, dim(values$flow.frames[[input$filesBackGating]])[1])
    if(input$annotation!= "" || input$annotation_id != "" || is.null(input$annotation) || is.null(input$annotation_id)){
      col[which(values$flow.frames[[input$filesBackGating]]@exprs[,input$annotation]==as.numeric(input$annotation_id))] <- 3
    }
    output$plotPreviewBG <- renderPlot({
      plot(values$flow.frames[[1]]@exprs[,c(input$markerX,input$markerY)],col=col, pch=".")
    })
  })
  
  observeEvent(input$annotation,{
    if(is.null(input$annotation) || is.null(values$flow.frames) || input$annotation=="")return(NULL)
    if(checkAnnot(values$flow.frames[[1]],input$annotation,th=100)){return(NULL)}
    updateSelectInput(session,"annotation_id",choices=unique(values$flow.frames[[1]]@exprs[,input$annotation]))
  })
  
  observeEvent(input$addAllMarkersBG,{
    if(is.null(values$params)) return(NULL)
    updateSelectInput(session, "markers","Select Markers",choices=values$params,selected=values$params
    )
  })
  
  output$selectMarkersBG <- renderUI({
    if(is.null(values$params)) return(NULL)
    selectInput("markers","Select Markers",choices=values$params,multiple=TRUE)
  })
  
  observeEvent(input$addAllMarkersBG,{
    if(is.null(values$params)) return(NULL)
    updateSelectInput(session, "markers","Select Markers",
                      choices=values$params,
                      selected=values$params
    )
  })
  
  output$selectFilesBackGating <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("filesBackGating","Select Files",choices=names(values$flow.frames),selected=names(values$flow.frmaes)[[1]])
  })
  
  output$UIannotation <- renderUI({
    if(is.null(values$params)) return(NULL)
    selectInput("annotation","Annotation Parameters", choices = c("",values$params))
  })
  output$UImarkerX <- renderUI({
    if(is.null(values$params)) return(NULL)
    selectInput("markerX","X Parameters", choices=values$params, selected=1)
  })
  output$UImarkerY <- renderUI({
    if(is.null(values$params)) return(NULL)
    selectInput("markerY","X Parameters", choices=values$params, selected=1)
  })
  

  observeEvent(input$runHypergate,{
    if(is.null(input$annotation_id) || is.null(input$markers)) return(NULL)
    progress <- Progress$new()
    progress$set(message="Hypergate in progress...", value=1)
    
    gate_vector <- rep(0, dim(values$flow.frames[[input$filesBackGating]])[1])
    gate_vector[which(values$flow.frames[[input$filesBackGating]]@exprs[,input$annotation]==input$annotation_id)] <- 1
    
    hg_output <- hypergate(
      xp = values$flow.frames[[input$filesBackGating]]@exprs[, input$markers],
      gate_vector = gate_vector, 
      level = 1, verbose = FALSE
    )
    
    values$hg_output <- hg_output
    
    gating_predicted <- subset_matrix_hg(hg_output, values$flow.frames[[input$filesBackGating]]@exprs[, input$markers])
    
    output$tableOutput <- renderTable({
      table(ifelse(gating_predicted, "Gated-in", "Gated-out"), 
            ifelse(gate_vector == 1, "Events of interest", "Others")
      )
    },colnames=FALSE)
    
    output$gatingStrat <- renderUI({
      l <- length(hg_output$active_channels)
      progress$set(message=l)
      iter<-0
      channels<-hg_output$active_channels
      channels<-sub("_max","",channels)
      channels<-sub("_min","",channels)
      ranges.global<-apply(values$flow.frames[[input$filesBackGating]]@exprs[,channels,drop=F],2,range)
      rownames(ranges.global)=c("min","max")
      active_events<-rep(T,nrow(values$flow.frames[[input$filesBackGating]]@exprs))
      highlight <- "red"
      cols<-rep("black",nrow(values$flow.frames[[input$filesBackGating]]@exprs))
      cols[gate_vector==1]=highlight
      
      parameters<-hg_output$pars.history
      active_parameters<-hg_output$active_channels##apply(parameters,2,function(x){x[length(x)]!=x[1]})
      parameters<-parameters[,active_parameters,drop=FALSE]
      parameters_order<-apply(parameters,2,function(x)min(which(x!=x[1])))
      parameters<-parameters[,order(parameters_order,decreasing=FALSE),drop=FALSE]
      parameters<-setNames(parameters[nrow(parameters),,drop=TRUE],colnames(parameters))
      
      channels<-sub("_max","",names(parameters))
      channels<-sub("_min","",channels)
      
      direction <- rep(2,length(parameters))
      direction[grep("_max",names(parameters))]=1
      
      plotOutputObject <- list()
      
      plotOutputObject <- lapply(seq(1,length(parameters),by=2), function(i){
        # for(i in seq(1,length(parameters),by=2)){
        if((i+1)<=length(parameters)){
          chan1<-channels[i]
          chan2<-channels[i+1]
          iter<-iter+1
          tmp <- active_events
          
          if(direction[i]==2){
            test1<-values$flow.frames[[input$filesBackGating]]@exprs[,chan1]>=parameters[i] ##If _min, events above parameter are selected
          } else {
            test1<-values$flow.frames[[input$filesBackGating]]@exprs[,chan1]<=parameters[i] ##Else events above parameter below
          }
          if(direction[i+1]==2){
            test2<-values$flow.frames[[input$filesBackGating]]@exprs[,chan2]>=parameters[i+1]
          } else {
            test2<-values$flow.frames[[input$filesBackGating]]@exprs[,chan2]<=parameters[i+1]
          }
          active_events<<-active_events&test1&test2
          
          return(renderPlot({
            plot(values$flow.frames[[input$filesBackGating]]@exprs[tmp,chan1],values$flow.frames[[input$filesBackGating]]@exprs[tmp,chan2],
                 xlab=names(values$params)[grep(chan1,values$params)],
                 ylab=names(values$params)[grep(chan2,values$params)],
                 xlim=ranges.global[,chan1],
                 ylim=ranges.global[,chan2],
                 bty="l",pch=20,cex=0.1,col=cols[tmp],main=NULL
            )
            segments(
              x0=parameters[i],
              y0=parameters[i+1],
              x1=ranges.global[direction[i],chan1],
              col="red"
            )
            segments(
              x0=parameters[i],
              y0=parameters[i+1],
              y1=ranges.global[direction[i+1],chan2],
              col="red"
            )
          }, outputArgs = list(width="250px", height ="250px")))
        }
      })
      
      if(length(parameters)%%2==1){ 
        i <- length(parameters)
        chan1<-channels[i]
        #iter<-iter+1
        if(direction[i]==2){
          test1=values$flow.frames[[input$filesBackGating]]@exprs[,chan1]>=parameters[i] ##If _min, events above parameter are selected
        } else {
          test1=values$flow.frames[[input$filesBackGating]]@exprs[,chan1]<=parameters[i] ##Else events above parameter below
        }
        active_events<<-active_events&test1
        
        plotOutputObject[[length(plotOutputObject)+1]] <- renderPlot({
          plot(values$flow.frames[[input$filesBackGating]]@exprs[active_events,chan1],main="",
               ylab=channels[i],xlab="Events index",
               ylim=ranges.global[,chan1],
               bty="l",pch=16,cex=0.1,
               col=cols[active_events]
          )
        }, outputArgs = list(width="240px", height ="260px"))
      }
      
      do.call(tagList, plotOutputObject)
      return(plotOutputObject)
    })
    
    contributions <- channels_contributions(hg_output,xp=values$flow.frames[[input$filesBackGating]]@exprs[,input$markers],gate_vector=gate_vector,level=1, beta=1)
    
    values$contributions <- contributions
    
    output$barPlot <- renderPlotly({
      channels<-hg_output$active_channels
      channels<-sub("_max","",channels)
      channels<-sub("_min","",channels)
      plot_ly(
        y=contributions,
        x=names(contributions),
        type="bar"
        ,name="F1-score deterioration when the parameter is ignored",width = 800
      )
    }) 
    
    progress$close()
  })
  
  observeEvent(input$runGateFinder,{
    if(is.null(input$annotation_id) || is.null(input$markers) || is.null(input$annotation)) return(NULL)
    progress <- Progress$new()
    progress$set(message="GateFinder in progress...", value=1)
    
    mat <- values$flow.frames[[1]]@exprs[,input$markers]
    colnames(mat) <- input$markers
    targetpop <- (values$flow.frames[[1]]@exprs[,input$annotation] == input$annotation_id)
    
    results<-GateFinder(mat, targetpop)
    
    output$gatingStrat2 <-renderPlot({
      if(length(results@gates)<=3){
        row <- 1 ; col <- length(results@gates);
      } else {
        row <- length(results@gates)/3; col <- 3;
      }
      GateFinder::plot(mat, results, c(row,col), targetpop)
    })
    
    output$linePlot <- renderPlot({
      plot(results)
    })
    
    progress$close()
    
  })

##### OneSENSE ###################################
  
  output$uiSelectOnsenseFile <- renderUI({
    if(is.null(values$flow.frames)) return(NULL)
    selectInput("selectOnesenseFile","Select Files",choices=names(values$flow.frames),multiple=F)
  })
  
  observeEvent(input$firstMarkers, {
    if(is.null(values$params)){return(NULL)}
    myParam <- values$params
    updateSelectInput(session, "markers1", choices = myParam)
    updateSelectInput(session, "markers2", choices = myParam)
    # updateCheckboxGroupInput(session, "markers3", choices = myParam)
  })
  
  observeEvent(input$addAllMarkers1,{
    updateSelectInput(session,"markers1",choices=values$params,selected=values$params)
  })
  
  observeEvent(input$addAllMarkers2,{
    updateSelectInput(session,"markers2",choices=values$params, selected=values$params)
  })
  
  observeEvent(input$clearMarkers1,{
    updateSelectInput(session,"markers1",choices=values$params, selected=NULL)
  })
  
  observeEvent(input$clearMarkers2,{
    updateSelectInput(session,"markers2",choices=values$params, selected=NULL)
  })
  
  observeEvent(input$addAllGrepMarkers1,{
    if(is.null(values$flow.frames) || is.null(input$grepMarkers1)) return(NULL)
    id <- grep(input$grepMarkers1,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session,"markers1",choices=values$params, selected=values$params[id])
  })
  
  observeEvent(input$addAllGrepMarkers2,{
    if(is.null(values$flow.frames) || is.null(input$grepMarkers2)) return(NULL)
    id <- grep(input$grepMarkers2,names(values$params))
    if(length(id)==0)return(NULL)
    updateSelectInput(session,"markers2",choices=values$params, selected=values$params[id])
  })
  
  observeEvent(input$updatenamescsv, {
    myParam <- values$params
    paras1 <- c(input$markers1)
    input1 <- rep("", length(myParam))
    input1[match(paras1, myParam)] <- "Y"
    names1 <- data.frame(Marker = myParam, input = input1)
    
    paras2 <- c(input$markers2)
    input2 <- rep("", length(myParam))
    input2[match(paras2, myParam)] <- "Y"
    names2 <- data.frame(Marker = myParam, input = input2)
    
    final <- cbind(names1, input2 = names2[, 2])
    
    row.names(final) <- myParam
    values$names.csv <- final
    showNotification("Marker Updated",type="message")
    print(final)
  })
  
  observeEvent(input$submit, {
    progress <- Progress$new()
    progress$set(message="Run FCStSNE",value=0.5)
    res <- FCStSNE.CIPHE(
      LoaderPATH = values$flow.frames[[1]],
      ceil = input$num, FNnames = values$names.csv , OutputSuffix = "Out",
      DotSNE = TRUE, DoOneSENSE = TRUE, Bins = input$bins, meth=input$oneSENSEMeth
    )
    
    values$res1.onesense <- res[[1]]
    output$testOneSense1 <- renderPlot({
      res[[2]]
    })
    
    progress$set(message="Run OnesMapperFlour",value=0.9)
    
    list.plot.output <- OneSmapperFlour.CIPHE(
      LoaderPATH = values$res1.onesense, Fname = values$names.csv,
      Bins = input$bins,
      doCoords = FALSE,
      doFreq = FALSE
    )
    values$list.plot.output <- list.plot.output
    print(length(list.plot.output))
    
    output$testOneSense2 <- renderPlotly({
      list.plot.output[[3]]
    })
    
    output$testOneSense3 <- renderPlotly({
      list.plot.output[[2]]
    })
    
    output$testOneSense4 <- renderPlotly({
      list.plot.output[[1]]
    })
    
    output$dotPlot <- renderPlotly({
      list.plot.output[[3]]
    })
    
    output$heatmap1 <- renderPlotly({
      list.plot.output[[2]]
    })
  
    output$heatmap2 <- renderPlotly({
      list.plot.output[[1]]
    })
    
    progress$close()
    showNotification("One SENSE Done",type = "message")
  })
  
  observe({
    if(input$tabs != "onesens"){return(NULL)}
    if(is.null(values$list.plot.output)) return(NULL)
    output$downloadHeatMap1 <- downloadHandler(
      filename=function(){
        return("hetmap1.csv")
      },
      content = function(file){
        write.csv(values$list.plot.output[[4]],file)
      }
    )
    output$downloadHeatMap2 <- downloadHandler(
      filename=function(){
        return("hetmap2.csv")
      },
      content = function(file){
        write.csv(values$list.plot.output[[5]],file)
      }
    )
  })
  
  observeEvent(input[["enrich_onesense"]],{
    if(is.null(values$res1.onesense))return(NULL)
    c <- TRUE
    for(i in c("input1","input2")){
      nw <- input[[paste0("dim_",i)]]
      if(length(grep(nw,colnames(values$flow.frames[[1]])))>0){
        c <- FALSE
      }
    }
    if(c){
      flow.frames <- values$flow.frames
      flow.frames <- lapply(c(1:length(flow.frames)), function(i){
        new_dim1 <- values$res1.onesense[[i]]@exprs[,"input"]
        new_dim2 <- values$res1.onesense[[i]]@exprs[,"input2"]
        fcs <- enrich.FCS.CIPHE(values$flow.frames[[i]],new_dim1,input[[paste0("dim_input1")]])
        fcs <- enrich.FCS.CIPHE(fcs,new_dim2,input[[paste0("dim_input2")]])
      })
      values$flow.frames <- flow.frames
      names(values$flow.frames) <- values$names.files
      getParams()
      showNotification("Enrich OneSENSE done", type="message")
      values$log <- rbind(values$log, paste0("Reduc Dim : oneSENSE"))
      removeModal()
    } else {
      showNotification("One Dimension Already Exists",type = "error")
    }
  })
  
  observeEvent(input$enrichOneSENSE,{
    if(is.null(values$res1.onesense))return(NULL)
    showModal(checkDimName(c("input1","input2"), values$flow.frames[[1]],button_id="enrich_onesense"))
  })
  
##### MANUAL ANNOTATION ##########################
  
  observe({
    if(input$tabs != "manuannot"){return(NULL)}
    if(is.null(values$flow.frames)) return(NULL)
    if(is.null(values$params)) return(NULL)
    updateSelectInput(session,"xAnnotPlot",choices = values$params,selected=values$params[[1]])
    updateSelectInput(session,"yAnnotPlot",choices = values$params,selected=values$params[[2]])
    updateSelectInput(session,"zAnnotPlot",choices = values$params,selected=values$params[[length(values$params)]])
    updateSelectInput(session,"xAnnotPlotSub",choices = values$params,selected=values$params[[1]])
    updateSelectInput(session,"yAnnotPlotSub",choices = values$params,selected=values$params[[2]])
    updateSelectInput(session,"zAnnotPlotSub",choices = values$params,selected=values$params[[length(values$params)]])
  })
  
  observeEvent(input$zAnnotPlot,{
    if(is.null(values$flow.frames))return(NULL)
    if(checkAnnot(values$flow.frames[[1]],input$zAnnotPlot,50)==TRUE){
      showNotification("More than 50 unique value in your Z Params", type = "warning")
    }
  })
  
  observeEvent(input$zAnnotPlotSub,{
    if(is.null(values$flow.frames))return(NULL)
    if(checkAnnot(values$flow.frames[[1]],input$zAnnotPlotSub,50)==TRUE){
      showNotification("More than 50 unique value in your Z Params", type = "warning")
    }
  })
  
  output$annotPlot <- renderPlotly({
    if(is.null(input$yAnnotPlot) || is.null(input$xAnnotPlot))return(NULL)
    data <- data.frame(values$flow.frames[[1]]@exprs[,c(input$xAnnotPlot,input$yAnnotPlot)])
    if(dim(data)[1]>100000){
      showNotification("More than 100 000 events", type = "warning")
      return(NULL)
    }
    z <- as.character(rep(1,dim(data)[1]))
    if(checkAnnot(values$flow.frames[[1]],input$zAnnotPlot,50)==FALSE){
      z <- as.character(values$flow.frames[[1]]@exprs[,input$zAnnotPlot])
    } 
    data <- cbind(data,z)
    colnames(data) <- c("X","Y","Z")
    plot_ly(data, x = ~X,y=~Y,color=~Z, mode='markers',marker=list(size=5))%>%layout(dragmode="select")
  })
  
  observe({
    if(input$tabs != "manuannot"){return(NULL)}
    if(is.null(values$flow.frames)) return(NULL)
    events <- event_data("plotly_selected")
    fcs <- values$flow.frames[[1]]
    if(is.null(events) || length(events)==0)return(NULL)
    id <- unlist(lapply(1:dim(events)[1], function(x) {
      return(which(round(fcs@exprs[,input$xAnnotPlot],3)==round(events[x,"x"],3) & round(fcs@exprs[,input$yAnnotPlot],3)==round(events[x,"y"],3)))
    }))
    sub <- fcs[id,]
    output$subAnnotPlot <- renderPlotly({
      data <- data.frame(sub@exprs[,c(input$xAnnotPlotSub,input$yAnnotPlotSub)])
      z <- as.character(rep(1,dim(data)[1]))
      if(checkAnnot(fcs,input$zAnnotPlotSub,50)==FALSE){
        z <- as.character(sub@exprs[,input$zAnnotPlotSub])
      }
      data <- cbind(data,z)
      colnames(data) <- c("X","Y","Z")
      plot_ly(data, x = ~X,y=~Y,color=~Z, mode='markers',marker=list(size=5))%>%layout(dragmode="select")
    })
  })
  
  observeEvent(input$addAnnotName,{
    values$manu.annot <- rep(0,dim(values$flow.frames[[1]])[1])
    names(values$manu.annot) <- input$annotName
    # shinyjs::disable("addAnnotName")
    output$uiAddPop <- renderUI({
      actionButton("addPop","Add Pop selected")
    })
  })
  
  observe({
    if(input$tabs != "manuannot"){return(NULL)}
    output$annotTable <- renderTable({
      if(is.null(values$manu.annot)) return(NULL)
      mat <- data.frame(table(values$manu.annot))
      colnames(mat) <- c("ID","Events")
      return(mat)
    })
  })
    
  observeEvent(input$addPop,{
    print("Add")
    events <- event_data("plotly_selected")
    fcs <- values$flow.frames[[1]]
    pop <- max(values$manu.annot)+1
    id <- unlist(lapply(1:dim(events)[1], function(x) {
      return(which(round(fcs@exprs[,input$xAnnotPlot],3)==round(events[x,"x"],3) & round(fcs@exprs[,input$yAnnotPlot],3)==round(events[x,"y"],3)))
    }))
    values$manu.annot[id] <- pop
    print(table(values$manu.annot))
  })
  
  observeEvent(input$enrichAnnotName,{
    if(is.null(values$manu.annot))return(NULL)
    values$flow.frames[[1]] <- enrich.FCS.CIPHE(values$flow.frames[[1]],values$manu.annot,input$annotName)
    showNotification("Enrich Manual Annotation", type="message")
    getParams()
  })

##### AUTO ANNOTATION ############################

  observeEvent(input$ref_file,{
    progress <- Progress$new()
    progress$set(message="Read reference file", value=0)
    paths <- input$ref_file$datapath
    names <- as.vector(input$ref_file$name)
    values$ref.flow.frames <- lapply(c(1:length(names)),function(i){
      read.FCS(paths[i],emptyValue=FALSE)
    })
    names(values$ref.flow.frames) <- names
    progress$close()
  })
  
  output$uiMarkersTransRef <- renderUI({
    if(is.null(values$ref.flow.frames)){return(NULL)}
    labels <- pData(parameters(values$ref.flow.frames[[1]]))[,2]
    params <- pData(parameters(values$ref.flow.frames[[1]]))[,1]
    labels[which(is.na(labels))] <- colnames(values$ref.flow.frames[[1]])[c(which(is.na(labels)))]
    labels[which(labels=="<NA>")] <- colnames(values$ref.flow.frames[[1]])[c(which(labels=="<NA>"))]
    names(params) <- labels
    selectInput("markersTransRef","Ref Markers Trans",choices=params,multiple = TRUE)
  })
  
  observeEvent(input$ref_prepross,{
    if(is.null(values$ref.flow.frames))return(NULL)
    progress <- Progress$new()
    progress$set(message="Preprocess ref...",value=0.5)
    ref.flow.frames <- values$ref.flow.frames
    
    if(input$ref_com && !is.null(ref.flow.frames[[1]]@description[["SPILL"]])){
      ref.flow.frames <- lapply(ref.flow.frames,function(i){
        return(compensate(i,i@description[["SPILL"]]))
      })
    }
    
    if(input$trans_ref != "None"){
     ref.flow.frames <- lapply(ref.flow.frames, function(i){
       if(input$trans_ref == "Logicle"){
         return(logicleTransform(i))
       }
       if(input$trans_ref == "Arcsinh"){
         return(arcsinhTransCIPHE(i, marker=NULL,arg = 5))
       }
     }) 
    }
    
    names(ref.flow.frames) <- names(values$ref.flow.frames)
    values$ref.flow.frames <- ref.flow.frames
    progress$close()
  })
  
  output$listRef <- renderTable({
    if(is.null(values$ref.flow.frames))return(NULL)
    filesname <- names(values$ref.flow.frames)
    return(data.frame(filesname))
  })

##### MANUAL GATING ##############################
  
  output$selectGateFile <- renderUI({
    if(is.null(values$flow.frames)){return(NULL)}
    selectInput("gateFile","Select File",choices=names(values$flow.frames),multiple=F)
  })

  output$selectGateStep <- renderUI({
    if(is.null(values$step.gate)){return(NULL)}
    # l <- unlist(lapply(values$setp.gate, function(j){return(names(j))}))
    selectInput("gateStep","Select Step",choices=names(values$step.gate),multiple=F)
  })

  output$selectXGate <- renderUI({
    if(is.null(values$params))return(NULL)
    selectInput("xGate","X",choices=values$params, multiple=F,selected=values$params[[1]])
  })

  output$selectYGate <- renderUI({
    if(is.null(values$params))return(NULL)
    selectInput("yGate","Y",choices=values$params, multiple=F,selected=values$params[[2]])
  })

  observe({
    if(input$tabs != "gate"){return(NULL)}
    if(length(values$step.ff)>0){return(NULL)}
    if(is.null(input$gateFile)){return(NULL)}
    values$step.ff <- list("Root" = values$flow.frames[[input$gateFile]])
    values$step.gate <- list("Root" = NULL)
  })

  output$plotGate <- renderPlot({
    if(length(values$step.ff)==0)return(NULL)
    if(is.null(input$gateStep)) return(NULL)
    plotDens(values$step.ff[[input$gateStep]],c(input$xGate,input$yGate),main="Roots")
    if(!is.null(values$polygon.temp)) {
      points(values$polygon.temp[,1], values$polygon.temp[,2],pch=22, lty=2, lwd=2,col="red")
      lines(values$polygon.temp,lty=2, lwd=2,col="red")
    }
    if(!is.null(values$step.gate)){
      for(i in values$step.gate){
        if(!is.null(i)){
          if(i[[2]]==input$gateStep){lines(i[[1]],lty=1, lwd=2,col="red")}
        }
      }
    }
  })

  observeEvent(input$point,{
    if(is.null(values$polygon.temp)){
      values$polygon.temp <- data.frame(input$point$x,input$point$y)
    } else {
      values$polygon.temp <- rbind(values$polygon.temp,c(input$point$x,input$point$y))
    }
  })
  
  observe({
    updateSelectInput(session,"gateStep",selected = names(values$step.gate)[length(values$step.gate)],choices=names(values$step.gate))
  })

  observeEvent(input$refresh,{
    if(input$tabs != "gate"){return(NULL)}
    if(is.null(input$gateFile)){return(NULL)}
    print("refresh")
    values$step.ff <- list("Root" = values$flow.frames[[input$gateFile]])
    values$step.gate <- list("Root" = NULL)
  })

  output$tableGateManual <- renderDataTable({
    d <- lapply(c(1:length(values$step.ff)), function(i){
      return(c(names(values$step.ff)[i],dim(values$step.ff[[i]])[1],round(dim(values$step.ff[[i]])[1]/dim(values$step.ff[[1]])[1]*100,2)))
    })
    mat <- do.call(rbind, d)
    colnames(mat) <- c("Gate","#Events","%Percentile Total")
    return(data.frame(mat))
  })

  observeEvent(input$createGate,{
    if(is.null(input$stepName) || length(input$stepName)==0) return(NULL)
    if(length(values$step.ff)==0)return(NULL)
    step.ff <- values$step.ff
    step.gate <- values$step.gate
    
    polygon <- values$polygon.temp
    polygon <- rbind(polygon,polygon[1,])
    colnames(polygon) <- c(input$xGate,input$yGate)

    fcs <- values$step.ff[[input$gateStep]]
    new.fcs <- Subset(fcs, polygonGate(gate=polygon))
    
    step.gate[[input$stepName]] <- list(polygon,input$gateStep)
    step.ff <- c(step.ff, list(new.fcs))
    names(step.ff)<- c(names(values$step.gate),input$stepName)

    values$step.ff <- step.ff
    values$step.gate <- step.gate
    names(values$step.ff) <- names(step.ff)
    names(values$step.gate) <- names(step.gate)

    values$polygon.temp <- NULL
  })
  
  observe({
    if(is.null(values$step.gate))return(NULL)
    if(is.null(values$flow.frames)) return(NULL)
    if(input$tabs != "gateview")return(NULL)
    progress <- Progress$new()
    progress$set(message="Gating in progress",value=0.5)
    gate.strat <- c()
    final <- list()
    for(i in values$flow.frames){
      gs <- list("Root"=i)
      for(j in c(1:length(values$step.gate))){
        if(!is.null(values$step.gate[[j]])){
          fcs <- gs[[values$step.gate[[j]][[2]]]]
          gate <- values$step.gate[[j]][[1]]
          new.fcs <<- Subset(fcs, polygonGate(gate=gate))
          gs[[names(values$step.gate)[j]]] <- new.fcs
        }
      }
      final <- c(final, list(gs))
    }
    values$gate.strat <- final
    progress$close()
  })
  
  output$selectGateStepView <- renderUI({
    if(is.null(values$gate.strat)) return(NULL)
    if(is.null(values$step.gate)){return(NULL)}
    selectInput("gateStepView","Select Step",choices=names(values$step.gate),multiple=TRUE)
  })
  
  observeEvent(input$dllSubPop,{
    if(is.null(values$gate.strat)) return(NULL)
    flow.frames <- c()
    names <- c()
    for(i in input$gateStepView){
      for(j in c(1:length(values$flow.frames))){
        name <- paste0(values$names.files[[j]],"_",i,".fcs")
        fcs <- values$gate.strat[[j]][[grep(i, names(values$gate.strat[[j]]))]]
        flow.frames <- c(flow.frames, list(fcs))
        print(name)
        names <- c(names,list(name))
      }
    }
    ## CREATE ZIP HERE
    output$linkBigSubPop <- renderText({
      return("Download")
    })
  })
  
  observeEvent(input$loadSubPop,{
    if(is.null(values$gate.strat)) return(NULL)
    if(is.null(input$gateStepView)) return(NULL)
    flow.frames <- c()
    names <- c()
    for(i in input$gateStepView){
      for(j in c(1:length(values$flow.frames))){
        name <- paste0(values$names.files[[j]],"_",i,".fcs")
        fcs <- values$gate.strat[[j]][[grep(i, names(values$gate.strat[[j]]))]]
        flow.frames <- c(flow.frames, list(fcs))
        print(name)
        names <- c(names,list(name))
      }
    }
    values$names.files <- unlist(names)
    values$flow.frames <- flow.frames
    names(values$flow.frames) <- values$names.files
    showNotification(ui="flowFrames replace by subPop", type = "message")
  })
  
  output$gateViewStep <- renderUI({
    if(is.null(values$gate.strat)) return(NULL)
    progress <- Progress$new()
    progress$set(message="Plot in progress...",value=0.9)
    # browser()
    outputGateView <- lapply(c(1:length(values$flow.frames)),function(ligne){
      renderPlot({
        par(mfrow=c(1,(length(values$gate.strat[[ligne]])-1)))
        parent <- unique(unlist(lapply(values$step.gate, function(a){return(a[[2]])})))
        
        for(g in parent){
          plotDens(values$gate.strat[[ligne]][[g]],channels=colnames(values$step.gate[[g+1]][[1]]),
                   main=g)
          for(i in values$step.gate){
            if(!is.null(i)){
              if(i[[2]]==names(values$gate.strat[[ligne]])[g]){
                lines(i[[1]],lty=1, lwd=2,col="red")
              }
            }
          }
        }
      },width=300*(length(values$gate.strat[[ligne]])-1),height = 300)
    })
    do.call(tagList,outputGateView)
    progress$close()
    return(outputGateView)
  })

}