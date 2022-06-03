
logiclTransformCIPHE <- function(flow.frame, value = NULL, marker = NULL){
  
  if(is.null(marker)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- marker
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value) || is.na(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])){
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      )
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x)
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  
  w.values <- (4.5-log10(262143/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  for(t in 1:length(markers.transform)){
    lgcl <- logicleTransform(w=w.values[t])
    flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
  }
  
  return(flow.frame)
}

inversLogiclTransformCIPHE <- function(flow.frame, value = NULL, marker = NULL){
  if(is.null(marker)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- marker
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value) || is.na(value)){
    if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) 
        as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
      ) 
    } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]])) {
      r.values <- unlist(lapply(list.index, function(x) 
        as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
      )   
    } else {
      r.values <- rep(90, length(list.index))
    }
  }
  else {
    r.values <- rep(value, length(list.index))
  }
  
  w.values <- (4.5-log10(262144/abs(r.values)))/2
  w.values[which(w.values<0)] <- 0.5
  w.values[which(is.infinite(w.values))] <- 0.5
  
  flow.frame.inv <- flow.frame
  
  for(t in 1:length(markers.transform)){
    invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
    flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
  }
  
  return(flow.frame.inv)
}

arcsinhTransCIPHE <- function(flow.frame, marker=NULL, arg=5){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  if(is.null(marker) || length(marker)<1){
    marker <- colnames(flow.frame)
  }
  # print(marker)
  mat <- mat[,marker]
  colnames(mat) <- marker
  mat <- asinh(mat/arg)
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

inversArcsinhTransCIPHE <- function(flow.frame, marker=NULL, arg){
  raw <- flow.frame@exprs
  mat <- flow.frame@exprs
  if(is.null(marker) || length(marker)<1){
    marker <- colnames(flow.frame)
  }
  # print(marker)
  mat <- mat[,marker]
  colnames(mat) <- marker
  mat <- sinh(mat)*arg
  raw[,marker] <- mat[,marker]
  flow.frame@exprs <- raw
  return(flow.frame)
}

flowVS.CIPHE <- function(flow.frames, markers=NULL){
  if(is.null(markers)) return(flow.frames)
  fs <- flowSet(flow.frames)
  png(paste0(tempdir(),"/temp.png"))
  res <- flowVS::estParamFlowVS(fs, channels=markers)
  dev.off()
  file.remove(paste0(tempdir(),"/temp.png"))
  flow.frames <- lapply(flow.frames, function(fcs){
    for(i in c(1:length(markers))){
      m <- markers[i]
      fcs@exprs[,m] <- asinh(fcs@exprs[,m]/res[i])
    }
    return(fcs)
  })
  return(list(flow.frames,res))
} 

deCompensateFlowFrame <- function(x, spillover){
  if(!is.null(spillover)){
    cols <- colnames(spillover)
    sel <- cols %in% colnames(x)
    if(!all(sel)) {
      stop(keyword(x)[["FILENAME"]], "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
           paste(cols[!sel], collapse=", "), call.=FALSE)
    }
    e <- exprs(x)
    e[, cols] <- e[, cols] %*% spillover
    exprs(x) = e
    return(x)
  } else {
    return(x)
  }
}

get_cluster_groups_table <- function(v, key){
  tags$table(class = "table table-hover table-striped",
             tags$tr(tags$td(
               v[1],
               tags$button(class = "btn btn-xs btn-warning pull-right", onClick = sprintf("Shiny.onInputChange('clusteringui_remove_clustering_group', {'key':'%s', 'x':Math.random()})", key),
                           tags$span(class = "glyphicon glyphicon-trash")
               )
             )),
             ifelse(length(v > 1),
                    tagList(lapply(tail(v, n = -1), function(x) {tags$tr(tags$td(x))})),
                    tagList()
             )
  )
}

renameLabeling <- function(flow.frame, label){
  
  for(i in 1:as.integer(flow.frame@description$'$PAR'))
  {
    index.name <- paste0("$P",i,"N") # Name
    index.label <- paste0("$P",i,"S") # Label
    par.name <- flow.frame@description[[index.name]]
    par.label <- flow.frame@description[[index.label]]
    new.par.name <- as.character(label[which(label[,"Parameter"]==par.name),"Name"])
    new.par.label <- as.character(label[which(label[,"Parameter"]==par.name),"Label"])
    
    # change name
    if (((length(new.par.name) == 0) && (typeof(new.par.name) == "character")) || (new.par.name == ""))
    {
      new.par.name <- par.name
      if(length(new.par.name)==0) new.par.name <- ""
    }
    if((!is.null(flow.frame@description[["SPILL"]])) && (class(flow.frame@description[['SPILL']]) == "matrix"))
    {
      colnames(flow.frame@description[["SPILL"]])[which(colnames(flow.frame@description[["SPILL"]])==par.name)] <- new.par.name
    }
    if((!is.null(flow.frame@description[["$TR"]])) && (class(flow.frame@description[["$TR"]])=="character"))
    {
      flow.frame@description[["$TR"]] <- gsub(par.name, new.par.name, flow.frame@description[["$TR"]], fixed=TRUE)
    }
    flow.frame@description[[index.name]] <- as.character(new.par.name)
    colnames(flow.frame)[which(colnames(flow.frame)==par.name)] <- as.character(new.par.name)
    
    # change label
    if(((length(new.par.label) == 0) && (typeof(new.par.label) == "character")) || (new.par.label == "")) 
    {
      new.par.label <- par.label
      if(length(new.par.label)==0) new.par.label <- ""
    }
    flow.frame@description[[index.label]] <- as.character(new.par.label)
    
    # change in parameter
    pData(flow.frame@parameters)[i,1] <- as.character(new.par.name)
    pData(flow.frame@parameters)[i,2] <- as.character(new.par.label)
    flow.frame@parameters@data[i,1] <- as.character(new.par.name)
    flow.frame@parameters@data[i,2] <- as.character(new.par.label)
  }
  
  return(flow.frame)
}

enrich.FCS.CIPHE <- function(original, new.column, nw.names=NULL){
  new_p <- parameters(original)[1,]
  
  ## Now, let's change it's name from $P1 to $P26 (or whatever the next new number is)
  new_p_number <- as.integer(dim(original)[2]+1)
  rownames(new_p) <- c(paste0("$P", new_p_number))
  
  ## Now, let's combine the original parameter with the new parameter 
  library('BiocGenerics') ## for the combine function
  allPars <-  BiocGenerics::combine(parameters(original), new_p)
  
  ## Fix the name and description of the newly added parameter, say we want to be calling it cluster_id
  
  if(is.null(nw.names)){
    new_p_name <- "cluster"
  } else {
    new_p_name <- nw.names
  }
  
  allPars@data$name[new_p_number] <- new_p_name
  allPars@data$desc[new_p_number] <- new_p_name
  
  new_exprs <- cbind(original@exprs, new.column)
  colnames(new_exprs) <- c(colnames(original@exprs),new_p_name)
  
  new_kw <- original@description
  new_kw["$PAR"] <- as.character(new_p_number)
  new_kw[paste0("$P",as.character(new_p_number),"N")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"S")] <- new_p_name
  new_kw[paste0("$P",as.character(new_p_number),"E")] <- "0,0"
  new_kw[paste0("$P",as.character(new_p_number),"G")] <- "1"
  new_kw[paste0("$P",as.character(new_p_number),"B")] <- new_kw["$P1B"]
  new_kw[paste0("$P",as.character(new_p_number),"R")] <- new_kw["$P1R"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmin")] <- new_kw["flowCore_$P1Rmin"]
  new_kw[paste0("flowCore_$P",as.character(new_p_number),"Rmax")] <- new_kw["flowCore_$P1Rmax"]
  
  ## Now, let's just combine it into a new flowFrame
  new_fcs <- new("flowFrame", exprs=new_exprs, parameters=allPars, description=new_kw)
  
  return(new_fcs)
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

concatenateCIPHE <- function(flow.frames, params="Flag") {
  ff.concat <- NULL
  n <- length(flow.frames)
  for(i in 1:n){
    ff.raw <- flow.frames[[i]]
    p <- matrix(i, nrow = nrow(ff.raw), ncol=1, dimnames = list(NULL, params))
    new.col <- as.vector(p)
    ff.raw <- enrich.FCS.CIPHE(ff.raw, new.col, nw.names=params)
    if(is.null(ff.concat)){
      ff.concat  <- ff.raw
    } else {
      exprs(ff.concat) <- rbind(exprs(ff.concat),exprs(ff.raw))
    }
  }
  return(ff.concat)
}

createFCSfromCSV <- function(csv){
  print("a")
  id <- which(unlist(lapply(colnames(csv),nchar))==0)
  if(length(id)>0){
    csv <- csv[,-id] #delete cols with no name
  }
  id <- which(unlist(sapply(csv,function(x){any(is.character(x))}))==TRUE)
  if(length(id)>0){
    csv <- csv[,-id]#delete cols with string
  }
  row.names(csv) <- as.numeric(c(1:dim(csv)[1]))
  print("b")
  # csv[which(is.na(csv))] <- 0
  cols <- colnames(csv)
  csv <- data.matrix(csv)
  print("c")
  # csv <- str_replace(csv,pattern = "\\,",replacement = "\\.")
  # colnames(csv) <- cols
  metadata <- data.frame(name=colnames(csv),desc=colnames(csv),stringsAsFactors = F)
  metadata$range <- apply(apply(csv,2,range),2,diff)
  metadata$minRange <- apply(csv,2,min)
  metadata$maxRange <- apply(csv,2,max)
  row.names(metadata) <- paste0("$P",row.names(metadata))
  
  data.ff <- new("flowFrame",exprs=csv,parameters=AnnotatedDataFrame(metadata))
  data.ff <- enrich.FCS.CIPHE(data.ff, c(1:dim(data.ff)[1]), "index")

  return(data.ff) 
}

updateKeywords <- function(flowFrame){
  flowFrame@exprs[which(is.na(flowFrame@exprs))] <- 0
  params = parameters(flowFrame)
  pdata = pData(params)
  flowFrame@exprs[which(is.na(exprs(flowFrame)))] <- 0
  max <- round(max(exprs(flowFrame)))
  if(max > 1000000000) max <- 1000000000
  for (i in 1:ncol(flowFrame)){
    s = paste("$P",i,"S",sep="");
    n = paste("$P",i,"N",sep="");
    r = paste("$P",i,"R",sep="");
    b = paste("$P",i,"B",sep="");
    e = paste("$P",i,"E",sep="");
    keyval=list();
    keyval[[s]] = colnames(flowFrame)[i];
    keyval[[n]] = colnames(flowFrame)[i];
    keyval[[r]] = as.character(max) #ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[b]] = 32;
    keyval[[e]] = "0,0";
    keyword(flowFrame) = keyval;
    
    pdata[i,"minRange"]=min(exprs(flowFrame)[,i])
    pdata[i,"maxRange"]=max(exprs(flowFrame)[,i])
    pdata[i,"range"]=max(exprs(flowFrame)[,i])
    # colnames(pdata)[i] <- paste("$P",i)
  }
  pData(params)=pdata
  parameters(flowFrame)=params

  return(flowFrame)
}

FlowRepositoryReadFCS <- function(data, i=c(1:length(data@fcs.files))){
  # data is a flowRepData
  # i is index of file from this repositorydataset you want ddl,
  # by default you download all fcs files from this repository
  # return a list opf flow.frames with name
  
  res <- lapply(i, function(j){ #download data in temp file
    return(download(fcs.files(data)[[j]],tempdir()))
  })
  
  flow.frames <- lapply(res, function(j){ #read fcs from temp file
    return(read.FCS(localpath(j)))
  })
  
  names <- lapply(res, function(j){ #get name from FlowRepository
    return(j@name)
  })
  
  names(flow.frames) <- names #name your list with files names
  
  lapply(res, function(j){ # clear temp file
    file.remove(localpath(j))
  })
  
  return(flow.frames)
}

plotDensCIPHE <- function(x,y,z=0,xlab="x",ylab="y",main="main"){
  if(length(x)>10000){
    pch <- "."
  } else {
    pch <- 20
  }
  if(length(z)==length(x)){
    palette <- colorRampPalette(c(rgb(0,0,1,0.3),rgb(1,1,0,0.3),rgb(1,0,0,0.3)),alpha=TRUE)
    col <- palette(20)[as.numeric(cut(z,breaks=20))]
    plot(x,y,color=cl,pch=pch,xlab=xlab, ylab=ylab, main=main)
  } else {
    if(length(x)>10000){
      plotDens(fcs, c())
    } else {
      
    }
  }  
}

read.FCS.CIPHE <- function(fcs){
  out <- tryCatch({
    read.FCS(fcs,emptyValue=FALSE)
  },
  error = function(cond){
    return(NULL)
  })
}

catch.createFCSfromCSV <- function(csv){
  if(dim(csv)[2]<2){return(NULL)}
  out <- tryCatch({
    createFCSfromCSV(csv)
  },
  error = function(cond){
    return(NULL)
  })
}

normalize.FCS.CIPHE <- function(fcs,markers=NULL,quantile){
  if(is.null(markers)){
    if(is.null(fcs@description[["SPILL"]])){
      markers <- colnames(markers)
    } else {
      markers <- colnames(fcs@description[["SPILL"]])
    }
  }
  
  for(i in markers){
    x <- fcs@exprs[,i]
    x <- x/quantile(x,probs = quantile)
    fcs@exprs[,i] <- x
  }
  
  return(fcs)
}

detach_package <- function(pkg, character.only = FALSE){
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

chooseDir <-  function(){
    OS <- Sys.info()["sysname"]
      if (OS=="Windows") {
      Dir <-
        choose.dir(default = "", caption = "Select a Folder for saving:")
      }
      else if (OS=="Linux") {
        Dir <- tk_choose.dir(default = "", caption = "Select a Folder for saving:")
        }
      else {
        Dir <- choose.mac.dir()
        }
    return(Dir)
}