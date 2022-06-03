## This file contain all function who needed to gating automate with flowDensity
# in different panel development 


## Load library
######################################################################################################
######################################################################################################
######################################################################################################

# Cytometry packages
suppressWarnings(suppressMessages(library("flowCore")))
suppressWarnings(suppressMessages(library("flowDensity")))
suppressWarnings(suppressMessages(library("flowType")))
suppressWarnings(suppressMessages(library("RchyOptimyx")))
# suppressWarnings(suppressMessages(library("flowMeans")))
# suppressWarnings(suppressMessages(library("flowMerge")))
suppressWarnings(suppressMessages(library("flowClust")))
suppressWarnings(suppressMessages(library("openCyto")))

# Tools packages
suppressWarnings(suppressMessages(library("stringr")))
suppressWarnings(suppressMessages(library("gtools")))
# suppressWarnings(suppressMessages(library("XML")))

# Plot packages
# suppressWarnings(suppressMessages(library("gplots")))
# suppressWarnings(suppressMessages(library("png")))
# suppressWarnings(suppressMessages(library("animation")))

# Cluster packages
suppressWarnings(suppressMessages(library("foreach")))
suppressWarnings(suppressMessages(library("doParallel")))

## Function
######################################################################################################
######################################################################################################
######################################################################################################



renameLabel <- function(flow.frame, csv)
{
	
	######################################################################################################
	#
	#	Rename the label of the FCS objects with the label table. 
	#
	#	Args : 
	#					flow.frame : flowFrame value is an FCS object
	#					label.table : label.table object loaded previously with loadMarker
	#
	#	Returns :
	#					flow.frame : fcs object with the new parameters name
	#
	######################################################################################################

	label <- read.csv(csv, check.names=FALSE)

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





getPeaks <- function(frame,chans,tinypeak.removal=tinypeak.removal)
{

  ######################################################################################################
  #
  # Finds the peaks in the given density
  #
  # Args :
  #		   		dens: density of the channel whose peaks are going to be found. It is a variable of class 'density'.
  # 		  	w: the length of the window where the function searches for the peaks. If all peaks required, use the default w=1.
  # Returns:
  #   			peaks in the density of the provided channel
  #
  ######################################################################################################

 	data <- exprs(frame)[,chans]
  dens <- density(data[which(!is.na(data))])
  dens <- smooth.spline(dens$x, dens$y, spar=0.4)
  dens$y[which(dens$y<0)] <- 0
  d <- dens$y
  w <- 1
  peaks <- c()
  peaks.ind <- c()

  for(i in 1:(length(d)-w))
  {
    if(d[i+w] > d[(i+w+1):(i+2*w)] && d[i+w] > d[i:(i+w-1)] && d[i+w] > tinypeak.removal*max(d)){ # also removes tiny artificial peaks less than ~%4 of the max peak
      peaks <- c(peaks, dens$x[i+w])
      peaks.ind <- c(peaks.ind, i+w)
    }
  }

  return(list(Peaks=peaks, Dens=dens,Ind=peaks.ind))
}



loadMarkers <- function(label.file, flow.frame, rename = FALSE, rewrite.label = FALSE)
{
	
	#######################################################################################################
	#
  # Load a label file to create a list of marker for FCS file and, if its asked, 
  #	this functiob can call the relabel function in shiny interface 
  #
  # Args : 
  #				 labelFile : a CSV or TXT file format labelFCSChannels 
  #										 3 columns with Parameters, Relabel, and name.
  #				 rename : if TRUE, use Label name if FALSE use Parameter name.
  #				 flow.frame : flowFrame value is an FCS object load by flowCore used to validate
  #											the marker found in the labelFile, return error if its not
  #
  # Returns:
  #					list.marker : a list of variables with for each marker,
  #												each key are the CD value used like marker$CDX in your gating
  #
  #######################################################################################################

  label <- read.csv(label.file, check.names=FALSE) # load csv file

  no.marker <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time") # universal parameter

  list.marker <- list() # create the list output

  markers <- colnames(flow.frame)[which(colnames(flow.frame)%in%no.marker==FALSE)] # load all marker execpted 

	for(i in 1:length(label[,2])) # assigne for each marker the marke.cdx <- 
	{
		m <- as.vector(label[i,2])
		if(m != "" && length(m)>0 && !is.na(m) && !is.null(m))
		{
			# list.marker[[m]] <- as.vector(label[i,3]) # append to list.marker the key avec value for each marker CDX
			list.marker[[m]] <- colnames(flow.frame)[which(m==pData(flow.frame@parameters)[,2])]
		}
	}

  return(list.marker)
}





qualityGate <- function(flow.frame, all.marker, live.marker, filter=FALSE, first.gate = TRUE, filter.size = TRUE, compensate.value=TRUE, transform.value = TRUE) 
{

	#######################################################################################################
	#
	#	First step for all gating pipeline, compensate, transformation, cleaning and extract
	# good cells from all, clog, live, and size gate
	#
	# Args :
	#					flow.frame : FCS objects flow.frame S4
	#					all.marker : the marker used to see All events with Time
	#					live.marker : the marker to found the live cells 
	#
	# Returns :
	#					flow.frame.output : a list with all flow.frame object and gate or plotting
	#
	#######################################################################################################
		
	filter.marker <- c(grep (colnames(flow.frame),pattern = "FSC*"),grep (colnames(flow.frame),pattern = "SSC*"))

	for(chan in filter.marker) # filter for physical marker before trans and comp with max threshold
  {
    stain.max <-max(exprs(flow.frame)[,chan])
    margins <- which(exprs(flow.frame)[, chan] >= stain.max)
    exprs(flow.frame) <- exprs(flow.frame)[-margins,]
  }

  if(compensate.value==TRUE){
		flow.frame <- compensate(flow.frame, flow.frame@description$SPILL) # compensation
	}
	no.transform <- c(colnames(flow.frame)[filter.marker],"Time","TIME") # select no transform marker
	markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE] # select transform markers
	#  lgcl <- estimateLogicle(flow.frame, markers.transform) # compute the logicle value
	if(transform.value==TRUE){
		flow.frame <- logiclTransformCiphe(flow.frame) # transform value
	}
	# lgcl <- logicleTransform(w = 0.5, t = 10000, m = 2.5)
 	# trans <- transformList(markers.transform, lgcl)
 	# flow.frame <-  transform(flow.frame,trans)

	for(chan in markers.transform) # filter under -0.5 for each marker after compensation and transformation
	{
		margins <- which(exprs(flow.frame)[,chan]< -1)
		if(length(margins)>0)
		{
			exprs(flow.frame) <- exprs(flow.frame)[-margins,]
		}
	}

	# for(chan in markers.transform) # filter under -0.5 for each marker after compensation and transformation
	# {
	# 	margins <- which(exprs(flow.frame)[,chan]> 5)
	# 	if(length(margins)>0)
	# 	{
	# 		exprs(flow.frame) <- exprs(flow.frame)[-margins,]
	# 	}
	# }

	if(first.gate == FALSE)
	{
		return(flow.frame)
	}

	flow.frame.output <- c(flow.frame) # initiate list output with flow.frame and cellpopulation objects

	flowType <- flowType(flow.frame, 
									c("FSC-A",live.marker), 
									Methods = 'kmeans',
									PartitionsPerMarker = c(2,2)
									)					

	live <- flowDensity(flow.frame, 
											c("FSC-A",live.marker), 
											position=c(TRUE,FALSE), 
											gate = c(flowType@Thresholds[[1]], flowType@Thresholds[[2]]),
											#gate = c(50000,3),
											ellip.gate = TRUE,
											scale=0.999
											)						
	fs.live <- getflowFrame(live)
	flow.frame.output <- c(flow.frame.output,live, fs.live)

  if(length(getPeaks(fs.live, "FSC-A", 0.3)$Peaks)>1 && filter.size == TRUE)
	{
		x.th <- deGate(fs.live, "FSC-A")
		size <- flowDensity(fs.live,
												c("FSC-A","SSC-A"),
												position=c(TRUE,TRUE),
												gates = c(x.th, NA),
												use.percentile = c(FALSE,TRUE),
												percentile = c(NA,0.001)
												)						
		fs.size <- getflowFrame(size)	
	} else {
		x.th <- 50000
		size <- flowDensity(fs.live,
												c("FSC-A","SSC-A"),
												position=c(TRUE,TRUE),
												gates = c(x.th, NA),
												use.percentile = c(FALSE,TRUE),
												percentile = c(NA,0.001)
												)						
		fs.size <- getflowFrame(size)	
	}
	flow.frame.output <- c(flow.frame.output,size, fs.size)

	# singletfsc <- flowDensity(fs.size,
 # 														c("FSC-A","FSC-H"),
 # 														use.percentile = c(TRUE,TRUE),
 #                            position = c(FALSE,FALSE),
 #                            percentile = c(0.9,0.9),
 #                            ellip.gate=TRUE,
 #                            scale=.999999999999999
 #                           )

	singletfsc <- flowDensity(fs.size,
    channels=c("FSC-A","FSC-H"),
    use.percentile = c(TRUE,TRUE),
    position = c(FALSE,FALSE),
    percentile = c(0.9,0.9),
    ellip.gate=TRUE,
    scale=.999999999999999
  )
  p <- openCyto:::.singletGate(fs.size, channels = c("FSC-A","FSC-H"))
  p@boundaries <- rbind(p@boundaries,p@boundaries[1,])
  singletfsc@filter <- p@boundaries
  fs.singletfsc<-applyGateObject(fs.size, singletfsc)

 	fs.singletfsc <- getflowFrame(singletfsc)
 	flow.frame.output <- c(flow.frame.output, singletfsc, fs.singletfsc)

 	singletssc <- flowDensity(fs.singletfsc,
 													c("SSC-W","SSC-A"),
 													position=c(FALSE,FALSE),
 													use.percentile=c(TRUE,TRUE),
 													percentile = c(0.95,0.9999999999)
 													)
 	fs.singletssc <- getflowFrame(singletssc)
 	flow.frame.output <- c(flow.frame.output,singletssc,fs.singletssc)

 	return(flow.frame.output)
}





clusteriseGating <- function(experiment.name, gating.panel, group, result=NULL) 
{

	#######################################################################################################
	#
	#	Core function used to clusterise the gating in multiple core
	#
	# Args : 
	#					experiment.name : The name of experiment, the name of the file with much information
	#					gating.panel : the name of the fonction used to make the gating
	#
	# Returns : 
	#					Create PDF, PNG and CSV file output but they are no object return
	#
	#######################################################################################################

	# create  path
	path <- paste0("/media/data/html/INPUT/DATA/",experiment.name)
	output.path <- paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/")

	# check if exist and delete
	list <- list.files("/media/data/html/OUTPUT/GATING", recursive=TRUE, full.names=TRUE)
	list.experiment <- list[grep(experiment.name, list)]
	unlink(list.experiment)

	# create repository
	suppressWarnings(suppressMessages(dir.create(output.path)))

	#load variable
	mat <- NULL
	row.names <- c()
	mice.table <- read.csv(paste0("/media/data/html/INPUT/METADATA/",experiment.name,"_MetaData.csv"), check.names = FALSE) 

	# load fcs files 
	list.files <- list.files(path, pattern=".fcs",recursive = TRUE, full.names=TRUE)
	list.files <- list.files[mixedorder(list.files)]


	## Gating sur les CRL
	if(group == "CTRL"){

		ctrl.index <- which(mice.table[,"GROUP"]=="CTRL")
		ctrl.list.files <- list.files[ctrl.index]

		# create cluster and cluterise
		nb.file <- length(ctrl.list.files)
		result <- c()
		cl <- makeCluster(nb.file)
		registerDoParallel(cl)
		print(paste0("Number of cluster - ",nb.file,"</br>"))

		result <- list()
		t1 <- Sys.time()
		result <- foreach(i=1:nb.file) %dopar% {
			gating.panel(ctrl.list.files[i], output.path, ctrl.index[i])
		}
		print(paste0(Sys.time()-t1,"</br>"))
		stopCluster(cl)
		return(result)
	}

	seuil <- result

	seuil <- matrix(unlist(seuil), ncol=length(seuil), byrow=FALSE)
	seuil <- apply(seuil, 1, function(x)  mean(x, na.rm=TRUE)) # median

	nb.file <- length(list.files)
	result <- c()
	cl <- makeCluster(nb.file)
	registerDoParallel(cl)
	print(paste0("Number of cluster - ",nb.file,"</br>"))

	## Gating Mutant
	t1 <- Sys.time()
	resultPartTwo <- foreach(i=1:nb.file) %dopar% {
		gating.panel(list.files[i], output.path, i, seuil)
	}

	print(paste0(Sys.time()-t1,"</br>"))
	stopCluster(cl)

	resultTotal <- c(resultPartTwo)

	# concat all result in matrix
	mat <- matrix(unlist(resultTotal), nrow=length(list.files), byrow=TRUE)
	colnames(mat) <- colnames(resultTotal[[1]])
	row.names(mat) <- unlist(lapply(resultTotal, FUN = function(x) { return(rownames(x))}))

	# mat <- mat[,c(1:dim(mat)[2])[-seq(3,dim(mat)[2],3)]] # A désactiver pour CYTOP01A EXT099
	# mat <- mat[,-c(grep("Parent",colnames(mat)))]

	mat <- cbind(mice.table,mat[c(match(rownames(mat),row.names(mice.table))),]) #row.names(mice.table)[1]
	suppressWarnings(suppressMessages(write.csv(mat, paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating.csv"), quote=FALSE,sep="\t")))

	pngfile <- list.files(path=output.path, pattern=".png", ,recursive=TRUE, full.names = TRUE)
	pngfile <- pngfile[mixedorder(pngfile)]
	pdf(paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating.pdf"), width = 12, height = 9)
	for(i in pngfile){
		plot.new()
		img <- readPNG(i)
		rasterImage(img,0,0,1,1)
	}
	dev.off()
	print(paste0(Sys.time()-t1,"</br>"))
	
}





clusteriseCombi <- function(experiment.name, combi.script, disco.script)
{

	#######################################################################################################
	#
	#	Core function used to clusterise the gating in multiple core
	#
	# Args : 
	#					experiment.name : The name of experiment, the name of the file with much information
	#
	# Returns : 
	#					Creates files with fcs invers, table with count and freq, plot all histograme with thresholds
	#
	#######################################################################################################

	## Meta Data
	meta.data <- read.csv(paste0("/media/data/html/INPUT/METADATA/",experiment.name,"_MetaData.csv"), check.names = FALSE, row.names = 1) 

	list.popFCS <- list.files(paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/FILE/"))
	list.popFCS <- list.popFCS[mixedorder(list.popFCS)]
	popFCS <- unlist(lapply(list.popFCS, FUN = function(x) {return(str_split(x,"_")[[1]][1])}))
	popFCS <- unique(popFCS)

	xml.file <- str_split(experiment.name,"_")[[1]][3]
	xmlPath <- paste0("/media/data/html/INPUT/TEMPLATE/",xml.file,".xml")

	popTable <- read.csv("/media/data/html/INPUT/DATABASE/populationName.csv", check.names = FALSE)
	indexPop <-  which(xml.file==popTable[,"PANEL"])
	indexPop <- indexPop[as.vector(unlist(lapply(popFCS, function(x) {which(x==popTable[indexPop,"VALUE"])})))] #which(popTable[,"VALUE"]%in%popFCS)
	popFCS2 <- as.vector(popTable[indexPop,"DIVA"])
	popFCS3 <- as.vector(popTable[indexPop,"PLOT"])

	unlink(paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/"), recursive=TRUE)

	list.fcs <- list.files(paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/FILE"), recursive=TRUE,full.names=TRUE, pattern=".fcs")
	list.fcs <- list.fcs[mixedorder(list.fcs)]

	outputPath <- paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/")
	dir.create(outputPath)
	dir.create(paste0(outputPath,"invFILE/"))

	## Create cluster by population
	cluster <- length(popFCS)
	cl <- makeCluster(cluster)
	registerDoParallel(cl)

	result <- c()
	t1 <- Sys.time()

	# Run mutil core
	result <- foreach(i=1:cluster) %dopar% 
	{
	  combi.script(list.fcs[grep(popFCS[i], list.fcs)], 
	                outputPath,
	                popFCS2[i],
									popFCS[i],
									xmlPath,
	                )
	}
	stopCluster(cl)
	print(paste0(Sys.time()-t1,"<br/>"))
	sub.group <- as.vector(meta.data[,"GROUP"])

	## Write CSV first Positive Population
	listCSV <- list.files(outputPath, pattern="Combi", full.names=TRUE, recursive=TRUE)
	listCSV <- listCSV[mixedorder(listCSV)]
	freq.deepOne <- NULL
		
	for(i in listCSV)
	{
		csv <- read.csv(i, row.names = 1)
		temp.indice <- grep("{1}[\\+]",row.names(csv))
		indice <- unlist(lapply(row.names(csv)[temp.indice], function(x){
			if(length(strsplit(x,"[+-]")[[1]])==1){
				return(x)
			}
		}))
			
		cols <- unlist(lapply(colnames(csv), function(x){
			return(strsplit(x,"CYTO")[[1]][2])
		}))
		colnames(csv) <- cols
		
		pop <- str_split(i,"_")[[1]][length(str_split(i,"_")[[1]])]
		pop <- str_split(pop,".csv")[[1]][1]
		
		temp <- csv[indice,]
		row.names(temp) <- paste0(popFCS3[which(pop==popFCS)], indice)
		
		if(is.null(freq.deepOne)) {
	    freq.deepOne <- temp
	  } else {
	    freq.deepOne <- rbind(freq.deepOne,temp)
	  }
	}
	
	## HERE ADD NAME MICE OR NUMBER MICE
	#####################################	#####################################
	colnames(freq.deepOne) <- c(rbind(as.vector(paste0(sub.group,"_",row.names(meta.data),"_POP")), as.vector(paste0(sub.group,"_",row.names(meta.data),"_FREQ"))))

	write.csv(freq.deepOne, paste0(outputPath,"DeepOne.csv"), quote=FALSE)

	print(paste0(Sys.time()-t1,"</br>"))

	pathCombi <-  paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/RSAVE/")
	pathTable <-  paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/")
	pathFCS <- paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/FILE/")
	pathDATA <- paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating.csv")
	outputPath <- paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/discoPheno/")
	outputLog <- paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/pdfLog/")
	outputVolcano <- paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/VolcanoPlot/")

	dir.create(outputPath)
	dir.create(outputLog)
	dir.create(outputVolcano)

	listCombi <- list.files(pathCombi, pattern=".Rdata", recursive = TRUE, full.names = TRUE)
	listTable <- list.files(pathTable, pattern=".csv", recursive = TRUE, full.names = TRUE)
	listFCS <- list.files(pathFCS, pattern=".fcs", recursive = TRUE, full.names = TRUE)
	listFCS <- listFCS[mixedorder(listFCS)]

	group <- unique(sub.group)

	all.group <- lapply(group, function(x){paste0(sub.group,"_",row.names(meta.data))[which(str_detect(paste0(sub.group,"_",row.names(meta.data)),as.vector(x)))]})
	
	names(all.group) <- group

	##### HARD CORE ##########
	sam.group <- rep.int(1, length(sub.group))
	j <- 2
	for(i in group){
		if(i == "CTRL"){
			sam.group[grep(i,sub.group)] <- 1
		} else {
			sam.group[grep(i,sub.group)] <- j
			
			j <- j+1
		}
	}

	##############################

	cluster <- length(popFCS)
	cl <- makeCluster(cluster)
	registerDoParallel(cl)

	t1 <- Sys.time()
	result <- foreach(i=1:cluster) %dopar% 
	{
	  disco.script(listCombi[grep(popFCS[i],listCombi)], 
	                  listTable[grep(popFCS[i],listTable)], 
										listFCS[grep(popFCS[i], listFCS)],
	                  outputPath,
										outputLog,
	                  popFCS[i], 
	                  TRUE, 
	                  0.5,
										all.group,
										sam.group,
										pathDATA,
										popFCS2[i],
										popFCS3[i],
										pathTable,
										outputVolcano
										)
	}
	stopCluster(cl)
	t2 <- Sys.time()-t1
	print(t2)

	################################# WRITE OUTPUT PDF ############################################

	pngfile <- list.files(path=outputPath, pattern=".png", recursive=TRUE, full.names = TRUE)
	pngfile <- pngfile[mixedorder(pngfile)]
	gatingfile <- list.files(path=paste0("/media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/"), pattern=".png",recursive=TRUE, full.names = TRUE)
	gatingfile <- gatingfile[mixedorder(gatingfile)]
	g <- length(all.group[[1]])-1
	histfie <- list.files
	mutants <- names(all.group)[which(names(all.group)!="CTRL")]

	t1 <- Sys.time()
	for(i in mutants){
		pdf(paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/",experiment.name,"_",i,"_discoPheno-V1.pdf"), width = 12, height = 9)
		
		plot(0:10, typ="n",xaxt="n",yaxt="n",bty="n",xlab="", ylab="")
		text(6,8, paste0("Experiment : ",experiment.name),cex=1.5)
		text(6,7, paste0("Population names:"),cex=2)
		text(6,6, paste0(paste(popFCS3, sep=" / ",collapse=" ")),cex=1)
		text(6,5, paste0("Date : ",format(Sys.time(), "%a %b %d %X %Y")),cex=2)
		text(6,4, paste0("CTRL vs ",i), cex=2)

		## Create volcano plot

		volcanoList <- list.files(outputVolcano,pattern=".png", full.name=TRUE, recursive=FALSE)[grep(i,list.files(outputVolcano,pattern=".png", full.name=TRUE, recursive=FALSE))]
		volcanoList <-  volcanoList[unlist(lapply(c(1:length(popFCS3)), function(x){grep(popFCS3[x],volcanoList)}))]
		x <- round(length(volcanoList)/2)
		y <- 3
		par(mar=c(0,0,0,0))
		plot(c(1:x),rep(y,x),ann=FALSE, xlab="", ylab="", type="n", axes=FALSE,xlim=c(-0.5,x+0.5),ylim=c(0.5,3.5))
		for(l in c(1:length(volcanoList))){
			xmax <- l
			ymax <- 3
			if(l>x){
				ymax <- 2
				xmax <- l-x
			}
			xmin <- xmax-1
			ymin <- ymax-1

			img <- readPNG(volcanoList[l])
			rasterImage(img,xmin,ymin,xmax,ymax, interpolate = FALSE)
		}
		
		temp <- pngfile[which(str_detect(pngfile, i)==TRUE)]

		g <- g+(length(all.group[[i]]))
		plot(c(0,1),c(0,1),ann=FALSE, xlab="", ylab="", type="n", axes=FALSE)
		img <- readPNG(gatingfile[g-2])
		rasterImage(img,0,0,1,1, interpolate = FALSE)
		
		for(j in popFCS3){
			list.png <- temp[which(str_detect(temp,paste0(str_replace(j,'\\+','\\\\+'),"_",i))==TRUE)]
			for(k in list.png){
				plot(c(0,1),c(0,1),ann=FALSE, xlab="", ylab="", type="n", axes=FALSE)
				img <- readPNG(k)
				rasterImage(img,0,0,1,1, interpolate = FALSE)
			}
		}
		dev.off()
	}

	t2 <- Sys.time()-t1
	print(t2)

	#################################### ZIP INVERS FILE ###########################################

	t1 <- Sys.time()

	pathInvFile <- paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/invFILE")
	listInvFile <- list.files(pathInvFile, pattern=".fcs", recursive=TRUE, full.names=TRUE)
	
	## RENAME ALL INVFILE
	zip(paste0("/media/data/html/OUTPUT/COMBI/",experiment.name,"_combi/",experiment.name,"_combi.zip"),listInvFile)
	t2 <- Sys.time()-t1
	print(t2)

}






foundThresholdXMLDiva <- function(xml.file, popFCS)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	popFCS <- gsub("\\+","\\\\+",popFCS)

	data <- xmlParse(xml.file)
	xml <- xmlToList(data)

	nbGates <- length(xml$worksheet_template$gates)
	threshold <- c()
	thparam <- c()
	thpipegate <- c()
	thname <- c()
	thscale <- c()
		
	for(gate in 1:nbGates)
	{
		nameGate <- unlist(xml$worksheet_template$gates[gate]$gate$.attrs[["fullname"]])

		if(length(strsplit(nameGate,popFCS)[[1]])>2)
		{
			thparam <- c(thparam, as.vector(unlist(xml$worksheet_template$gates[gate]$gate$region$.attrs["xparm"])))
			threshold <- c(threshold, as.integer(unlist(xml$worksheet_template$gates[gate]$gate$region$points[[1]]["x"])))
			thpipegate <- c(thpipegate, as.vector(nameGate))
			thname <- c(thname, as.vector(unlist(strsplit(nameGate,popFCS))[3]))
			thscale <- c(thscale, as.numeric(unlist(xml$worksheet_template$gates[gate]$gate$x_parameter_scale_value)))
		}
	}

	return(list(thparam, thname, threshold, thscale))
}







notSubPop <- function(flow.frame, cell.population)
{

	#######################################################################################################
	#
	#	Extract a sub population from a parent population
	#
	#	Args :
	#					flow.frame : flow frame object parents which extract a sub population
	#					cell.population : cell population from flow.frame which extract 
	#	
	#	Returns:
	#					fs.not : returna flow.frame object without cell.population 
	#
	########################################################################################################
	
	if(cell.population@cell.count > 1){
		res <- which(is.na(exprs(cell.population@flow.frame)[,1])==FALSE)
		#res <- which(exprs(cell.population@flow.frame)[,1])==exprs(exprs(flow.frame))
		fs.not <- flow.frame[-res]
		return(fs.not)
	} else {
		return(getflowFrameCIPHE(cell.population))
	}
}






gatingGIF <- function(expeiment.name, time=1)
{

	#######################################################################################################
	#
	#	Convert gating auto 01 in GIF file to quick visual 
	#
	#	Args :
	#					experiment.name : name of experiment to find the png file
	#					time : time beteween each picture in the gif file
	#	
	#	Returns:
	#					Create GIF file 
	#
	########################################################################################################

	system(paste0("convert -delay 200 /media/data/html/OUTPUT/GATING/",experiment.name,"_autogating/\\*_automate_01.png /media/data/html/OUTPUT/GATING/",experiment.name,"_autogating.gif"))

}






selectOutlier <- function(vector, probability)
{

	#######################################################################################################
	#
	#	After standardisation of values with gaussien function,  
	#
	#	Args :
	#					vector : vector of value
	#					probatility : probability fitting with gaussian distribution
	#	
	#	Returns:
	#					list.outlier : list with two key, upper and lower index value from vector
	#
	########################################################################################################

	list.outlier <- list()

	norm.vector <- scale(vector)

	t <- qnorm(probability)

	z <- ((-t)*sd(vector))+mean(vector)

	list.outlier[["upper"]] <- which(norm.vector>z)
	list.outlier[["lower"]] <- which(norm.vector<(-z))

	return(list.outlier)

}






sumPopGate <- function(flow.frame.one, flow.frame.two)
{

	#######################################################################################################
	#
	#	Concat two flow.frame or cellpopulation in a gating strategye  
	#
	#	Args :
	#					flow.frame.one : first.flow.frame
	#					flow.frame.two : second.flow.frame
	#	
	#	Returns:
	#					flow.frame : output flow frame S4 object 
	#
	########################################################################################################

	flow.frame <- flow.frame.one

	exprs(flow.frame) <- rbind(exprs(flow.frame.one), exprs(flow.frame.two))

	return(flow.frame)
}


delPopGate <- function(flow.frame.one, flow.frame.two,m)
{

	#######################################################################################################
	#
	#	Concat two flow.frame or cellpopulation in a gating strategye  
	#
	#	Args :
	#					flow.frame.one : first.flow.frame
	#					flow.frame.two : second.flow.frame
	#					m : markers
	#	
	#	Returns:
	#					flow.frame : output flow frame S4 object 
	#
	########################################################################################################

	flow.frame <- flow.frame.one

	exprs(flow.frame) <- exprs(flow.frame.one)[-which(exprs(flow.frame.one)[,m]%in%exprs(flow.frame.two)[,m]),]

	return(flow.frame)
}





logicleThreshold <- function(threshold, flow.frame, thparam, x_scale=NULL)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	list.index <- names(unlist(lapply(thparam, function(x) return(which(flow.frame@description==x)))))
	list.index <- gsub("N","", list.index)
	list.index <- gsub("\\$P","", list.index)
	
	# if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
	# {
	# 	r.values <- unlist(lapply(list.index, function(x) 
	# 		as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
	# 	)
	# } else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
	# {
	# 	r.values <- unlist(lapply(list.index, function(x) 
	# 		as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
	# 	)
	# } else 
	# {
	# 	r.values <- rep(90, length(list.index))
	# }
	
	# w.values <- (4.5-log10(262144/abs(r.values)))/2
	# w.values[which(w.values<0)] <- 0.5
	# w.values[which(is.infinite(w.values))] <- 0.5

	# # threshold <- unlist(lapply(c(1:length(list.index)), FUN=function(x) 
	# # 	lnTransform(r=x_scale[x])@.Data(threshold[x]))
	# # 	#log10(threshold[x])*x_scale[x]
	# # ))

	# #threshold <- c(10^(threshold/750))

	# #threshold <- threshold*64
	threal <- unlist(lapply(c(1:length(list.index)), FUN=function(x) 
	# 	logicleTransform(w=w.values[x])@.Data(r.values[x])
	# 	#invLgcl <- inverseLogicleTransform(trans = biexponentialTransform()@.Data(r.values[x]))
			return(threshold[x]/(4095/4.5))
	))
	
	l <- list()
	for (k in threal)
	{
		l <- c(l, c(k))
	}
	threal <- l

	return(threal)

}







logiclTransformCiphe <- function(flow.frame)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
	markers.transform <- colnames(flow.frame@description[["SPILL"]])

	list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
	list.index <- gsub("N","", list.index)
	list.index <- gsub("\\$P","", list.index)
		
	if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
		)	
	} else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
		)	
	} else 
	{
		r.values <- rep(90, length(list.index))
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







logiclTransformDIVA <- function(flow.frame)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
	markers.transform <- colnames(flow.frame@description[["SPILL"]])

	list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
	list.index <- gsub("N","", list.index)
	list.index <- gsub("\\$P","", list.index)
		
	if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
		)	
	} else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
	{
		r.values <- unlist(lapply(list.index, function(x) 
			as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
		)	
	} else 
	{
		r.values <- rep(90, length(list.index))
	}
	
	w.values <- (4096-log10(262144/abs(r.values)))/2
	w.values[which(w.values<0)] <- 0.5
	w.values[which(is.infinite(w.values))] <- 0.5

	for(t in 1:length(markers.transform)){
		lgcl <- logicleTransform(w=w.values[t], t=262144, m=4096)
		flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
	}

	return(flow.frame)
}





# inversLogiclTransf <- function(flow.frame, rewrite = FALSE, path.rewrite=NULL)
# {
# 	#######################################################################################################
# 	#
# 	#	Rewrite flow.frame object with initial value, invers logicle bi exponential transformation
# 	#
# 	#	Args :
# 	#					flow.frame : flow.frame object from an FCS file
# 	#					rewrite : write the flow.frame with a path
# 	#					path.rewrite : path to rewrite the fcs file (with name.fcs)
# 	#
# 	#	Returns:
# 	#					flow.frame.inv : output flow frame S4 object with inverse logicle transform value
# 	#
# 	########################################################################################################

# 	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
# 	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
# 	markers.transform <- colnames(flow.frame@description[["SPILL"]])
#   # lgcl <- logicleTransform( w = 0.5, t= 10000, m =2.5)
#   #  invLgcl <- inverseLogicleTransform(trans = lgcl)
#   #  trans <- transformList(markers.transform, invLgcl)
#   #  flow.frame.inv <- transform(flow.frame, trans)

# 	list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
# 	list.index <- gsub("N","", list.index)
# 	list.index <- gsub("\\$P","", list.index)
		
# 	if(!is.null(flow.frame@description[[paste0("$P",list.index[1],"MS")]]))
# 	{
# 		r.values <- unlist(lapply(list.index, function(x) 
# 			as.integer(flow.frame@description[[paste0("$P",x,"MS")]]))
# 		)	
# 	} else if(!is.null(flow.frame@description[[paste0("P",list.index[1],"MS")]]))
# 	{
# 		r.values <- unlist(lapply(list.index, function(x) 
# 			as.integer(flow.frame@description[[paste0("P",x,"MS")]]))
# 		)		
# 	} else 
# 	{
# 		r.values <- rep(90, length(list.index))
# 	}
	
# 	w.values <- (4.5-log10(262144/abs(r.values)))/2
# 	w.values[which(w.values<0)] <- 0.5
# 	w.values[which(is.infinite(w.values))] <- 0.5
		
# 	flow.frame.inv <- flow.frame
	
# 	for(t in 1:length(markers.transform)){
# 		invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
# 		flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
# 	}

# 	if(rewrite == TRUE && path.rewrite != NULL)
# 	{
# 		write.FCS(flow.frame.inv, path.rewrite)
# 	}

# 	return(flow.frame.inv)
# }


inversLogiclTransformCIPHE <- function(flow.frame, value = NULL, markers = NULL)
{

	if(dim(flow.frame)[1]<5){
		return(flow.frame)
	}

   if(is.null(markers)){
    if(is.null(flow.frame@description[["SPILL"]])){
      markers.transform <- colnames(flow.frame)
    } else {
      markers.transform <- colnames(flow.frame@description[["SPILL"]])
    }
  } else {
    markers.transform <- markers
  }
  
  list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
  list.index <- gsub("N","", list.index)
  list.index <- gsub("\\$P","", list.index)
  
  if(is.null(value)){
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




inversLogiclTransfValue <- function(value, flow.frame, param){

	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
	markers.transform <- colnames(flow.frame@description[["SPILL"]])
	index <- names(which(flow.frame@description==param))
	index <- gsub("N","", index)
	index <- gsub("\\$P","", index)
		
	if(!is.null(flow.frame@description[[paste0("$P",index,"MS")]]))
	{
		r.value <- as.integer(flow.frame@description[[paste0("$P",index,"MS")]])
	} else if(!is.null(flow.frame@description[[paste0("P",index,"MS")]]))
	{
		r.value <- as.integer(flow.frame@description[[paste0("P",index,"MS")]])
	} else 
	{
		r.value <- rep(90, length(list.index))
	}
	
	w.value <- (4.5-log10(262144/abs(r.value)))/2
	w.value[which(w.value<0)] <- 0.5
	w.value[which(is.infinite(w.value))] <- 0.5

	inv.value <- inverseLogicleTransform(trans = logicleTransform(w=w.value))@.Data(value)

	return(inv.value)

}






logiclTransfValue <- function(value, flow.frame, param){

	# no.transform <- c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time","Flag")
	# markers.transform <- colnames(flow.frame)[colnames(flow.frame)%in%no.transform == FALSE]
	markers.transform <- colnames(flow.frame@description[["SPILL"]])

	index <- names(which(flow.frame@description==param))
	index <- gsub("N","", index)
	index <- gsub("\\$P","", index)
		
	if(!is.null(flow.frame@description[[paste0("$P",index,"MS")]]))
	{
		r.value <- as.integer(flow.frame@description[[paste0("$P",index,"MS")]])
	} else if(!is.null(flow.frame@description[[paste0("P",index,"MS")]]))
	{
		r.value <- as.integer(flow.frame@description[[paste0("P",index,"MS")]])
	} else 
	{
		r.value <- rep(90, length(list.index))
	}
	
	w.value <- (4.5-log10(262144/abs(r.value)))/2
	w.value[which(w.value<0)] <- 0.5
	w.value[which(is.infinite(w.value))] <- 0.5
	
	inv.value <- logicleTransform(w=w.value)@.Data(value)

	return(inv.value)

}






writeLogCSV <- function (experiment){
	
}






decompensate <- function(flow.frame)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	markers.transform <- colnames(flow.frame@description[["SPILL"]])
  mat <- apply(exprs(flow.frame)[,markers.transform], 1, FUN = function(x){lapply(names(x),function(y){sum(x*flow.frame@description[["SPILL"]][,y])})})
  mat <- matrix(unlist(mat), ncol=length(markers.transform), byrow=TRUE)
  exprs(flow.frame)[,markers.transform] <- mat
  return(flow.frame)

}






meansPop <- function(flow.frame,parameter)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	inv.flow.frame <- inversLogiclTransformCIPHE(flow.frame)
	# inv.flow.frame <- flow.frame
	value <- mean(exprs(inv.flow.frame)[,parameter])
	return(value)

}






outputTable <- function(table)
{

	#######################################################################################################
	#
	# list.lines = list of vetor with pop, parents, if parameters 
	#								means pos and sd for this pop and this parameter
	#								pop, parents,
	#
	#
	#######################################################################################################

	names <- names(table)
	lines <- length(names)
	
	col.names <- as.vector(unlist(lapply(names, function(x) 
		return(c(paste0(x," Events"),paste0(x," Parent"),paste0(x," Total"))))))

	return(col.names)

}





outputPop <- function(pop.vectors,fs.raw=NULL)
{
	tableOutput <- list()
	for(i in c(1:length(pop.vectors))){
		pop.name <- names(pop.vectors)[i]
		pop.name <- unlist(str_split(pop.name,"\\."))
		pop.name <- paste(pop.name,collapse=" ")
		count <- pop.vectors[[i]]@cell.count
		freq <-  round(count/dim(pop.vectors[[i]]@flow.frame)[1],3)*100
		if(is.null(fs.raw)){
			tableOutput[[pop.name]]=c(count,freq,"###")
		} else {
			freq.total <- round(count/dim(fs.raw)[1],3)*100
			tableOutput[[pop.name]]=c(count,freq,freq.total)
		}
	}
	return(tableOutput)
}




notSubPopGate <- function(flowframe, gate)
{
	not.gate <- gate
	not.gate@flow.frame <- flowframe
	exprs(not.gate@flow.frame)[gate@index,] <- NA
	not.gate@index <- which(!is.na(exprs(not.gate@flow.frame)[,1]))
	not.gate@cell.count <- length(not.gate@index)
	return(not.gate)
}


applyGateObject <- function(flow.frame, gate, marker1=NULL, marker2=NULL)
{

	#######################################################################################################
	#
	#
	#
	#
	#
	#######################################################################################################

	markers <- colnames(gate@filter)
	not.sub.gate <- notSubFrame(flow.frame, channels=markers, filter=gate@filter)
	#not.flow.frame <- not.sub.gate@flow.frame

	flow.frame.output <- flow.frame
	exprs(flow.frame.output) <- exprs(flow.frame.output)[-not.sub.gate@index,]

	return(flow.frame.output)
}




###############################################################
### ALL FLOWDENSITY FUNCTION WITH CONTROL IF EMPTY, > 2 , ....
###############################################################



concatenateCIPHE <- function(flow.frames) {
  ff.concat <- NULL
  n <- length(flow.frames)
  for(i in 1:n){
    ff.raw <- flow.frames[[i]]
   	ff.raw <- cbind2(ff.raw, matrix(i, nrow = nrow(ff.raw), ncol=1, dimnames = list(NULL, "Flag")))
   	if(is.null(ff.concat)){
      ff.concat  <- ff.raw
    } else {
      exprs(ff.concat) <- rbind(exprs(ff.concat),exprs(ff.raw))
    }
  }
  return(ff.concat)
}




getflowFrameCIPHE <- function(gate) {
	if(gate@cell.count==1 && length(gate@cell.count)>0){
		flow.frame <- gate@flow.frame
		exprs(flow.frame) <- exprs(flow.frame)[-c(1:dim(flow.frame)[1]),]
	} else {
		flow.frame <- getflowFrame(gate)
	}
	return(flow.frame)
}




flowDensityCIPHE <- function(flow.frame, 
														channels, 
														position, 
														gate = c(NA,NA),
														use.percentile = c(FALSE, FALSE),
														percentile = c(NA,NA),
														ellip.gate=FALSE,
														scale = NA
														) {
	result <- NULL
	try(
		result <- flowDensity(flow.frame,
												 channels = channels, 
												 position = position, 
												 gate = gate, 
												 use.percentile = use.percentile, 
												 percentile = percentile, 
												 ellip.gate = ellip.gate, 
												 scale = scale
												)
		, silent=TRUE)

	if(is.null(result)){
		result <- new("CellPopulation")
		result@flow.frame <- flow.frame
		result@filter <- matrix(c(NA,NA), nrow=1)
		result@cell.count <- 0
	}
	return(result)
}



mindesity2CIPHE <- function(flow.frame, marker, min=NULL, max=NULL){
	if(dim(flow.frame)[1]<10){
		return(0)
	} else {
		return(mindensity2(flow.frame, marker, min=min, max=max)@min)
	}
}




plotDensCIPHE <- function(flow.frame,
													channels,
													xlim=c(-0.5,4.5),
													ylim=c(-0.5,4.5),
													main=NULL,
													cex.lab=3.5, 
													cex.main=3.5, 
													cex.axis=2.5
													){
	xlim <- c(min(exprs(flow.frame)[,channels[1]])*0.9,max(exprs(flow.frame)[,channels[1]])*1.1)
	ylim <- c(min(exprs(flow.frame)[,channels[2]])*0.9,max(exprs(flow.frame)[,channels[2]])*1.1)
	x <- pData(flow.frame@parameters)[grep(channels[1],colnames(flow.frame)),2]
	y <- pData(flow.frame@parameters)[grep(channels[2],colnames(flow.frame)),2]
	if(dim(flow.frame)[1]<3000) {
		plot(exprs(flow.frame)[,channels[1]],exprs(flow.frame)[,channels[2]],pch=".",col="blue",xlim=xlim, ylim=ylim, main=main, cex.lab = cex.lab, cex.main=cex.main, cex.axis=cex.axis,xlab=x,ylab=y)
	}else if(dim(flow.frame)[1]>3000){
		plotDens(obj=flow.frame, channels=channels, xlim=xlim, ylim=ylim, main=main, cex.lab = cex.lab, cex.main=cex.main, cex.axis=cex.axis,xlab=x,ylab=y)
 	}else {
		plot(NULL,NULL,xlim=xlim, ylim=ylim, main=main, cex.lab = cex.lab, cex.main=cex.main, cex.axis=cex.axis,xlab=x,ylab=y)
	}
}




qualityGate_1 <- function(input)
{
  flowType <- flowType(input, 
                  c("FSC-A",colnames(input)[14]), 
                  Methods = 'kmeans',
                  PartitionsPerMarker = c(2,2)
                  )         

	reduce <- input[sample(c(1:dim(input)[1]),10000),]
  live <- flowDensityCIPHE(reduce, 
                      channels = c(1,14), 
                      position=c(TRUE,FALSE), 
                      gate = c(flowType@Thresholds[[1]], flowType@Thresholds[[2]]),
                      #gate = c(50000,2.1),
                      ellip.gate = TRUE,
                      scale=0.99999
                      )
 	output <- applyGateObject(input, live)
  return(list(live,output))
}


qualityGate_2 <- function(input, gate=NULL)
{
	reduce <- input[sample(c(1:dim(input)[1]),10000),]
  singletfsc <- flowDensityCIPHE(reduce,
                           channels = c("FSC-A","FSC-H"),
                            use.percentile = c(TRUE,TRUE),
                            position = c(FALSE,FALSE),
                            percentile = c(0.9,0.9),
                            ellip.gate=TRUE,
                            scale=.999999999999999
                           )
  output <- applyGateObject(input, singletfsc)
  return(list(singletfsc,output))
}


qualityGate_3 <- function(input, gate=NULL) 
{
	reduce <- input[sample(c(1:dim(input)[1]),10000),]
  singletssc <- flowDensityCIPHE(reduce,
                         channels = c("SSC-W","SSC-H"),
                          position=c(FALSE,FALSE),
                          use.percentile=c(FALSE,TRUE),
                          percentile = c(0.9,0.9999)
                          )
  output <- applyGateObject(input, singletssc)
  return(list(singletssc,output))
}


qualityGate_4 <- function(input)
{
	reduce <- input[sample(c(1:dim(input)[1]),10000),]
  x.th <- 50000
  size <- flowDensityCIPHE(reduce,
                        channels = c("FSC-A","SSC-A"),
                        position=c(TRUE,TRUE),
                        gate = c(x.th, NA),
                        use.percentile = c(FALSE,TRUE),
                        percentile = c(NA,0.001)
                        )
  output <- applyGateObject(input, size)        
  return(list(size,output))
}




concatenateCIPHE <- function(flow.frames) {
  ff.concat <- NULL
  for(i in 1:n){
    progress$inc(1/n, detail=paste0("File : ",i))
    ff.raw <- flow.frames[[i]]
   	ff.raw <- cbind2(ff.raw, matrix(i, nrow = nrow(ff.raw), ncol=1, dimnames = list(NULL, "Flag")))
   	if(is.null(ff.concat)){
      ff.concat  <- ff.raw
    } else {
      exprs(ff.concat) <- rbind(exprs(ff.concat),exprs(ff.raw))
    }
  }
  return(ff.concat)
}

cytofCore.updateFlowFrameKeywords <- function(flowFrame){
  
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  params = parameters(flowFrame)
  pdata = pData(params)
  for (i in 1:ncol(flowFrame)){
    s = paste("$P",i,"S",sep="");
    n = paste("$P",i,"N",sep="");
    r = paste("$P",i,"R",sep="");
    b = paste("$P",i,"B",sep="");
    e = paste("$P",i,"E",sep="");
    fcmax1 <- paste("flowCore_$P",i,"Rmax",sep="");
    fcmin1 <- paste("flowCore_$P",i,"Rmin",sep="");
    fcmax <- paste("flowCore_P",i,"Rmax",sep="");
    fcmin <- paste("flowCore_P",i,"Rmin",sep="");
    keyval=list();
    label <- pData(flowFrame@parameters)[,"desc"][i]
    if(is.na(label)) {label <- colnames(flowFrame)[i] }
    keyval[[s]] = label
    keyval[[n]] = colnames(flowFrame)[i]         
    keyval[[r]] = ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[b]] = 32;
    keyval[[e]] = "0,0";
    keyval[[fcmax1]] <- ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[fcmin1]] <- ceiling(min(exprs(flowFrame)[,i]))
    keyval[[fcmax]] <- ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
    keyval[[fcmin]] <- ceiling(min(exprs(flowFrame)[,i]))
    keyword(flowFrame) = keyval;
    
    pdata[i,"minRange"]=min(exprs(flowFrame)[,i])
    pdata[i,"maxRange"]=max(exprs(flowFrame)[,i])
    
  }
  pData(params)=pdata
  parameters(flowFrame)=params
  row.names(flowFrame@parameters) <- paste0("$P",c(1:length(row.names(flowFrame@parameters))))
  # keyval[["$DATATYPE"]] <- "F"
  return(flowFrame)
}