library(flowCore)
library(flowViz)
library(MASS)
library(flowDensity)
library(openCyto)
library(data.table)
library(mixtools)

# Cleaning gate:live,Size,singletSSC,sFSC
cleaning.CIPHE <- function(flow.frame, cleaning.step = c(1,1,1,1), cleaning.param.live = "Live/Dead",
                           cleaning.min.live = 1.5, cleaning.max.live = 2,
                           cleaning.param.size = "FSC-A",
                           cleaning.min.size = 60000, cleaning.max.size = 100000,
                           cleaning.x.param.singletFSC = "FSC-A", cleaning.y.param.singletFSC = "FSC-H",
                           cleaning.x.param.singletSSC = "SSC-W", cleaning.y.param.singletSSC = "SSC-H",
                           cleaning.singletSSC.bot.left.corner = c(0,0), cleaning.singletSSC.top.right.corner = c(85000,250000)) {
    # Verify if its a flowframe
    if(class(flow.frame)!= "flowFrame" ) {
        print("Votre fichier n'est pas un flowFrame")
        return(NULL)
    } else {

        # list of results to return
        list.cleaning <- list()
        result <- flow.frame
        if(dim(result)[1]<1000) {
            print("a")
            list.cleaning$live <- list(dim = dim(result), flowFrame = result, gate = 2)
        } else{
            cleaning.param.live <- get.marker.CIPHE(cleaning.param.live, flow.frame)
            cleaning.param.size <- get.marker.CIPHE(cleaning.param.size, flow.frame)
            cleaning.x.param.singletFSC <- get.marker.CIPHE(cleaning.x.param.singletFSC, flow.frame)
            cleaning.y.param.singletFSC <- get.marker.CIPHE(cleaning.y.param.singletFSC, flow.frame)
            cleaning.x.param.singletSSC <- get.marker.CIPHE(cleaning.x.param.singletSSC, flow.frame)
            cleaning.y.param.singletSSC <- get.marker.CIPHE(cleaning.y.param.singletSSC, flow.frame)
            if(cleaning.step[1] == 1) {
                # Gate Live
                gate.live <- openCyto:::.mindensity(result, channels = cleaning.param.live,
                                                    min = cleaning.min.live, max = cleaning.max.live)
                draw.live <- gate.live@min
                row.selected.live <- which(result@exprs[,cleaning.param.live] <= gate.live@min)
                data.gated.lived <- result@exprs[row.selected.live,]
                flow.frame.live <- result
                flow.frame.live@exprs <- data.gated.lived
                list.cleaning$live <- list(dim = dim(result), flowFrame = flow.frame.live, gate = draw.live)
                result <- flow.frame.live
            }
        }

         if(dim(result)[1]<1000) {
            print("b")
            list.cleaning$size <- list(dim = dim(result), flowFrame = result, gate = 50000)
        } else {
            if(cleaning.step[2] == 1) {
                # Gate Size
                gate.size <- openCyto:::.mindensity(result, channels = cleaning.param.size,
                                                    min = cleaning.min.size, max = cleaning.max.size)
                draw.size <- gate.size@min
                row.selected.size <- which(result[,cleaning.param.size] >= gate.size@min)
                data.gated.size <- result@exprs[row.selected.size,]
                flow.frame.size <- result
                flow.frame.size@exprs <- data.gated.size
                list.cleaning$size <- list(dim = dim(result), flowFrame = flow.frame.size, gate = draw.size)
                result <- flow.frame.size
            }
        }

        if(dim(result)[1]<1000) {
            print("c")
            list.cleaning$singlets.FSC <- list(dim = dim(result), flowFrame = result, gate = 50000)
        } else {
            if(cleaning.step[3] == 1) {
                # Gate SingletFSC
                gate.singlet.FSC <- openCyto:::.singletGate(result, channels = c(cleaning.x.param.singletFSC,cleaning.y.param.singletFSC))
                draw.singlets.FSC <- matrix(c(gate.singlet.FSC@boundaries[1,1],gate.singlet.FSC@boundaries[2,1],
                                              gate.singlet.FSC@boundaries[3,1], gate.singlet.FSC@boundaries[4,1],
                                              gate.singlet.FSC@boundaries[1,1], gate.singlet.FSC@boundaries[1,2],
                                              gate.singlet.FSC@boundaries[2,2], gate.singlet.FSC@boundaries[3,2],
                                              gate.singlet.FSC@boundaries[4,2], gate.singlet.FSC@boundaries[1,2]),
                                            ncol = 2)
                data.not.selected.singletFSC <- notSubFrame(obj = result,
                                                            channels = colnames(gate.singlet.FSC@boundaries),
                                                            filter = gate.singlet.FSC@boundaries)
                flow.frames.singlet.FSC <- result
                flow.frames.singlet.FSC@exprs <- result@exprs[-data.not.selected.singletFSC@index,]
                list.cleaning$singlets.FSC <- list(dim = dim(result), flowFrame = flow.frames.singlet.FSC,
                                                   gate = draw.singlets.FSC)
                result <- flow.frames.singlet.FSC
            }
        }
        
         if(dim(result)[1]<1000) {
            print("d")
            list.cleaning$singlets.SSC <- list(dim = dim(result), flowFrame = result, gate = 50000)
        } else {
            if(cleaning.step[4] == 1) {
                # Gate SingletSSC
                gate.SSC <- openCyto:::.boundary(result, channels = c(cleaning.x.param.singletSSC,cleaning.y.param.singletSSC),
                                                 min = cleaning.singletSSC.bot.left.corner, max = cleaning.singletSSC.top.right.corner)
                draw.singletSSC <- rbind(gate.SSC@min, gate.SSC@max)
                rect.gate <- matrix(c(gate.SSC@min[1],gate.SSC@min[2],gate.SSC@min[1],
                                      gate.SSC@max[2],gate.SSC@max[1],gate.SSC@max[2],
                                      gate.SSC@max[1],gate.SSC@min[2], gate.SSC@min[1], gate.SSC@min[2]),
                                    nrow = 5, ncol=2, byrow = TRUE)
                data.not.selected.SSC <- notSubFrame(obj = result, channels = c(cleaning.x.param.singletSSC,cleaning.y.param.singletSSC), filter = rect.gate)
                flow.frames.SSC <- result
                flow.frames.SSC@exprs <- result@exprs[-data.not.selected.SSC@index,]
                list.cleaning$singlets.SSC <- list(dim = dim(result), flowFrame = flow.frames.SSC,
                                                   gate = rect.gate) 
            }
        }
        return(list.cleaning)
    }
}

plot.cleaning.CIPHE <- function(flow.frame,list.of.flowframes, param.live, param.size = "FSC-A",
                                param.singletFSC = c("FSC-A","FSC-H"),param.singletSSC = c("SSC-W","SSC-H"), parmar = TRUE) {
    if(parmar == TRUE) {par(mfrow = c(1,4))}
    channels.live <- get.marker.CIPHE(param.live, flow.frame)
    channels.size <- get.marker.CIPHE(param.size, flow.frame)
    channels.singletFSCX <- get.marker.CIPHE(param.singletFSC[1], flow.frame)
    channels.singletFSCY <- get.marker.CIPHE(param.singletFSC[2], flow.frame)
    channels.singletSSCX <- get.marker.CIPHE(param.singletSSC[1], flow.frame)
    channels.singletSSCY <- get.marker.CIPHE(param.singletSSC[2], flow.frame)
    if(is.null(list.of.flowframes$live)) {
        list.of.flowframes$live$flowFrame <- flow.frame
        list.of.flowframes$live$gate <- NULL
    }
    if(is.null(list.of.flowframes$size)) {
        list.of.flowframes$size <- list.of.flowframes$live
        list.of.flowframes$size$gate <- NULL
    }
    if(is.null(list.of.flowframes$singlets.SSC)) {
        list.of.flowframes$singlets.SSC <- list.of.flowframes$size
        list.of.flowframes$singlets.SSC$gate <- NULL
    }
    if(is.null(list.of.flowframes$singlets.FSC)) {
        list.of.flowframes$singlets.FSC <- list.of.flowframes$singlets.SSC
        list.of.flowframes$singlets.FSC$gate <- NULL
    }
    name.plot <- "Root"
    par(mar = c(5,5,2,2))
    channels <- c("SSC-A", channels.live)
    xlim <- c(-0.5,4.5)
    ylim <- c(-0.5,4.5)
    if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)}
    if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)}
    flow.frame.tmp <- flow.frame
    exprs <- flow.frame.tmp@exprs[sample(1:dim(flow.frame.tmp)[1],(0.5*dim(flow.frame.tmp)[1])),]
    flow.frame.sampled <- flow.frame.tmp
    flow.frame.sampled@exprs <- exprs
    plotDens(flow.frame.sampled, channels = channels, main = name.plot, xlim = xlim, ylim = ylim,
             xlab = "SSC-A", ylab = param.live)
    # plotDens.CIPHE(flow.frames,channels = c("SSC-A",param.live), main = name.plot, xlim = c(-0.5,262143), ylim=c(-0.5,4.5))
    if(!is.null(list.of.flowframes$live$gate)) {
        abline(h=list.of.flowframes$live$gate, lwd = 2)
        name.plot <- "Live"
    }
    
    par(mar = c(5,5,2,2))
    channels <- c(channels.size,"SSC-A")
    xlim <- c(-0.5,4.5)
    ylim <- c(-0.5,4.5)
    if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)}
    if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)}
    flow.frame.tmp <- list.of.flowframes$live$flowFrame
    exprs <- flow.frame.tmp@exprs[sample(1:dim(flow.frame.tmp)[1],(0.5*dim(flow.frame.tmp)[1])),]
    flow.frame.sampled <- flow.frame.tmp
    flow.frame.sampled@exprs <- exprs
    plotDens(flow.frame.sampled, channels = channels, main = name.plot, xlim = xlim, ylim = ylim,
             xlab = param.size, ylab = "SSC-A")
    # plotDens.CIPHE(list.of.flowframes$live$flowFrame, channels = c(param.size,"SSC-A"), main = name.plot)
    if(!is.null(list.of.flowframes$size$gate)) {
        abline(v=list.of.flowframes$size$gate, lwd = 2)
        name.plot <- "Size"
    }
    
    par(mar = c(5,5,2,2))
    channels <- c(channels.singletFSCX,channels.singletFSCY)
    xlim <- c(-0.5,4.5)
    ylim <- c(-0.5,4.5)
    if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)}
    if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)}
    flow.frame.tmp <- list.of.flowframes$size$flowFrame
    exprs <- flow.frame.tmp@exprs[sample(1:dim(flow.frame.tmp)[1],(0.5*dim(flow.frame.tmp)[1])),]
    flow.frame.sampled <- flow.frame.tmp
    flow.frame.sampled@exprs <- exprs
    plotDens(flow.frame.sampled, channels = channels, main = name.plot,
             xlim = xlim, ylim = ylim, xlab = param.singletFSC[1], ylab = param.singletFSC[2])
    # plotDens.CIPHE(list.of.flowframes$size$flowFrame,channels = param.singletFSC, main = name.plot)
    if(!is.null(list.of.flowframes$singlets.FSC$gate)) {
        lines(list.of.flowframes$singlets.FSC$gate)
        name.plot <- "Singlets.FSC"
    }
    
    par(mar = c(5,5,2,2))
    channels <- c(channels.singletSSCX,channels.singletSSCY)
    xlim <- c(-0.5,4.5)
    ylim <- c(-0.5,4.5)
    if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)}
    if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)}
    flow.frame.tmp <- list.of.flowframes$singlets.FSC$flowFrame
    exprs <- flow.frame.tmp@exprs[sample(1:dim(flow.frame.tmp)[1],(0.5*dim(flow.frame.tmp)[1])),]
    flow.frame.sampled <- flow.frame.tmp
    flow.frame.sampled@exprs <- exprs
    plotDens(flow.frame.sampled, channels = channels, main = name.plot,
             cex.main = 1.5, xlim = xlim, ylim = ylim, xlab = param.singletSSC[1], ylab = param.singletSSC[2])
    # plotDens.CIPHE(list.of.flowframes$singlets.FSC$flowFrame, channels = param.singletSSC, main = name.plot)
    if(!is.null(list.of.flowframes$singlets.SSC$gate)) {
        lines(list.of.flowframes$singlets.SSC$gate)
        name.plot <- "Singlets.SSC"
    }
}

