library(flowCore)
library(flowViz)
library(MASS)
library(flowDensity)
library(openCyto)
library(data.table)
library(mixtools)


# function of gating ####

hist.gate.CIPHE <- function(flow.frame, channel, x.gate.range=c(0.5,4.5), threshold=NULL){
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channel, flow.frame)[1]
    if(dim(flow.frame)[1]>10){
        th <- openCyto::mindensity2(flow.frame, channel = x.param, min=x.gate.range[[1]], max=x.gate.range[[2]])@min
    } else {
        th <- 2
    }
    final.list$threshold <- list(data.1 = th, data.2 = th)
    names(final.list$threshold[[1]]) <- channel
    names(final.list$threshold[[2]]) <- channel
    p <- flow.frame[which(flow.frame[,x.param]>th),] 
    n <- flow.frame[which(flow.frame[,x.param]<th),] 
    final.list$output.fcs <- list(data.1 = p, data.2 = n)
    final.list$channels <- channel
    names(final.list$threshold) <- channel
    return(final.list)
}

apply.list.hist.CIPHE <- function(flow.frame, list.threshold) {
    progress <- Progress$new()
    progress$set(message="apply hists...")
    final.list <- list()
    final.list$input.fcs <- flow.frame
    progress$set(message=dim(flow.frame))
    params <- get.marker.CIPHE(names(list.threshold)[1], flow.frame)
    final.list$channels <- names(list.threshold)[1]
    progress$set(message=names(list.threshold))
    final.list$threshold <- list.threshold
    progress$set(message=params)
    p <- flow.frame[which(flow.frame[,params]>final.list$threshold),] 
    n <- flow.frame[which(flow.frame[,params]<final.list$threshold),]
    progress$set(message=dim(p))
    progress$set(message=dim(n))
    final.list$output.fcs <- list(p, n)
    names(final.list$threshold) <- names(list.threshold)[1]
    progress$close()
    return(final.list)
}

vertical.quantile.CIPHE <- function(flow.frame, channels, probs, th.limit, min = NULL, max = NULL){
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    gate <- openCyto:::.quantileGate(flow.frame, channels = x.param, probs = probs)@min
    
    if(th.limit[1] == th.limit[2]) {
        final.list$output.fcs <- list(data.1 = flow.frame, data.2 = flow.frame)
        matrix <- matrix(c(th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[2],th.limit[2],th.limit[2],th.limit[2],th.limit[2]),
                         ncol = 2)
        final.list$threshold <- list(data.1 = matrix, data.2 = matrix)
        colnames(final.list$threshold$data.1) <- channels
        final.list$channels <- channels
    } else {
        x <- th.limit[1]
        y <- th.limit[2]
        matrix1 <- matrix(c(y,gate,gate,y,y,y,y,x,x,y), ncol = 2)
        matrix2 <- matrix(c(gate,x,x,gate,gate,y,y,x,x,y), ncol = 2)
        colnames(matrix1) <- channels
        colnames(matrix2) <- channels
        polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T), filter = matrix1)
        polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T), filter = matrix2)
        final.list$threshold <- list(data.1 = matrix1, data.2 = matrix2)
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2))
        final.list$channels <- colnames(final.list$threshold[[1]])
    }
    return(final.list)
}

horizontal.quantile.CIPHE <- function(flow.frame, channels, probs, th.limit, min = NULL, max = NULL){
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    gate <- openCyto:::.quantileGate(flow.frame, channels = y.param, probs = probs)@min
    if(th.limit[1] == th.limit[2]) {
        final.list$output.fcs <- list(data.1 = flow.frame, data.2 = flow.frame)
        matrix <- matrix(c(th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[2],th.limit[2],th.limit[2],th.limit[2],th.limit[2]),
                         ncol = 2)
        final.list$threshold <- list(data.1 = matrix, data.2 = matrix)
        colnames(final.list$threshold$data.1) <- channels
        final.list$channels <- channels
    } else {
        x <- th.limit[1]
        y <- th.limit[2]
        matrix1 <- matrix(c(y,x,x,y,y,y,y,gate,gate,y), ncol = 2)
        matrix2 <- matrix(c(y,x,x,y,y,gate,gate,x,x,gate), ncol = 2)
        colnames(matrix1) <- channels
        colnames(matrix2) <- channels
        polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T), filter = matrix1)
        polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T), filter = matrix2)
        final.list$threshold <- list(data.1 = matrix1, data.2 = matrix2)
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2))
        final.list$channels <- colnames(final.list$threshold[[1]])
    }
    return(final.list)
}

mindensity.2D.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL, th.limit = NULL) {
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    vertical.gate <- openCyto:::.mindensity(flow.frame, channels = x.param, gate_range = x.gate.range, min = min, max = max)@min
    horizontal.gate <- openCyto:::.mindensity(flow.frame, channels = y.param, gate_range = y.gate.range, min = min, max = max)@min
    if(th.limit[1] == th.limit[2]) {
        final.list$output.fcs <- list(data.1 = flow.frame, data.2 = flow.frame, data.3 = flow.frame, data.4 = flow.frame)
        matrix <- matrix(c(th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[2],th.limit[2],th.limit[2],th.limit[2],th.limit[2]),
                          ncol = 2)
        final.list$threshold <- list(data.1 = matrix, data.2 = matrix, data.3 = matrix, data.4 = matrix)
        colnames(final.list$threshold$data.1) <- channels
        final.list$channels <- channels
    } else {
        x <- th.limit[1]
        y <- th.limit[2]
        matrix1 <- matrix(c(y,vertical.gate,vertical.gate,y,y,y,y,horizontal.gate,horizontal.gate,y), ncol = 2)
        matrix2 <- matrix(c(vertical.gate,x,x,vertical.gate,vertical.gate,y,y,horizontal.gate,horizontal.gate,y), ncol = 2)
        matrix3 <- matrix(c(vertical.gate,x,x,vertical.gate,vertical.gate,horizontal.gate,horizontal.gate,x,x,horizontal.gate), ncol = 2)
        matrix4 <- matrix(c(y,vertical.gate,vertical.gate,y,y,horizontal.gate,horizontal.gate,x,x,horizontal.gate), ncol = 2)
        colnames(matrix1) <- channels
        colnames(matrix2) <- channels
        colnames(matrix3) <- channels
        colnames(matrix4) <- channels
        polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T), filter = matrix1)
        polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T), filter = matrix2)
        polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F), filter = matrix3)
        polygon4 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F), filter = matrix4)
        final.list$threshold <- list(data.1 = matrix1, data.2 = matrix2, data.3 = matrix3, data.4 = matrix4)
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                      data.3 = get.flowFrame.CIPHE(polygon3), data.4 = get.flowFrame.CIPHE(polygon4))
        final.list$channels <- colnames(final.list$threshold[[1]])
    }
    return(final.list)
}

mindensity.vertical.CIPHE <- function(flow.frame, channels, x.gate.range, min = NULL, max = NULL, th.limit = NULL) {
    final.list <- list()
    final.list$input.fcs <- flow.frame
    final.list$output.fcs <- list()
    final.list$threshold <- list()
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    if(th.limit[1] == th.limit[2]) {
        final.list$output.fcs <- list(data.1 = flow.frame, data.2 = flow.frame)
        matrix <- matrix(c(th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[2],th.limit[2],th.limit[2],th.limit[2],th.limit[2]),
                          ncol = 2)
        final.list$threshold <- list(data.1 = matrix, data.2 = matrix)
        colnames(final.list$threshold$data.1) <- channels
        final.list$channels <- channels
    } else {
        vertical.gate <- openCyto:::.mindensity(flow.frame, channels = x.param, gate_range = x.gate.range, min = min, max = max)@min
        if(!is.null(th.limit) || length(th.limit) == 2) {
            x <- th.limit[1]
            y <- th.limit[2]
            if(x > vertical.gate) vertical.gate <- x
            if(y < vertical.gate) vertical.gate <- y
            matrix1 <- matrix(c(y,vertical.gate,vertical.gate,y,y,y,y,x,x,y), ncol = 2)
            matrix2 <- matrix(c(vertical.gate,x,x,vertical.gate,vertical.gate,y,y,x,x,y), ncol = 2)
            colnames(matrix1) <- channels
            colnames(matrix2) <- channels
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA), filter = matrix1)
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,NA), filter = matrix2)
            final.list$threshold <- list(data.1 = matrix1, data.2 = matrix2)
        } else {
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA), gates = c(vertical.gate,NA))
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,NA), gates = c(vertical.gate,NA))
            final.list$threshold <- list(data.1 = polygon1@filter, data.2 = polygon2@filter)
        }
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2))
        final.list$channels <- colnames(final.list$threshold[[1]])
    }
    return(final.list)
}

mindensity.horizontal.CIPHE <- function(flow.frame, channels, y.gate.range, min = NULL, max = NULL, th.limit = NULL) {
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    if(th.limit[1] == th.limit[2]) {
        final.list$output.fcs <- list(data.1 = flow.frame, data.2 = flow.frame)
        matrix <- matrix(c(th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[1],th.limit[2],th.limit[2],th.limit[2],th.limit[2],th.limit[2]),
                          ncol = 2)
        final.list$threshold <- list(data.1 = matrix, data.2 = matrix)
        colnames(final.list$threshold$data.1) <- channels
        final.list$channels <- channels
    } else {
        horizontal.gate <- openCyto:::.mindensity(flow.frame, channels = y.param, gate_range = y.gate.range, min = min , max = max)@min
        if(!is.null(th.limit) || length(th.limit) == 2) {
            x <- th.limit[1]
            y <- th.limit[2]
            if(x > horizontal.gate) horizontal.gate <- x
            if(y < horizontal.gate) horizontal.gate <- y
            matrix1 <- matrix(c(y,x,x,y,y,y,y,horizontal.gate,horizontal.gate,y),ncol = 2)
            matrix2 <- matrix(c(y,x,x,y,y,horizontal.gate,horizontal.gate,x,x,horizontal.gate), ncol = 2)
            colnames(matrix1) <- channels
            colnames(matrix2) <- channels
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,T), filter = matrix1)
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,F), filter = matrix2)
            final.list$threshold <- list(data.1 = matrix1, data.2 = matrix2)
        } else {
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,T), gates = c(NA,horizontal.gate))
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,F), gates = c(NA,horizontal.gate))
            final.list$threshold <- list(data.1 = polygon1@filter, data.2 = polygon2@filter)
        }
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2))
        final.list$channels <- colnames(final.list$threshold[[1]])
    }
    return(final.list)
}

mindensity.3data.xpos.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL,
                                        th.limit = NULL) {
    final.list <- list()
    list.4.data <- list()
    list.4.data <- mindensity.2D.CIPHE(flow.frame = flow.frame, channels = channels, x.gate.range = x.gate.range,
                                       y.gate.range = y.gate.range, min = min , max = max, th.limit = th.limit)
    list.2.data <- mindensity.vertical.CIPHE(flow.frame, channels = channels, x.gate.range = x.gate.range, th.limit = th.limit)
    final.list$input.fcs <- flow.frame
    final.list$threshold <- list(data.1 = list.2.data$threshold$data.1, data.2 = list.4.data$threshold$data.2,
                                 data.3 = list.4.data$threshold$data.3)
    final.list$output.fcs <- list(data.1 = list.2.data$output.fcs$data.1, data.2 = list.4.data$output.fcs$data.2,
                                  data.3 = list.4.data$output.fcs$data.3)
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

mindensity.3data.xneg.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL,
                                        th.limit = th.limit) {
    final.list <- list()
    list.4.data <- list()
    list.4.data <- mindensity.2D.CIPHE(flow.frame = flow.frame, channels = channels, x.gate.range = x.gate.range,
                                       y.gate.range = y.gate.range, min = min , max = max, th.limit = th.limit)
    list.2.data <- mindensity.vertical.CIPHE(flow.frame, channels = channels, x.gate.range = x.gate.range, th.limit = th.limit)
    final.list$input.fcs <- flow.frame
    final.list$threshold <- list(data.1 = list.4.data$threshold$data.1, data.2 = list.2.data$threshold$data.2,
                                 data.3 = list.4.data$threshold$data.4)
    final.list$output.fcs <- list(data.1 = list.4.data$output.fcs$data.1, data.2 = list.2.data$output.fcs$data.2,
                                  data.3 = list.4.data$output.fcs$data.4)
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

mindensity.3data.ypos.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL,
                                        th.limit = NULL) {
    final.list <- list()
    list.4.data <- list()
    list.4.data <- mindensity.2D.CIPHE(flow.frame = flow.frame, channels = channels, x.gate.range = x.gate.range,
                                       y.gate.range = y.gate.range, min = min , max = max, th.limit = th.limit)
    list.2.data <- mindensity.horizontal.CIPHE(flow.frame, channels = channels,y.gate.range = y.gate.range, th.limit = th.limit)
    final.list$input.fcs <- flow.frame
    final.list$threshold <- list(data.1 = list.2.data$threshold$data.1, data.2 = list.4.data$threshold$data.3,
                                 data.3 = list.4.data$threshold$data.4)
    final.list$output.fcs <- list(data.1 = list.2.data$output.fcs$data.1, data.2 = list.4.data$output.fcs$data.3,
                                  data.3 = list.4.data$output.fcs$data.4)
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

mindensity.3data.yneg.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL,
                                        th.limit = NULL) {
    final.list <- list()
    list.4.data <- list()
    list.2.data <- list()
    list.4.data <- mindensity.2D.CIPHE(flow.frame = flow.frame, channels = channels, x.gate.range = x.gate.range,
                                       y.gate.range = y.gate.range, min = min , max = max, th.limit = th.limit)
    list.2.data <- mindensity.horizontal.CIPHE(flow.frame, channels = channels,y.gate.range = y.gate.range, th.limit = th.limit)
    final.list$input.fcs <- flow.frame
    final.list$threshold <- list(data.1 = list.4.data$threshold$data.1, data.2 = list.4.data$threshold$data.2,
                                 data.3 = list.2.data$threshold$data.2)
    final.list$output.fcs <- list(data.1 = list.4.data$output.fcs$data.1, data.2 = list.4.data$output.fcs$data.2,
                                  data.3 = list.2.data$output.fcs$data.2)
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

mindensity.5data.xpos.CIPHE <- function(flow.frame, channels, x.gate.range, y.gate.range, min = NULL, max = NULL,
                                        th.limit = NULL) {
    final.list <- list()
    list.4.data <- list()
    list.2.data <- list()
    list.4.data <- mindensity.2D.CIPHE(flow.frame = flow.frame, channels = channels, x.gate.range = x.gate.range,
                                       y.gate.range = y.gate.range, min = min , max = max, th.limit = th.limit)
    list.2.data <- mindensity.vertical.CIPHE(flow.frame, channels = channels, x.gate.range = x.gate.range, th.limit = th.limit)
    final.list$input.fcs <- flow.frame
    final.list$threshold <- list(data.1 = list.4.data$threshold$data.1, data.2 = list.4.data$threshold$data.2,
                                 data.3 = list.4.data$threshold$data.3, data.4 = list.4.data$threshold$data.4,
                                 data.5 = list.2.data$threshold$data.1)
    final.list$output.fcs <- list(data.1 = list.4.data$output.fcs$data.1, data.2 = list.4.data$output.fcs$data.2,
                                  data.3 = list.4.data$output.fcs$data.3, data.4 = list.4.data$output.fcs$data.4,
                                  data.5 = list.2.data$output.fcs$data.1)
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

ellipse.gate.CIPHE <- function(flow.frame, channels, position, percentile = NULL, scale, x.gate.range = NULL, y.gate.range = NULL,
                               filter = FALSE, gate = NULL) {
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    params <- colnames(flow.frame)
    labels <- pData(flow.frame@parameters)[,2]
    labels[which(is.na(labels))] <- params[which(is.na(labels))]
    labels[which(labels=="<NA>")] <- params[which(labels=="<NA>")]
    names(params) <- labels
    if(channels[1]%in%labels){x.param <- params[[channels[1]]]} else {x.param <- colnames(flow.frame)[1]}
    if(channels[2]%in%labels){y.param <- params[[channels[2]]]} else {y.param <- colnames(flow.frame)[1]}

    gate.range <- c(-0.5,4.5)
    if(x.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
    if(!is.null(x.gate.range)) {gate.range <- x.gate.range}
    vertical.gate <- openCyto:::.mindensity(fr = flow.frame, channels = x.param, gate_range = gate.range)@min
    gate.range <- c(-0.5,4.5)
    if(y.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
    if(!is.null(y.gate.range)) {gate.range <- y.gate.range}
    horizontal.gate <- openCyto:::.mindensity(fr = flow.frame, channels = y.param, gate_range = gate.range)@min

        if(is.null(gate)) {
            if(dim(flow.frame)[1] <10) {
                final.list$output.fcs <- list(data.1 = flow.frame, not.data = flow.frame)
                final.list$threshold <- list(data.1 = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2),
                                             not.data = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2))
            } else {
                if(is.null(filter)){filter<-FALSE}
                gate <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = position,
                                          ellip.gate = TRUE, gate = c(vertical.gate, horizontal.gate), scale = scale, filter = filter)
                final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(gate)
                not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
                                              gates = gate@gates, filter = gate@filter)
                final.list$output.fcs$not.data <- get.flowFrame.CIPHE(not.data)
                final.list$threshold$data.1 <- gate@filter
                colnames(final.list$threshold$data.1) <- channels
            }
        } else {
            if(dim(flow.frame)[1] <10) {
                final.list$output.fcs <- list(data.1 = flow.frame, not.data = flow.frame)
                final.list$threshold <- list(data.1 = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2),
                                             not.data = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2))
            } else {
                if(is.null(filter)){filter<-FALSE}
                gate <- flowDensity.CIPHE(flow.frame =  flow.frame, channels = c(x.param,y.param), position = position, 
                                          ellip.gate = TRUE, scale = scale, filter = filter, gate = gate)
                final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(gate)
                not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
                                              gates = gate@gates, filter = gate@filter)
                final.list$output.fcs$not.data <- get.flowFrame.CIPHE(not.data)
                final.list$threshold$data.1 <- gate@filter
                colnames(final.list$threshold$data.1) <- channels
            }
        }

    colnames(final.list$threshold$data.1) <- channels
    final.list$channels <- channels
    return(final.list)
}

# ellipse.gate.CIPHE <- function(flow.frame, channels, position, percentile = NULL, scale, x.gate.range = NULL, y.gate.range = NULL,
#                                filter = NULL, gate = NULL) {
#     final.list <- list()
#     final.list$input.fcs <- flow.frame
#     if(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")) {x.param <- channels[1]}
#     else {
#         x.param <- get.labels.CIPHE(channels[1], flow.frame)
#         if(length(x.param) == 0) x.param <- as.vector(get.markers.CIPHE(flow.frame)[1])
#     }
#     if(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")) {y.param <- channels[2]}
#     else {
#         y.param <- get.labels.CIPHE(channels[2], flow.frame)
#         if(length(y.param) == 0) y.param <- as.vector(get.markers.CIPHE(flow.frame)[1])
#     }
#     if(is.null(gate)) {
#         gate.range <- c(-0.5,4.5)
#         if (x.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
#         if(!is.null(x.gate.range)) {gate.range <- x.gate.range}
#         vertical.gate <- openCyto:::.mindensity(fr = flow.frame, channels = x.param, gate_range = gate.range)@min
#         gate.range <- c(-0.5,4.5)
#         if (y.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
#         if(!is.null(y.gate.range)) {gate.range <- y.gate.range}
#         horizontal.gate <- openCyto:::.mindensity(fr = flow.frame, channels = y.param, gate_range = gate.range)@min
#     }
#     if(x.param == y.param) {
#         final.list$output.fcs <- list(data.1 = flow.frame, not.data = flow.frame)
#         final.list$threshold$data.1 <- matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2)
#         colnames(final.list$threshold$data.1) <- channels
#         final.list$channels <- channels
#     } else {
#         if(is.null(gate)) {
#             if(dim(flow.frame)[1] <10) {
#                 final.list$output.fcs <- list(data.1 = flow.frame, not.data = flow.frame)
#                 final.list$threshold <- list(data.1 = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2),
#                                              not.data = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2))
#             } else {
#                 gate <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = position,
#                                           ellip.gate = TRUE, gate = c(vertical.gate, horizontal.gate), scale = scale, filter = filter)
#                 final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(gate)
#                 not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
#                                               gates = gate@gates, filter = gate@filter)
#                 final.list$output.fcs$not.data <- get.flowFrame.CIPHE(not.data)
#                 final.list$threshold$data.1 <- gate@filter
#                 colnames(final.list$threshold$data.1) <- channels
#             }
#         } else {
#             if(dim(flow.frame)[1] <10) {
#                 final.list$output.fcs <- list(data.1 = flow.frame, not.data = flow.frame)
#                 final.list$threshold <- list(data.1 = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2),
#                                              not.data = matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2))
#             } else {
#                 gate <- flowDensity.CIPHE(flow.frame =  flow.frame, channels = c(x.param,y.param), position = position, 
#                                           ellip.gate = TRUE, scale = scale, filter = filter, gate = gate)
#                 final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(gate)
#                 not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
#                                               gates = gate@gates, filter = gate@filter)
#                 final.list$output.fcs$not.data <- get.flowFrame.CIPHE(not.data)
#                 final.list$threshold$data.1 <- gate@filter
#                 colnames(final.list$threshold$data.1) <- channels
#             }
#         }
#     }
#     colnames(final.list$threshold$data.1) <- channels
#     final.list$channels <- channels
#     return(final.list)
# }

ellipse.4.gate.CIPHE <- function(flow.frame, channels, scale, filter = FALSE, gate = NULL, x.gate.range = NULL, y.gate.range = NULL) {
    final.list <- list()
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    gate.range <- c(-0.5,4.5)
    if (x.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
    if(!is.null(x.gate.range)) {gate.range <- x.gate.range}
    vertical.gate <- openCyto:::.mindensity(fr = flow.frame, channels = x.param, gate_range = gate.range)@min
    gate.range <- c(-0.5,4.5)
    if (y.param%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") == TRUE) gate.range <- c(-0.5,250000)
    if(!is.null(y.gate.range)) {gate.range <- y.gate.range} 
    horizontal.gate <- openCyto:::.mindensity(fr = flow.frame, channels = y.param, gate_range = gate.range)@min
    if(is.null(gate)) {
        gate1.T.T <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[1]], gate = c(vertical.gate,horizontal.gate))
        gate2.F.T <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[2]], gate = c(vertical.gate,horizontal.gate))
        gate3.F.F <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[3]], gate = c(vertical.gate,horizontal.gate))
        gate4.T.F <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[4]], gate = c(vertical.gate,horizontal.gate))
        final.list$input.fcs <- flow.frame
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(gate1.T.T), data.2 = get.flowFrame.CIPHE(gate2.F.T),
                                      data.3 = get.flowFrame.CIPHE(gate3.F.F), data.4 = get.flowFrame.CIPHE(gate4.T.F))
        final.list$threshold <- list(data.1 = gate1.T.T@filter, data.2 = gate2.F.T@filter,
                                     data.3 = gate3.F.F@filter, data.4 = gate4.T.F@filter) 
    } else {
        gate1.T.T <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[1]],
                                       gate = gate)
        gate2.F.T <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[2]],
                                       gate = gate)
        gate3.F.F <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[3]],
                                       gate = gate)
        gate4.T.F <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                       ellip.gate = TRUE, scale = scale, filter = filter[[4]],
                                       gate = gate)
        final.list$input.fcs <- flow.frame
        final.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(gate1.T.T), data.2 = get.flowFrame.CIPHE(gate2.F.T),
                                      data.3 = get.flowFrame.CIPHE(gate3.F.F), data.4 = get.flowFrame.CIPHE(gate4.T.F))
        final.list$threshold <- list(data.1 = gate1.T.T@filter, data.2 = gate2.F.T@filter,
                                     data.3 = gate3.F.F@filter, data.4 = gate4.T.F@filter)
    }
    colnames(final.list$threshold$data.1) <- channels
    colnames(final.list$threshold$data.2) <- channels
    colnames(final.list$threshold$data.3) <- channels
    colnames(final.list$threshold$data.4) <- channels
    final.list$channels <- channels
    return(final.list)
}

poly.ellipse.gate.CIPHE <- function(flow.frame, channels, scale, filter = FALSE,
                                    doublepos = NULL, negpos = NULL, doubleneg = NULL, posneg = NULL) {
    final.list <- list()
    if(!is.null(doublepos)) {
        gate.double.positive <- flowDensity.CIPHE(flow.frame = flow.frame, channels = channels, position = c(T,T), ellip.gate = TRUE, scale = scale)
    }
    if(!is.null(negpos)) {
        gate.neg.pos <- flowDensity.CIPHE(flow.frame = flow.frame, channels = channels, position = c(F,T), ellip.gate = TRUE, scale = scale)
    }
    if(!is.null(doubleneg)) {
        gate.double.neg <- flowDensity.CIPHE(flow.frame = flow.frame, channels = channels, position = c(F,F), ellip.gate = TRUE, scale = scale)
    }
    if(!is.null(posneg)) {
        gate.pos.neg <- flowDensity.CIPHE(flow.frame = flow.frame, channels = channels , position = c(T,F), ellip.gate = TRUE, scale = scale)
    }
    final.list$input.fcs <- flow.frame
    return(final.list)
}

tailgate.CIPHE <- function(flow.frame, channel, min = NULL, max = NULL, tol = NULL) {
    gate <- openCyto:::.tailgate(flow.frame, channels = channel, min = min , max = max, tol = tol, num_peaks = 1)
    return(gate)
}

boundary.CIPHE <- function(flow.frame, channels, probs = NULL, bot.left.corner = c(0,0), top.right.corner = c(2.5e5,2.5e5), polygon = NULL) {
    final.list <- list()
    final.list$threshold <- list()
    final.list$output.fcs <- list()
    final.list$input.fcs <- flow.frame
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    if(is.null(polygon)) {
        gate <- openCyto:::.boundary(flow.frame, channels = c(x.param, y.param), min = c(as.numeric(bot.left.corner)), max = c(as.numeric(top.right.corner)),
                                     probs = probs)
        matrix <- matrix(c(gate@max[1],gate@min[1],gate@min[1],gate@max[1],gate@max[1],
                           gate@max[2],gate@max[2],gate@min[2],gate@min[2],gate@max[2]), ncol = 2)
        polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T), filter = matrix)
        final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(polygon1)
        not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
                                      gates = polygon1@gates, filter = polygon1@filter)
        final.list$output.fcs$data.2 <- get.flowFrame.CIPHE(not.data)
        final.list$threshold$data.1 <- matrix
        colnames(final.list$threshold$data.1) <- channels
    } else {
        polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T), filter = polygon)
        final.list$output.fcs$data.1 <- get.flowFrame.CIPHE(polygon1)
        not.data <- notSubFrame.CIPHE(obj = flow.frame, channels = c(x.param,y.param),
                                      gates = polygon1@gates, filter = polygon1@filter)
        final.list$output.fcs$data.2 <- get.flowFrame.CIPHE(not.data)
        final.list$threshold$data.1 <- polygon
    }
    final.list$channels <- channels
    return(final.list)
}

transitional.CIPHE <- function(flow.frame, channel, nb.clusters = 2, target = NULL, quantile = NULL, translation = NULL,
                               min = NULL, max = NULL, angle = NULL) {
    gate <- openCyto:::.flowClust.2d(flow.frame, channels = channel, K = nb.clusters, target = target,
                                     quantile = quantile, translation = translation, transitional = TRUE, min = min,
                                     max = max, transitional_angle = pi * angle)
    return(gate)
}

notSubFrame.CIPHE <- function(obj, channels, gates, filter) {

    res <- tryCatch({
        notSubFrame(obj = obj, channels = channels, gates = gates, filter = filter)
    }, error = function(e) {
        new("CellPopulation", flow.frame = obj, proportion = 100, cell.count = dim(obj)[1], channels = channels,
            filter = matrix(c(1,1,1,1,1,1,1,1,1,1), ncol = 2, dimnames = list(c(1:5),
                     c(channels[1],channels[2]))))
    })
    return(res)
}

flowDensity.CIPHE <- function(flow.frame, channels, position, ellip.gate = FALSE, scale = NULL, gate = NULL, filter = FALSE) {
    if(is.null(gate)) {
        res <- tryCatch({
            flowDensity(obj = flow.frame, channels = channels , position = position, percentile = c(0.95,0.95), use.percentile = c(T,T),
                        ellip.gate = ellip.gate, scale = scale, filter = filter)
        }, error = function(e) {
            new("CellPopulation", flow.frame = flow.frame, proportion = 100, cell.count = dim(flow.frame)[1], channels = channels,
                filter = matrix(c(1,1,1,1,1,1,1,1,1,1), ncol = 2, dimnames = list(c(1:5), 
                         c(channels[1],channels[2]))))
        })
    } else {
        res <- tryCatch({
            flowDensity(obj = flow.frame, channels = channels , position = position,
                        ellip.gate = ellip.gate, scale = scale, gate = gate, filter = filter, percentile = c(0.95,0.95), use.percentile = c(T,T))
        }, error = function(e) {
            new("CellPopulation", flow.frame = flow.frame, proportion = 100, cell.count = dim(flow.frame)[1], channels = channels,
                gates = gate, filter = matrix(c(1,1,1,1,1,1,1,1,1,1), ncol = 2,dimnames = list(c(1:5), 
                                       c(channels[1],channels[2]))))
        }) 
    }
    return(res)
    
}

apply.list.polygon.CIPHE <- function(flow.frame, list.threshold) {
    final.list <- list()
    final.list$input.fcs <- flow.frame
    final.list$channels[1] <- get.marker.CIPHE(colnames(list.threshold[[1]])[1], flow.frame)
    final.list$channels[2] <- get.marker.CIPHE(colnames(list.threshold[[1]])[2], flow.frame)
    threshold <- lapply(c(1:length(list.threshold)), function(i) {
        flowDensity.CIPHE(flow.frame = flow.frame, channels = final.list$channels, position = c(T,NA), filter = list.threshold[[i]])
    })
    final.list$output.fcs <- lapply(c(1:length(list.threshold)), function(i) {i <- get.flowFrame.CIPHE(threshold[[i]])})
    final.list$threshold <- lapply(c(1:length(list.threshold)), function(i) {threshold[[i]]@filter})
    final.list$channels <- colnames(final.list$threshold[[1]])
    return(final.list)
}

# transformation ####
logiclTransformCiphe <- function(flow.frame) {
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

data.processing.CIPHE <- function(flow.frame) {
    # source("myHelpFunction_V2.R")
  mat <- flow.frame@description[["SPILL"]]
  fcs.comp <- compensate(flow.frame, mat)
  fcs.trans <- logiclTransformCiphe(fcs.comp)
  return(fcs.trans)
}

preproces <- function(flow.frame){
    flow.frame <- compensate(flow.frame, flow.frame@description[["SPILL"]])
    list.index <- names(unlist(lapply(markers.transform, function(x) return(which(flow.frame@description==x)))))
    list.index <- gsub("N","", list.index)
    list.index <- gsub("\\$P","", list.index)
              
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
              
    w.values <- (4.5-log10(262143/abs(r.values)))/2
    w.values[which(w.values<0)] <- 0.5
    w.values[which(is.infinite(w.values))] <- 0.5
                  
    for(t in 1:length(markers.transform)){
        lgcl <- logicleTransform(w=w.values[t])
        flow.frame <- transform(flow.frame, transformList(markers.transform[t],lgcl))
    }
    return(flow.frame)
}

deCompensateFlowFrame <- function(x, spillover) {
    cols <- colnames(spillover)
    sel <- cols %in% colnames(x)
    if(!all(sel)) {
        stop(keyword(x)[["FILENAME"]], "\\nThe following parameters in the spillover matrix are not present in the flowFrame:\\n",
             paste(cols[!sel], collapse=", "), call.=FALSE)
    }
    e <- exprs(x)
    e[, cols] <- e[, cols] %*% spillover
    exprs(x) = e
    x
}

inversLogiclTransf <- function(flow.frame, rewrite = FALSE, path.rewrite=NULL) {
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
    
    w.values <- (4.5-log10(262144/abs(r.values)))/2
    w.values[which(w.values<0)] <- 0.5
    w.values[which(is.infinite(w.values))] <- 0.5
    
    flow.frame.inv <- flow.frame
    
    for(t in 1:length(markers.transform)){
        invLgcl <- inverseLogicleTransform(trans = logicleTransform(w=w.values[t]))
        flow.frame.inv <- transform(flow.frame.inv, transformList(markers.transform[t],invLgcl))
    }
    
    if(rewrite == TRUE && is.null(path.rewrite))
    {
        write.FCS(flow.frame.inv, path.rewrite)
    }
    return(flow.frame.inv)
}


get.flowFrame.CIPHE <- function(gate) {
    if(gate@cell.count <2){
        flow.frame <- gate@flow.frame
        exprs(flow.frame) <- exprs(flow.frame)[-c(1:dim(flow.frame)[1]),]
    } else {flow.frame <- getflowFrame(gate)}
    return(flow.frame)
}

get.markers.CIPHE <- function(fcs) {
    labels <- pData(fcs@parameters)[,2]
    names <- pData(fcs@parameters)[,1]
    names(names) <- labels
    markers <- names[!is.na(labels)]
    markers <- markers[which(markers!="NA")]
    return(markers)
}

get.markers.2.CIPHE <- function(fcs) {
    labels <- pData(fcs@parameters)[,2]
    labels <- unlist(lapply(c(1:length(labels)), function(i) {
        if(is.na(labels[i])) {return(pData(fcs@parameters)[,1][i])}
        else {return(labels[i])}
    }))
    names(labels) <- pData(fcs@parameters)[,1]
    return(labels)
}

get.names.CIPHE <- function(fcs) {
    labels <- pData(fcs@parameters)[,2]
    names <- labels[which(labels!="NA")]
    return(names)
}

get.labels.CIPHE <- function(marker, fcs) {
    markers <- get.markers.CIPHE(fcs)
    label <- markers[which(names(markers) == marker)]
    return(label)
}

get.labels.2.CIPHE <- function(marker, fcs) {
    markers <- get.markers.CIPHE(fcs)
    label <- names(markers[which(markers == marker)])
    return(label)
} 

get.marker.CIPHE <- function(label, fcs) {
    markers <- get.markers.2.CIPHE(fcs)
    marker <- names(markers[which(markers == label)])
    return(marker)
}    

check.flow.frames.CIPHE <- function(list) {
    if(dim(list$input.fcs)[1] <10) return(TRUE)
    else return(FALSE)
}

check.flow.frames.output.CIPHE <- function(list) {
    final.list <- list
    final.list$output.fcs <- lapply(list$output.fcs, function(i) {
        if(dim(i)[1] <10) i <- list$input.fcs
        else i <- i
    })
    return(final.list)
}

check.plot.CIPHE <- function(list) {
    if(dim(list$input.fcs)[1] <10) return(TRUE)
    else return(FALSE)
    if(list$channels != 2) return(TRUE)
    else return(FALSE)
}

draw.without.threshold.CIPHE <- function(list, nb, channels) {
    final.list <- lapply(c(1:nb), function(i){
        matrix(c(4.5,-0.5,-0.5,4.5,4.5,4.5,4.5,-0.5,-0.5,4.5), ncol = 2, dimnames = list(letters[1:5],channels))
    })
    return(final.list)
}

shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {

    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}

#Plot ####
## ON OFF POUR LES NOM DES GATES (variuable binaire)
plot.framework.CIPHE <- function(list, xlim = c(-0.5,4.5), ylim = c(-0.5,4.5), 
    percent.sample = 0.1, plot.name.gate = FALSE) {
    par(mar = c(5,5,2,2))
    if(check.plot.CIPHE(list) == TRUE) {
        plot(1, col = "white",xlim = xlim, ylim = ylim, xlab = list$channels[1], ylab = list$channels[2],main = list$name)
        text(2, 2, "Less than 10 evts")
    } else {
        if(length(list$channels)<2) {
            markers <- get.markers.2.CIPHE(list$input.fcs)
            x.param <- names(markers[which(markers == list$channels[1])])
            hist(list$input.fcs@exprs[,x.param], breaks = 100, border = "cyan", col = "cyan", xlab = list$channels, main = list$name,
                 xlim = c(-0.5,4.5))
            abline(v=list$draw, lty=2, lwd=2)
        } else {
            if(length(which(list$channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)} 
            if(length(which(list$channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)} 
            plotDens.CIPHE(flow.frame = list$input.fcs, channels = list$channels, main = list$name,
                           xlim = xlim, ylim = ylim, percent.sample = percent.sample)
        
            if(!is.null(list$draw)){
                lines(list$draw, lwd = 2)
                if(plot.name.gate==TRUE){
                    shadowtext(mean(list$draw[,1]),mean(list$draw[,2]),
                        names(list$output.fcs)[1],cex=0.75
                    )
                    shadowtext(mean(list$draw[,1]),mean(list$draw[,2]-0.4),
                        round(dim(list$output.fcs[[1]])[1]/dim(list$input.fcs)[1]*100,1),cex=0.75
                    )
                }
            }
            if(!is.null(list$draw2)){
                lines(list$draw2, lwd = 2)
                if(plot.name.gate==TRUE){
                    shadowtext(mean(list$draw2[,1]),mean(list$draw2[,2]),
                        names(list$output.fcs)[2],cex=0.75
                    )
                    shadowtext(mean(list$draw2[,1]),mean(list$draw2[,2]-0.4),
                        round(dim(list$output.fcs[[2]])[1]/dim(list$input.fcs)[1]*100,1),cex=0.75
                    )
                }
            }
            if(!is.null(list$draw3)){
                lines(list$draw3, lwd = 2)
                if(plot.name.gate==TRUE){
                    shadowtext(mean(list$draw3[,1]),mean(list$draw3[,2]),
                        names(list$output.fcs)[3],cex=0.75
                    )
                    shadowtext(mean(list$draw3[,1]),mean(list$draw3[,2]-0.4),
                        round(dim(list$output.fcs[[3]])[1]/dim(list$input.fcs)[1]*100,1),cex=0.75
                    )
                }
            }
            if(!is.null(list$draw4)){
                lines(list$draw4, lwd = 2)
                if(plot.name.gate==TRUE){
                    shadowtext(mean(list$draw4[,1]),mean(list$draw4[,2]),
                        names(list$output.fcs)[4],cex=0.75
                    )
                    shadowtext(mean(list$draw4[,1]),mean(list$draw4[,2]-0.4),
                        round(dim(list$output.fcs[[4]])[1]/dim(list$input.fcs)[1]*100,1),cex=0.75
                    )
                }
            }
        }
    }
}

plot.framework.CIPHE.LowQ <- function(list, xlim = c(-0.5,4.5), ylim = c(-0.5,4.5)){
    par(mar = c(5,5,2,2))
    if(check.plot.CIPHE(list) == TRUE) {
        plot(1, col = "white",xlim = xlim, ylim = ylim, xlab = list$channels[1], ylab = list$channels[2],main = list$name)
        text(2, 2, "Less than 10 evts")
    } else {
        if(length(list$channels) <2) {
            hist(list$input.fcs@exprs[,list$channels], breaks = 100, border = "cyan", col = "cyan", xlab = list$channels, main = list$name,
                 xlim = c(-0.5,4.5))
        } else {
            if(length(which(list$channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)} 
            if(length(which(list$channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)} 
            plot(list$input.fcs@exprs[,c(7,8)],col="blue", main = list$name, xlim = xlim, ylim = ylim,
                xlab = names(list$channels)[1], ylab = names(list$channels)[2], pch=".", cex=0.5, cex.lab = 1, cex.main = 1, cex.axis = 1)
        }
        if(!is.null(list$draw)) lines(list$draw, lwd = 2)
        if(!is.null(list$draw2)) lines(list$draw2, lwd = 2)
        if(!is.null(list$draw3)) lines(list$draw3, lwd = 2)
        if(!is.null(list$draw4)) lines(list$draw4, lwd = 2)
    }
}

plotDens.CIPHE <- function(flow.frame, channels, xlim = c(-0.5,4.5), ylim = c(-0.5,4.5),
                           main = NULL, cex.lab = 1, cex.main = 1, cex.axis = 1, percent.sample = percent.sample) {
    if(is.null(dim(flow.frame))) return(NULL)
    markers <- get.markers.2.CIPHE(flow.frame)
    x.param <- names(markers[which(markers == channels[1])])
    y.param <- names(markers[which(markers == channels[2])])

    if(length(x.param) == 0) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    if(length(y.param) == 0) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    if(dim(flow.frame)[1] <=10000) {
        plot(exprs(flow.frame)[,x.param], exprs(flow.frame)[,y.param], pch = ".", col = "blue", main = main,
             cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, xlab = channels[1], ylab = channels[2], xlim = xlim, ylim = ylim)
    }else if(dim(flow.frame)[1] >10000){
        exprs <- flow.frame@exprs[sample(1:dim(flow.frame)[1],(percent.sample*dim(flow.frame)[1])),]
        flow.frame.sampled <- flow.frame
        flow.frame.sampled@exprs <- exprs
        if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)} 
        if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)} 
        plotDens(obj = flow.frame.sampled, channels = c(x.param, y.param), main = main, cex.lab = cex.lab, cex.main = cex.main, 
                 cex.axis = cex.axis,xlab = channels[1], ylab = channels[2], xlim = xlim, ylim = ylim, pch = ".")
    }else if(dim(flow.frame)[1]<500){
        if(length(which(channels[1]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){xlim = c(0,250000)} 
        if(length(which(channels[2]%in%c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")))){ylim = c(0,250000)} 
        plotDens(obj = flow.frame, channels = c(x.param, y.param), main = main, cex.lab = cex.lab, cex.main = cex.main, 
                 cex.axis = cex.axis,xlab = channels[1], ylab = channels[2], xlim = xlim, ylim = ylim, pch = ".")
    } else {
        plot(NULL, NULL, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, xlab = channels[1], ylab = channels[2],
             xlim = xlim, ylim = ylim, main = "Less than 50 evts") 
    }
}

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

addDensityFlagCIPHE <- function(flow.frame) {
    column <- c()
    for(i in unique(exprs(flow.frame)[,"Flag"])){
        column <- c(column,seq(i-0.3,i+0.3,length.out=length(which(exprs(flow.frame)[,"Flag"]==i))))
    }
    flow.frame <- cbind2(flow.frame, matrix(column, nrow = length(column), ncol=1, dimnames = list(NULL, "FlagDens")))
    return(flow.frame)
}

get_combi_plot <- function(flow.frames){
    fcs <- concatenateCIPHE(flow.frames)
    fcs <- addDensityFlagCIPHE(fcs)
    # id <- as.vector(which(!is.na(as.vector(pData(fcs@parameters)[,2]))))
    # id <- as.vector(which(as.vector(pData(fcs@parameters)[,2])!="NA"))
    id <- as.vector(which(colnames(fcs)%in%colnames(fcs@description[["SPILL"]])))
    plot_output_hist <- lapply(id,function(x) {
        plotname <- paste0("plotHist_",x[1],"_",x[2])
        plotname2 <- paste0("plotBin_",x[1],"_",x[2])
        
        plot_object <- plotOutput(plotname)
        plot_object <- renderPlot({
            par(mar=c(5,2,2,0))
            hist(exprs(fcs)[sample(c(1:dim(fcs)[1]),round(0.1*dim(fcs)[1])),x],
                 main=pData(fcs@parameters)[x,2], ylab="Counts", xlim=c(-0.5,5), breaks=100,xlab=NULL)  
        },height=180, width=180, outputArgs=list(width="33%"))

        plot_object2 <- plotOutput(plotname2)
        plot_object2 <-  renderPlot({
            par(mar=c(0,3,0,0))
            plotDens(fcs[sample(c(1:dim(fcs)[1]),round(0.1*dim(fcs)[1])),],
                     c(colnames(fcs)[x],"FlagDens"), xlim=c(-0.5,5),ylab="Files",main=NULL,xlab=NULL,yaxt="n")
            axis(2, at=c(1:length(flow.frames)),labels=c(1:length(flow.frames)), las=2)
        },height=180, width=380,outputArgs = list(width="60%"))

        return(list(plot_object, plot_object2))
    })

    plot_output_hist <- unlist(plot_output_hist, recursive = F)

    return(plot_output_hist)
}



## Args
# flow.frame = le fichier fcs
# channels = choix des deux dimensions (labels) ex c("CD4", "CD8)
# nb.clusters = nombre de clusters total
# target= position de la population visée ex:c(3,1)
# translation = Aire de la gate entre 0 et 1 ex:0.5
# quantile = Taille de la barre diagonale entre 0 et 1 ex:0.5
# angle = angle de la 1ere gate HD, HG, BD, BG
# transionntal = nombre de gate à réaliser 2,3 ou 4
# transitionnal.3.gate = position de la seconde gate HD, HG, BD, BG
# th.limit = limites du graphique

## Output
# return.list
# input.fcs = le flow.frame de départ
# output.fcs = toutes les populations recherchés séparés selon la/les gates
# channels = listes des deux dimensions où les populations ont été recherchés
# threshold = formes géométriques des gates (polygones) colnames = channels

transitional.gate.CIPHE <- function(flow.frame, channels, nb.clusters = 2, target = NULL, quantile = NULL, translation = NULL, angle,
                                    transitionnal = 2, transitionnal.3.gate = NULL, th.limit = c(-0.5,4.5)) {
    return.list <- list()
    return.list$input.fcs <- flow.frame
    return.list$output.fcs <- list()
    return.list$threshold <- list()
    x.param <- get.marker.CIPHE(channels[1], flow.frame)[1]
    if(length(x.param) == 0 || is.na(x.param) == TRUE) {x.param <- get.markers.CIPHE(flow.frame)[1]}
    y.param <- get.marker.CIPHE(channels[2], flow.frame)[1]
    if(length(y.param) == 0 || is.na(y.param) == TRUE) {y.param <- get.markers.CIPHE(flow.frame)[1]}
    angle.value <- switch(angle, "HD"= 0.25, "HG" = 0.75, "BG" = 1.25, "BD" = 1.75)
    gate <- openCyto:::.flowClust.2d(flow.frame, channels = c(x.param, y.param), K = nb.clusters, target = target,
                                     quantile = quantile, translation = translation, transitional = TRUE,
                                     transitional_angle = pi * angle.value)
    xmin <- ymin <- th.limit[1]
    xmax <- ymax <- th.limit[2]
    if(transitionnal == 2) {
        return.list$threshold <- list(data.1 = gate@boundaries)
        polygon1 <- flowDensity.CIPHE(flow.frame, channels = c(x.param, y.param), position = c(T,NA), filter = gate@boundaries)
        not.data <- notSubFrame.CIPHE(flow.frame, channels = c(x.param, y.param), gates = polygon1@gates, filter = polygon1@filter)
        return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(not.data))
        colnames(return.list$threshold[[1]]) <- channels
        return.list$channels <- channels
    }
    if(transitionnal == 3) {
        if(angle == "HG" & transitionnal.3.gate == "BG") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,x2,x2,x1,x1,xmax,xmax,ymax,ymax,y2,y1,ymin,ymin,ymax), ncol = 2)
            data.2 <- matrix(c(x2,xmin,xmin,x1,x2,x2,ymax,ymax,y1,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x1,x1,y1,y1,ymin,ymin,y1), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "HG" & transitionnal.3.gate[2] == "HD") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,x2,x2,xmax,xmax,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.2 <- matrix(c(x2,xmin,xmin,x1,x2,x2,ymax,ymax,y1,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x1,xmin,xmin,xmax,xmax,y2,y2,y1,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "HD" & transitionnal.3.gate[2] == "BD") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,x1,x1,x2,xmax,xmax,ymax,ymax,y2,y1,y1,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x2,x2,x1,x1,ymax,ymax,ymin,ymin,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x2,xmax,xmax,y1,y1,ymin,ymin,y1), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "HD" & transitionnal.3.gate[2] == "HG") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,x1,x1,x2,xmax,xmax,ymax,ymax,y2,y1,y1,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x1,x1,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x1,xmin,xmin,xmax,xmax,xmax,y1,y1,y2,y2,ymin,ymin,y1), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BD" & transitionnal.3.gate[2] == "HD") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,x2,x2,xmax,xmax,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.2 <- matrix(c(x2,xmin,xmin,x1,x1,x2,x2,ymax,ymax,ymin,ymin,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x1,x1,xmax,xmax,y2,y2,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BD" & transitionnal.3.gate[2] == "BG") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,xmin,xmin,x1,x2,xmax,xmax,ymax,ymax,y1,y1,y2,y2,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x1,x1,y1,y1,ymin,ymin,y1), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x1,x1,xmax,xmax,y2,y2,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,F),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,F),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BG" & transitionnal.3.gate[2] == "HG") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,x1,x1,x2,x2,xmax,xmax,ymax,ymax,y2,y1,ymin,ymin,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x1,x1,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x2,x2,x1,y2,y2,ymin,ymin,y1,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(T,NA),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BG" & transitionnal.3.gate[2] == "BD") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,xmin,xmin,x1,x2,xmax,xmax,ymax,ymax,y2,y2,y1,y1,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x2,x2,x1,y2,y2,ymin,ymin,y1,y2), ncol = 2)
            data.3 <- matrix(c(xmax,x2,x2,xmax,xmax,y1,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(NA,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param, y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3))
            for(i in c(1:3)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
    }
    if(transitionnal == 4){
        if(angle == "HD") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,x1,x1,x2,xmax,xmax,ymax,ymax,y2,y1,y1,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x1,x1,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x2,x2,x1,y2,y2,ymin,ymin,y1,y2), ncol = 2)
            data.4 <- matrix(c(xmax,x2,x2,xmax,xmax,y1,y1,ymin,ymin,y1), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3, data.4 = data.4)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            polygon4 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                          filter = return.list$threshold[[4]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3), data.4 = get.flowFrame.CIPHE(polygon4))
            for(i in c(1:4)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "HG") {
            x1 <- gate@boundaries[1,1]
            x2 <- gate@boundaries[5,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,x2,x2,xmax,xmax,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.2 <- matrix(c(x2,xmin,xmin,x1,x2,x2,ymax,ymax,y1,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x1,x1,y1,y1,ymin,ymin,y1), ncol = 2)
            data.4 <- matrix(c(xmax,x2,x1,x1,xmax,xmax,y2,y2,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3, data.4 = data.4)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            polygon4 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                          filter = return.list$threshold[[4]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3), data.4 = get.flowFrame.CIPHE(polygon4))
            for(i in c(1:4)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BD") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[5,2]
            y2 <- gate@boundaries[1,2]
            data.1 <- matrix(c(xmax,x2,x2,xmax,xmax,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.2 <- matrix(c(x2,xmin,xmin,x1,x2,x2,ymax,ymax,y1,y1,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x1,x1,y1,y1,ymin,ymin,y1), ncol = 2)
            data.4 <- matrix(c(xmax,x2,x1,x1,xmax,xmax,y2,y2,y1,ymin,ymin,y2), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3, data.4 = data.4)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            polygon4 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                          filter = return.list$threshold[[4]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3), data.4 = get.flowFrame.CIPHE(polygon4))
            for(i in c(1:4)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
        if(angle == "BG") {
            x1 <- gate@boundaries[5,1]
            x2 <- gate@boundaries[1,1]
            y1 <- gate@boundaries[1,2]
            y2 <- gate@boundaries[5,2]
            data.1 <- matrix(c(xmax,x1,x1,x2,xmax,xmax,ymax,ymax,y2,y1,y1,ymax), ncol = 2)
            data.2 <- matrix(c(x1,xmin,xmin,x1,x1,ymax,ymax,y2,y2,ymax), ncol = 2)
            data.3 <- matrix(c(x1,xmin,xmin,x2,x2,x1,y2,y2,ymin,ymin,y1,y2), ncol = 2)
            data.4 <- matrix(c(xmax,x2,x2,xmax,xmax,y1,y1,ymin,ymin,y1), ncol = 2)
            return.list$threshold <- list(data.1 = data.1, data.2 = data.2, data.3 = data.3, data.4 = data.4)
            polygon1 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,T),
                                          filter = return.list$threshold[[1]])
            polygon2 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,T),
                                          filter = return.list$threshold[[2]])
            polygon3 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(F,F),
                                          filter = return.list$threshold[[3]])
            polygon4 <- flowDensity.CIPHE(flow.frame = flow.frame, channels = c(x.param,y.param), position = c(T,F),
                                          filter = return.list$threshold[[4]])
            return.list$output.fcs <- list(data.1 = get.flowFrame.CIPHE(polygon1), data.2 = get.flowFrame.CIPHE(polygon2),
                                           data.3 = get.flowFrame.CIPHE(polygon3), data.4 = get.flowFrame.CIPHE(polygon4))
            for(i in c(1:4)) {colnames(return.list$threshold[[i]]) <- channels}
            return.list$channels <- channels
        }
    }
    return(return.list)
}