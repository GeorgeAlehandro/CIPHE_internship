
FCStSNE.CIPHE <- function(LoaderPATH, #flow.frames
                    ceil = 5000,
                    FNnames, #marker.matrice
                    OutputSuffix = "Out",
                    DotSNE = TRUE,
                    DoOneSENSE = TRUE,
                    Bins = 250,
                    meth = "Rtsne"
                    )
{
    fs <- flowSet(LoaderPATH)
    FcsFileNames <- names(LoaderPATH)#rownames(keyword(fs, "FILENAME"))
    NumBC <- length(fs)  #3
    FFdata <- NULL  #FlowFrame data

    for (FFs in 1:NumBC) {
        # iterate through each FCS file
        FFt <- exprs(fs[[FFs]])
        ## Downsample ##
        if (nrow(FFt) <= ceil)
            FFa <- FFt else FFa <- FFt[sample(nrow(FFt), ceil,replace = FALSE), ]
            colnames(FFa) <- fs[[FFs]]@parameters$name
            empties <- which(is.na(colnames(FFa)) | colnames(FFa) == " ")
            colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
            fs[[FFs]]@parameters$name <- colnames(FFa)
            FFa <- cbind(FFa, rep(FFs, dim(FFa)[1]))
            colnames(FFa)[dim(FFa)[2]] <- "InFile"
            # Concatenate
            FFdata <- rbind(FFdata, FFa)
    }
    message("FCS Files Read")
    keeptable <- FNnames
    # print(keeptable[,2])
    keeprowbool <- sapply(keeptable[,2], function(x) any(x == "Y"))
    keeprows <- subset(keeptable, keeprowbool)
    # print(keeprows[, 1])
    data <- FFdata[, which(colnames(FFdata) %in% keeprows[, 1])]
    # lgcl <- logicleTransform(w = 0.25, t = 16409, m = 4.5, a = 0)
    # ilgcl <- inverseLogicleTransform(trans = lgcl)
    # data1 <- apply(data, 2, lgcl)
    data1 <- data
    # FFdata1 <- apply(FFdata, 2, lgcl)
    FFdata1 <- FFdata
    score <- NULL
    tSNEmat <- NULL
    if (DotSNE) {
        data1[is.na(data1)] <- 0
        if(meth == "Rtsne"){
            data1 <- unique(data1)
            data1[is.na(data1)] <- 0
            tSNEdata3 <- Rtsne(data1, dims = 2)
            tSNEmat <- tSNEdata3$Y
        }
        if(meth == "UMAP"){
            tSNEdata3 <- data.matrix(umap(data1))
            tSNEmat <- tSNEdata3
        }
        # tSNEmat <- tSNEdata3$Y
        colnames(tSNEmat) <- c("tSNE1", "tSNE2")
        temp.data.1 <- data.frame(tSNEmat[, 1], tSNEmat[, 2])
        colnames(temp.data.1) <- c("tSNE1","tSNE2")
        p1 <- ggplot(temp.data.1, aes(x=tSNE1,y=tSNE2))+ geom_point()
        # plot(tSNEmat[, 1], tSNEmat[, 2], pch = ".",xlab = "tSNE1", ylab = "tSNE2", cex = 0.1)
    } else tSNEmat <- NULL


    if (DoOneSENSE) {
        message("Doing oneSENSE")
        Xx1DtSNEmat <- NULL
        #     if (dim(keeptable)[2] == 4) {
        #     for (factor in 2:(dim(keeptable)[2])) {
        #         # for loop from 2 to 4
        #         OneDtSNEname <- colnames(keeptable)[factor]
        #         keeprowbool <- sapply(keeptable[, factor],function(x) any(x == "Y"))
        #         keeprows <- subset(keeptable, keeprowbool)
        #         dataX <- FFdata1[,which(colnames(FFdata1) %in% keeprows[, 1])]
        #         print(dataX)
        #         tSNEdata3 <- Rtsne(as.matrix(dataX), dims = 1, check_duplicates = FALSE)
        #         tSNEmat1 <- tSNEdata3$Y
        #         colnames(tSNEmat1) <- OneDtSNEname
        #         hist(tSNEmat1, 100)
        #         Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
        #     }
        #     p <- scatterplot3d(x = Xx1DtSNEmat[, 1], y = Xx1DtSNEmat[, 2],z = Xx1DtSNEmat[, 3], pch = ".",
        #                     xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1],sep = ""),
        #                     ylab = paste("tSNE2 ", colnames(Xx1DtSNEmat)[2],sep = ""),
        #                     zlab = paste("tSNE3 ", colnames(Xx1DtSNEmat)[3], sep = ""))
        # } else {
        #     # loop from 2 to 4
        for (factor in 2:(dim(keeptable)[2])) {
            OneDtSNEname <- colnames(keeptable)[factor]
            keeprowbool <- sapply(keeptable[, factor],function(x) any(x == "Y"))
            keeprows <- subset(keeptable, keeprowbool)
            # print(keeprows[, 1])
            dataX <- FFdata1[, which(colnames(FFdata1) %in% keeprows[, 1])]
            if(meth == "Rtsne"){
                dataX <- unique(dataX)
                dataX[is.na(dataX)] <- 0
                tSNEdata3 <- Rtsne(dataX, dims = 1, check_duplicates = FALSE)
                tSNEmat1 <- tSNEdata3$Y
            }
            if(meth == "UMAP"){
                tSNEdata3 <- umap(dataX)
                tSNEmat1 <- data.matrix(tSNEdata3[,1])
            }
            print(OneDtSNEname)
            print(dim(tSNEmat1))
            colnames(tSNEmat1) <- OneDtSNEname
            hist(tSNEmat1, 100, plot = F)
            Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
        }
        temp.data.1 <- data.frame(Xx1DtSNEmat)
        colnames(temp.data.1) <- c("tSNE1","tSNE2")

        p2 <- ggplot(temp.data.1, aes(x=tSNE1,y=tSNE2))+ geom_point()
            # p <- plot(Xx1DtSNEmat[, 1], Xx1DtSNEmat[, 2], pch = ".",
            #     xlab = paste("tSNE1 ", colnames(Xx1DtSNEmat)[1], sep = ""),ylab = paste("tSNE2 ",
            #     colnames(Xx1DtSNEmat)[2], sep = ""), cex = 1)

    } else Xx1DtSNEmat <- NULL

    NXx1 <- apply(Xx1DtSNEmat, 2,function(x) ((x - min(x))/(max(x) - min(x))) * 10000)
    score2 <- cbind(score, tSNEmat)

    Nscore <- apply(score2, 2,function(x) ((x - min(x))/(max(x) - min(x))) * 3.7)
    # ilgcl <- inverseLogicleTransform(trans = lgcl)
    NIscore <- Nscore
    colnames(NIscore) <- colnames(score2)
    NIscore <- cbind(NIscore, NXx1)
    # output new FCS files
    new.flow.frames <- c()
    for (FFs in 1:NumBC) {
        newFF <- fs[[1]]
        newBaseData <- FFdata[FFdata[, dim(FFdata)[2]] == FFs, -dim(FFdata)[2]]
        colnames(newBaseData) <- colnames(exprs(newFF))
        exprs(newFF) <- newBaseData
        subsetNIscore <- NIscore[FFdata[, dim(FFdata)[2]] == FFs, ]
        newFF <- cbind2(newFF, subsetNIscore)
        newFF@parameters$desc <- colnames(cbind(newBaseData, subsetNIscore))
        # suppressWarnings(dir.create(paste0(LoaderPATH, "_", OutputSuffix)))
        # print(FcsFileNames)
        # print(FFs)
        BaseFN <- sapply(strsplit(FcsFileNames[FFs], split = "\\."), "[", 1)
        # FNresult <- paste0(LoaderPATH, "_", OutputSuffix,
                        # "/", BaseFN, "_", OutputSuffix, ".fcs")
        newFF@description$"$FIL" <- paste0(BaseFN, "_", OutputSuffix, ".fcs")
        newFF@description$FILENAME <- paste0(BaseFN, "_", OutputSuffix, ".fcs")
        identifier(newFF) <- paste0(BaseFN, "_", OutputSuffix)
        new.flow.frames <- c(new.flow.frames, list(newFF))
        # suppressWarnings(write.FCS(newFF, FNresult))
    }
    return(list(new.flow.frames,p1,p2))
}

OneSmapperFlour.CIPHE <- function(LoaderPATH ,#flow.frames
                            Fname, #markers table
                            Bins = 250,
                            doCoords = FALSE,
                            doFreq = FALSE)
{
    message("Running OneSmapperFlour")
    fs <- flowSet(LoaderPATH)
    FcsFileNames <- rownames(keyword(fs, "FILENAME"))
    FNumBC <- length(fs)
    FFdata <- NULL
    for (FFs in 1:FNumBC) {
        FFa <- exprs(fs[[FFs]])

        # Fixup column names
        colnames(FFa) <- fs[[FFs]]@parameters$desc
        empties <- which(is.na(colnames(FFa)) | colnames(FFa) == " ")
        colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
        fs[[FFs]]@parameters$desc <- colnames(FFa)
        FFa <- cbind(FFa, rep(FFs, dim(FFa)[1]))
        colnames(FFa)[dim(FFa)[2]] <- "InFile"
        # Concatenate
        FFdata <- rbind(FFdata, FFa)
    }
    Xx1DtSNEmat <- NULL
    keeptable <- Fname
    p <- c()
    list.table <- c()
    for (factor in 2:(dim(keeptable)[2])) {
        # edited code
        # print(keeptable)
        keeprowbool <- sapply(keeptable[, factor], function(x) any(x == "Y"))
        keepnames <- keeptable[keeprowbool, 1]
        keeprows <- subset(keeptable, keeprowbool)
        data <- FFdata[, which(colnames(FFdata) %in% keeprows[, 1])]
        data <- cbind(data,
                    FFdata[,
                    which(colnames(FFdata) %in% colnames(keeptable[-1]))])
        data1 <- data
        OneDtSNEname <- colnames(keeptable)[factor]
        dataX <- data1[, which(colnames(data1) %in% keeprows[, 1])]
        tSNEmat1 <- as.matrix(data[, OneDtSNEname])
        colnames(tSNEmat1) <- OneDtSNEname
        hist(tSNEmat1, 100,plot=FALSE)
        Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)

        # Make heatplot annotation


        tSNEbins <- cut(tSNEmat1, breaks = Bins, labels = 1:Bins)
        OneSENSEmed <- apply(dataX, 2,function(x) (tapply(x, tSNEbins, FUN = median)))
        OneSENSEmed[is.na(OneSENSEmed)] <- 0
        # pdfFN <- paste(dirname(LoaderPATH),
        #                 paste0(OneDtSNEname, "_Freq.pdf"),
        #                 sep = .Platform$file.sep)
        # pdf(file = pdfFN, width = 14, height = 6)
        breaks = seq(0, 3, by = 0.005)

        my_palette <-colorRampPalette(c("blue", "white", "red"))(n = length(breaks) - 1)
        png("temp.png",1000,1000)
        hmap <- heatmap.2(t(OneSENSEmed), col = my_palette, breaks = breaks,
                        margins = c(10, 20), Colv = FALSE, dendrogram = "row",
                        cexCol = 1, cexRow = 1, scale = "none", key = TRUE,
                        trace = "none", density.info = c("none"), keysize = 1)
        dev.off()
        file.remove("temp.png")
        hmapclu = t(OneSENSEmed)[hmap$rowInd, hmap$colInd]
        hmapclut <- t(hmapclu)

        testcol <- seq(from = min(Xx1DtSNEmat),to = (max(Xx1DtSNEmat) - (max(Xx1DtSNEmat)/Bins)),by = (max(Xx1DtSNEmat)/Bins))
        
        temp.map <- cbind(hmapclut,testcol)
        
        if (OneDtSNEname == "input") {
            p1 <- plot_ly(z = hmapclu, y = rownames(hmapclu), x = testcol, colors = my_palette,
                type = "heatmap") %>% layout(title = paste(OneDtSNEname,
                "Median heatplot", sep = " "))
        } else {
            p1 <- plot_ly(z = hmapclut, x = colnames(hmapclut), y = testcol,colors = my_palette,
                type = "heatmap") %>% layout(title = paste(OneDtSNEname,
                "Median heatplot", sep = " "), showlegend = FALSE)
        }
        p <- c(p,list(p1))
        list.table <- c(list.table, list(temp.map))
    }

    b <- plot_ly(data.frame(Xx1DtSNEmat),x = Xx1DtSNEmat[, 1],y = Xx1DtSNEmat[, 2],
        marker = list(size = 5,color = 'rgba(0, 0, 0, .2)'),
        type = "scatter",
        symbol = "circle-dot"
    )

    p <- c(p, list(b),list.table)
    
    return(p)
    # suppressWarnings(combined <- subplot(OneSplot, p2, p1,
    #                                     nrows = 2,
    #                                     shareY = TRUE,
    #                                     shareX = TRUE))
    # combined
    # export(combined, file = paste(dirname(LoaderPATH),
    #         "groupone.png",
    #         sep = .Platform$file.sep))
    # browseURL(paste(dirname(LoaderPATH),
    #                 "groupone.png",
    #                 sep = .Platform$file.sep))
}

getParameters <- function(rawFCSdir) 
{
    fcsFile <- list.files(path = rawFCSdir,
                            pattern = ".fcs$", full.names = TRUE)
    fcs <- suppressWarnings(read.FCS(fcsFile[1]))
    pd <- fcs@parameters@data
    markers <- paste("<", pd$name, ">:", pd$desc, sep = "")
    markers_desc <- pd$desc
    gParam = list(markers = markers, markers_desc = markers_desc)

    if (length(markers) == 0) {
        stop("No markers found in the FCS file!")
    }
    return(gParam)
}

OneSmapperFreq1 <- function(LoaderPATH = "fcs_Out") 
{
    fs <- read.flowSet(path = LoaderPATH, pattern = ".fcs$")
    FcsFileNames <- rownames(keyword(fs, "FILENAME"))
    FNumBC <- length(fs)
    FFdata <- NULL
    for (FFs in 1:FNumBC) {
        FFa <- exprs(fs[[FFs]])
        # Fixup column names
        colnames(FFa) <- fs[[FFs]]@parameters$desc
        empties <- which(is.na(colnames(FFa)) | colnames(FFa) == " ")
        colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
        fs[[FFs]]@parameters$desc <- colnames(FFa)
        # Add file label
        FFa <- cbind(FFa, rep(FFs, dim(FFa)[1]))
        colnames(FFa)[dim(FFa)[2]] <- "InFile"
        # Concatenate
        FFdata <- rbind(FFdata, FFa)
    }
    return(FFdata)
}

getCoords <- function(LoaderPATH = LoaderPATH, FFdata = FFdata)
{
    lgcl <- logicleTransform(w = 0.25, t = 16409, m = 4.5, a = 0)
    gckeeptable <- read.csv(paste(dirname(LoaderPATH),
                                "names.csv",
                                sep = .Platform$file.sep))
    gckeeprowbool <- apply(gckeeptable[, c(2, 3)],
                                1,
                                function(x) any(x == "Y"))
    gckeepnames <- gckeeptable[gckeeprowbool, 1]
    gckeeprows <- subset(gckeeptable, gckeeprowbool)
    gcdata <- FFdata[, which(colnames(FFdata) %in% gckeeprows[, 1])]
    gcdata <- cbind(gcdata,
                FFdata[, which(
                    colnames(FFdata) %in% colnames(gckeeptable[-1]))])
    data1 <- apply(gcdata, 2, lgcl)
    coordsVar <- list(keepnames = gckeepnames, data1 = data1)
    return(coordsVar)
}

OneSmapperFreq2 <- function(LoaderPATH = "fcs", Bins = 250, FFdata)
{
    lgcl <- logicleTransform(w = 0.05, t = 16409, m = 4.5, a = 0)
    Xx1DtSNEmat <- NULL
    keeptable <- read.csv(paste(dirname(LoaderPATH),
                                "names.csv",
                                sep = .Platform$file.sep))
    for (factor in 2:(dim(keeptable)[2])) {
    keeprowbool <- sapply(keeptable[, factor],
                            function(x) any(x == "Y"))
    keepnames <- keeptable[keeprowbool, 1]
    keeprows <- subset(keeptable, keeprowbool)
    data <- FFdata[, which(colnames(FFdata) %in% keeprows[, 1])]
    data <- cbind(data,
            FFdata[, which(colnames(FFdata) %in% colnames(keeptable[-1]))])
    data1 <- apply(data, 2, lgcl)  #logicle transform
    OneDtSNEname <- colnames(keeptable)[factor]
    keeprowbool <- sapply(keeptable[, factor],
                        function(x) any(x == "Y"))
    keepnames <- keeptable[keeprowbool, 1]
    keeprows <- subset(keeptable, keeprowbool)
    dataX <- data1[, which(colnames(data1) %in% keeprows[, 1])]
    tSNEmat1 <- as.matrix(data[, OneDtSNEname])
    colnames(tSNEmat1) <- OneDtSNEname
    hist(tSNEmat1, 100)
    Xx1DtSNEmat <- cbind(Xx1DtSNEmat, tSNEmat1)
    # Make heatplot annotation

    dataGNorm <- dataX
    tSNEbins <- cut(tSNEmat1, breaks = Bins, labels = 1:Bins)
    OneSENSEpp <- matrix(dataX, nrow = Bins, ncol = dim(dataX)[2])
    colnames(OneSENSEpp) <- colnames(dataX)
    Coords <- read.csv(paste(dirname(LoaderPATH),
                            "Coords.csv",
                            sep = .Platform$file.sep))
    for (pname in colnames(dataX)) {
        overthresh <- function(group) {
        # what is this group
        percpos <- (
            sum(group > Coords[which(Coords$X == pname),
                                2])/length(group)) * 100
        # print(percpos)
        return(percpos)
    }

    OneSENSEpp[, pname] <- tapply(dataX[, pname], tSNEbins, FUN = overthresh)
    }

    OneSENSEpp[is.na(OneSENSEpp)] <- 0


    pdfFN <- paste(dirname(LoaderPATH),
                    paste0(OneDtSNEname, "_Freq.pdf"),
                    sep = .Platform$file.sep)
    pdf(file = pdfFN, width = 14, height = 6)
    breaks = seq(0, 100, by = 0.05)
    my_palette <- colorRampPalette(c("blue", "white",
                                    "red"))(n = length(breaks) - 1)



    fhmap <- heatmap.2(t(OneSENSEpp),
                        col = my_palette,
                        breaks = breaks,
                        margins = c(10, 20),
                        Colv = FALSE,
                        dendrogram = "row",
                        cexCol = 1,
                        cexRow = 1,
                        scale = "none",
                        key = TRUE,
                        trace = "none",
                        density.info = c("none"),
                        keysize = 1)
    fhmapclu = t(OneSENSEpp)[fhmap$rowInd, fhmap$colInd]
    fhmapclut <- t(fhmapclu)
    ftestcol <- seq(from = min(Xx1DtSNEmat),
                    to = (max(Xx1DtSNEmat) - (max(Xx1DtSNEmat)/Bins)),
                    by = (max(Xx1DtSNEmat)/Bins))
    if (OneDtSNEname == "input") {
        suppressWarnings(d1 <- plot_ly(z = fhmapclu,
                                    y = rownames(fhmapclu),
                                    x = ftestcol,
                                    colors = my_palette,
                                    type = "heatmap") %>%
                        layout(title = paste(OneDtSNEname,
                                            "Frequency heatplot", sep = " ")))
        d1
    } else {
        suppressWarnings(d2 <- plot_ly(z = fhmapclut,
                                    x = colnames(fhmapclut),
                                    y = ftestcol,
                                    colors = my_palette,
                                    type = "heatmap") %>%
                        layout(title = paste(OneDtSNEname,
                                            "Frequency heatplot", sep = " ")))
        d2
    }
    dev.off()

    }
    OneSplot <- plot_ly(data.frame(Xx1DtSNEmat),
                        x = Xx1DtSNEmat[, 1],
                        y = Xx1DtSNEmat[, 2],
                        marker = list(size = 5,
                                    color = 'rgba(0, 0, 0, .2)'
                        ),
                        type = "scatter",
                        symbol = "circle-dot")
    suppressWarnings(fcombined <- subplot(OneSplot, d2, d1, nrows = 2,
                                shareY = TRUE, shareX = TRUE))
    fcombined
    export(fcombined, file = paste(dirname(LoaderPATH),
                                "grouptwo.png",
                                sep = .Platform$file.sep))

    browseURL(paste(dirname(LoaderPATH),
                "grouptwo.png",
                sep = .Platform$file.sep))
}
