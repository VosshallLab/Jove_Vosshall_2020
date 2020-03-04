vnorm <- function(ligandInfo,mean_GCaMP,pumps=NULL,exportXLSX=NULL,windowBefore=5,windowAfter=0,voltage=5,matchByTime=FALSE){
  # ligandInfo <- "/Users/tcarroll/Downloads/Re__summary_of_today's_meeting/19.9.10.F2_biopen003_events_ligand_info.csv"
  # mean_GCaMP <- "/Users/tcarroll/Downloads/Re__summary_of_today\'s_meeting/19.9.10.F2_biopen003_mean_GCaMP.txt"
  # pumps = NULL
  # exportXLSX="~/Downloads/veronicaRes1.csv"
  # windowBefore=5
  # windowAfter=0
  # voltage=5
  # matchByTime=FALSE
  require(magrittr)
  if(matchByTime){
    require(lubridate)
    ligandInfo <- read.delim(ligandInfo,sep=",")
    ligandInfo <-  ligandInfo[,grepl("^PUMP|Frame|Time",colnames(ligandInfo))]
    milliseconds_1 <- period_to_seconds(ms(gsub("\\..*","",ligandInfo$Time)))*10^4
    milliseconds_ligandInfo <- milliseconds_1+signif(as.numeric(gsub(".*\\.","",ligandInfo$Time)),3)/10
    ligandInfo$Time <- milliseconds_ligandInfo
    
    mean_GCaMP <- read.delim(mean_GCaMP,comment.char = ";",sep="\t")
    mean_GCaMP_values <- mean_GCaMP[,grepl("^X",colnames(mean_GCaMP))]
    milliseconds_2 <- as.numeric(gsub("\\..*","",mean_GCaMP[,grepl("^Time",colnames(mean_GCaMP))]))*10^4
    milliseconds_mean_GCaMP_values <- milliseconds_2+as.numeric(gsub(".*\\.","",formatC(mean_GCaMP[,grepl("^Time",colnames(mean_GCaMP))],format = "f",digits = 3)))
    mean_GCaMP_values$Time  <- milliseconds_mean_GCaMP_values
    mean_GCaMP_And_ligandInfo <- merge(mean_GCaMP_values,ligandInfo,by="Time",all=TRUE)
  }else{
    ligandInfo <- read.delim(ligandInfo,sep=",")
    ligandInfo <- eventsToPumps(ligandInfo)
    if(any(colnames(ligandInfo) %in% "Frame")){
      ligandInfo <-  ligandInfo[,grepl("^PUMP|Frame|Time",colnames(ligandInfo))]
    }else{
      ligandInfo <-  ligandInfo[,grepl("^PUMP|Index|Time",colnames(ligandInfo))]
      colnames(ligandInfo)[colnames(ligandInfo) %in% "Index"] <- "Frame"
    }
    
    colnames(ligandInfo)[grep("PUMP",colnames(ligandInfo))] <- gsub("V|\\.","",colnames(ligandInfo)[grep("PUMP",colnames(ligandInfo))])
    mean_GCaMP <- read.delim(mean_GCaMP,comment.char = ";",sep="\t")
    mean_GCaMP_values <- mean_GCaMP[,grepl("^X",colnames(mean_GCaMP))]
    mean_GCaMP_And_ligandInfo <- cbind(ligandInfo,mean_GCaMP_values)
  }
  
  stopifnot(nrow(mean_GCaMP_And_ligandInfo) == nrow(ligandInfo))
  
  colnames(mean_GCaMP_And_ligandInfo)[grep("^Time",colnames(mean_GCaMP_And_ligandInfo))] <- "Time"
  
  totalPumpColumns <- colnames(ligandInfo)[grep("PUMP",colnames(ligandInfo))]
  if(!is.null(pumps)){
    pumps <- paste0("PUMP",pumps)
  }else{
    pumps <- totalPumpColumns
  }
  
  resA <- list()
  namesA <- c()
  for(p in 1:length(pumps)){
    blocks <- mean_GCaMP_And_ligandInfo[,pumps[p]] == voltage
    if(any(blocks)){
      result <- diff(blocks)
      if(result[result != 0][1] == -1){
        resultStart <- c(1,which(result == 1)+1)
      }else{
        resultStart <- which(result == 1)+1
      }
      if(tail(result[result != 0],1) == 1){
        resultEnd <- c(which(result == -1),length(result))
      }else{
        resultEnd <- which(result == -1)
      }
      
      for(i in 1:length(resultStart)){
        
        baseline <- colMeans(mean_GCaMP_And_ligandInfo[seq(resultStart[i]-windowBefore,resultStart[i]-1),
                                                       grep("^X\\.*",colnames(mean_GCaMP_And_ligandInfo))])
        
        baselineFrames <- mean_GCaMP_And_ligandInfo[seq(resultStart[i]-windowBefore,resultStart[i]-1),
                                                    "Frame"]
        tobaseLineFrame <- paste0("baselineFrame_",min(baselineFrames),"_to_Frame",max(baselineFrames))
        
        toNormTable <- mean_GCaMP_And_ligandInfo[seq(resultStart[i],resultEnd[i]+windowAfter),
                                                 grep("^X\\.*",colnames(mean_GCaMP_And_ligandInfo))]
        
        toNormTable_WOI <- sapply(colnames(toNormTable),function(x)(toNormTable[,x]-baseline[x])/baseline[x])
        colnames(toNormTable_WOI) <- gsub("^X\\.","N_",colnames(toNormTable))
        
        baselineDF <- data.frame(matrix(baseline,ncol=length(baseline),byrow = TRUE,
                                        nrow=nrow(toNormTable))
        )
        colnames(baselineDF) <- gsub("^X\\.","Baseline_",names(baseline))
        
        tf_df <- mean_GCaMP_And_ligandInfo[seq(resultStart[i],resultEnd[i]+windowAfter),c("Time","Frame")]
        
        colnames(toNormTable) <- gsub("^X\\.","Raw_",colnames(toNormTable))
        pumpDF <- data.frame(pumpTested=rep(pumps[p],nrow(toNormTable)),Opening=paste0("Open_",i),
                             baseLineFrame <- rep(tobaseLineFrame,nrow(toNormTable)),
                             frameTested=rep(paste0("Frame_",min(tf_df$Frame),"_to_Frame",max(tf_df$Frame)),nrow(toNormTable)))
        
        newDF <- cbind(pumpDF,tf_df,toNormTable_WOI,baselineDF,toNormTable)
        newDF <- newDF[order(newDF$Frame,decreasing = FALSE),]
        colnames(newDF)[3] <- "baselineTested"
        
        
        resA <- c(resA,list(newDF))
        namesA <- c(namesA,paste0(pumps[p],"__","Frame_",min(newDF$Frame),"_to_Frame",max(newDF$Frame)))
      }
    }
    
  }
  names(resA) <- namesA
  all <- list(do.call(rbind,resA))
  names(all) <- "All_pumps_and_frames"
  all <- c(all,resA)
  if(!is.null(exportXLSX)){
    if(grepl("\\.xlsx$",exportXLSX)){
      require(rio)
      rio::export(all,file = exportXLSX)
    }else{
      write.table(all[[1]],file = exportXLSX,sep=",",row.names = FALSE,quote=FALSE)
    }
  }
  all[[1]] <- all[[1]][order(all[[1]]$Frame,decreasing=FALSE),]
  return(all)
}

normPerColumn <- function(column,before=1,after=1){
  maxPos <- which.max(column)
  mean(column[seq(maxPos-before,maxPos+after)],na.rm=TRUE)
}
maxPerColumn <- function(column){
  which.max(column)
}

calcMeans <- function(files,pumps=NULL,openings=NULL,avWindowBefore=1,avWindowAfter=1,outdir=getwd(),prefix=""){
  Movies <- lapply(files,read.delim,sep=",")
  dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
  names(Movies) <- gsub("\\.csv","",basename(files))
  if(is.null(pumps)) stop("Must Specify a/some pump/s of interest")
  if(is.null(openings)){
    message("Setting number of opening to analyse to 1")
    openings <- 1
  }
  if(length(openings) == 1){
    openings <- rep(openings,length(pumps))
  }
  frameMovieList <- NULL
  for(m in 1:length(Movies)){
    MOI <- Movies[[m]]
    frameList <- NULL
    for(p in 1:length(pumps)){
      selectedFrame <- MOI[MOI$pumpTested == paste0("PUMP",pumps[p]) & MOI$Opening == paste0("Open_",openings[p]),]
      signalColumns <- selectedFrame[,grepl("^N_\\d.*",colnames(selectedFrame))]
      normSignalColumns <- apply(signalColumns,2,normPerColumn,before=avWindowBefore,after=avWindowAfter)
      whereMaxSignalColumns <- apply(signalColumns,2,maxPerColumn)
      whereMaxSignalFrame <- t(data.frame(whereMaxSignalColumns))
      colnames(whereMaxSignalFrame) <- paste0("WhereMax_",colnames(whereMaxSignalFrame))
      normSignalFrame <- t(data.frame(normSignalColumns))
      newSelectedFrame <- cbind(data.frame(Movie=names(Movies)[m],selectedFrame[1,1:4,drop=FALSE]),normSignalFrame,whereMaxSignalFrame)
      frameList <- c(frameList,list(newSelectedFrame))
    }
    movieFrame <- do.call(rbind,frameList)
    frameMovieList <- c(frameMovieList,list(movieFrame))
  }
  names(frameMovieList) <- names(Movies)
  allMovieFrame <- do.call(rbind,frameMovieList)
  allMovieFrame <- allMovieFrame[order(allMovieFrame$pumpTested,allMovieFrame$Opening,allMovieFrame$Movie),]
  resultsList <- c(list(allConcatenated=allMovieFrame),frameMovieList)
  
  
  require(dplyr)
  require(magrittr)
  summarisedPumpFrame <- allMovieFrame %>% 
    group_by(pumpTested,Opening) %>% 
    summarise_at(vars(starts_with("N_")),list(Mean=mean,Stdev=var)) %>% 
    ungroup %>% 
    as.data.frame
  summarisedPumpFrame$PercentOfCellsAbove0.2=summarisedPumpFrame %>% dplyr::select(ends_with("_Mean")) %>% apply(1,function(x) 100*(sum(x > 0.2)/length(x)))
  resultsList <- c(list(SummarisedAcrossMovies=summarisedPumpFrame),resultsList)
  
  for(i in 1:length(resultsList)){
    outfileMain <- file.path(outdir,paste0(prefix,names(resultsList)[i],"_perCellPumpMean.csv"))
    write.table(resultsList[[i]],file=outfileMain,sep = ",",quote = FALSE,row.names = FALSE)
  }
  
  rio::export(resultsList,file = file.path(outdir,paste0(prefix,"_perCellPumpMean.xlsx")),format = "xlsx")
  
  
  return(resultsList)
}

plotAsHeatmap  <- function(input,
                           annotation,
                           outputfile=NULL,
                           boxWidth=40,
                           boxHeight=40,
                           minValue=0,maxValue=NULL,
                           width=10,
                           height=3,
                           colors=c("Blue","Red"),
                           maxInScale=NULL,
                           minInScale=NULL,
                           log2Transform=FALSE,
                           pseudoResponseforZero=0.0000001,
                           clusterCol=FALSE,
                           neuron_kclusters=NA,
                           neuron_cutclusters=NA,
                           cluster_colours=NULL,
                           ...
){
  require(gplots)

  # input <- list(myComMovies,myComMovies)
  # input <- myComMovies
  # annotation=annotation4
  # outputfile="myTest_Joined.pdf"
  # width = 20
  # height = 3
  # # colorMin = "cornflowerblue"
  # # colorMax = "firebrick3"
  # colors=c("Blue","Red")
  # maxInScale=NULL
  # minInScale=NULL
  # boxSize=40
  # boxWidth=40
  # boxHeight=40
  # minValue=0
  # maxValue=NULL
  # log2Transform=FALSE
  # pseudoResponseforZero=0.0000001
  # clusterCol=FALSE

  require(pheatmap)
  colAnnotation <- NULL
  
  if(any(names(input) %in% "allConcatenated")){
    myHeatmaps <- input$SummarisedAcrossMovies %>% dplyr::select(ends_with("Mean")) %>% as.matrix
    rownames(myHeatmaps) <- namesToReplace <- paste0(input$SummarisedAcrossMovies$pumpTested,":",input$SummarisedAcrossMovies$Opening)
    namesToRelabel <- unlist(annotation)
    names(namesToRelabel) <- gsub("annotation\\.","",names(namesToRelabel))
    rownames(myHeatmaps) <- names(namesToRelabel[match(rownames(myHeatmaps),namesToRelabel)])
    myHeatmaps <- myHeatmaps
    colnames(myHeatmaps) <- paste0("Movie_1","_","Neuron_",seq(1,ncol(myHeatmaps)))
    colAnnotation <- data.frame(Movie=rep("Movie_1",ncol(myHeatmaps)),row.names = colnames(myHeatmaps))
    
  }else{
    tempMMat <- list()
    for(i in 1:length(input)){
      myHeatmaps <- input[[i]]$SummarisedAcrossMovies %>% dplyr::select(ends_with("Mean")) %>% as.matrix
      rownames(myHeatmaps) <- namesToReplace <- paste0(input[[i]]$SummarisedAcrossMovies$pumpTested,":",input[[i]]$SummarisedAcrossMovies$Opening)
      namesToRelabel <- unlist(annotation)
      names(namesToRelabel) <- gsub("annotation\\.","",names(namesToRelabel))
      rownames(myHeatmaps) <- names(namesToRelabel[match(rownames(myHeatmaps),namesToRelabel)])
      colnames(myHeatmaps) <- paste0("Movie_",i,"_","Neuron_",seq(1,ncol(myHeatmaps)))
      tempColAnnotation <- data.frame(Movie=rep(paste0("Movie_",i),ncol(myHeatmaps)),row.names = colnames(myHeatmaps))
      tempMMat[[i]] <- myHeatmaps
      if(i == 1){
        fullHeatmap <- myHeatmaps
        colAnnotation <- tempColAnnotation
      }else{
        fullHeatmap <- merge(fullHeatmap,myHeatmaps,by=0)
        tempNames <- fullHeatmap[,1]
        fullHeatmap <- fullHeatmap[,-1]
        rownames(fullHeatmap) <- tempNames
        colAnnotation <- rbind(colAnnotation,tempColAnnotation)
      }
    }
    myHeatmaps <- fullHeatmap
  }
  
  blks <- list()
  for(i in 1:length(annotation)){
    blks[[i]] <- myHeatmaps[rownames(myHeatmaps) %in% names(annotation[[i]]),,drop=FALSE]
  }
  tempMat <- do.call(rbind,blks)
  tempMat <- tempMat[match(names(namesToRelabel),rownames(tempMat)),]
  gapsToPlace <- unlist(lapply(blks,nrow))
  gapsToPlace <-   gapsToPlace[-length(gapsToPlace)]
  tempMatOriginal <- tempMat
  
  if(is.null(minInScale)) minInScale <- 0
  if(is.null(maxInScale)) maxInScale <- ifelse(is.null(maxValue),max(tempMat),min(maxValue,max(tempMat)))
  
  tempMat[tempMat <= max(minValue,minInScale)] <- max(minValue,minInScale)
  if(is.null(maxValue)) maxValue <- maxInScale
  tempMat[tempMat > min(maxValue,maxInScale)] <-  min(maxValue,maxInScale)
  
  
  breaks = seq(minInScale, maxInScale,length.out = 100)
  scaleBreaks <- breaks[c(1,round(length(breaks)/4),round(length(breaks)/2),(round(length(breaks)/4))*3,length(breaks))]
  scaleBreakValues <- scaleBreaks
  
  if(length(colors) == 1)colors <- c("White",colors)
  colors <- colorRampPalette(colors)(length(breaks))
  
  if(log2Transform){
    tempMat[tempMat < pseudoResponseforZero] <- pseudoResponseforZero
    tempMat <- log2(tempMat)
    breaks = seq(ifelse(minInScale < pseudoResponseforZero,log2(pseudoResponseforZero),log2(minInScale)), log2(maxInScale),length.out = 100)
    scaleBreaks <- breaks[c(1,round(length(breaks)/4),round(length(breaks)/2),(round(length(breaks)/4))*3,length(breaks))]
    scaleBreakValues <- 2^scaleBreaks 
    
  }
  
  if(!clusterCol){
    labels_col <- c(1,rep("",ncol(myHeatmaps)-2),ncol(myHeatmaps))
  }else{
    labels_col <- seq(1,ncol(myHeatmaps))
  }

  
  if(!is.null(outputfile)){
    pdf(outputfile,width=width,height=height)
    require(tidyverse)
    require(stats)
    if(!is.na(neuron_cutclusters)){
      cutMembers <- as.data.frame(cutree(hclust(dist(t(tempMat))),neuron_cutclusters))
      colnames(cutMembers) <- "Cluster"
      cutMembers$Cluster <- factor(cutMembers$Cluster)
      colAnnotation <- merge(colAnnotation,cutMembers,by=0) %>% column_to_rownames("Row.names")
      if(!is.null(cluster_colours)){
        ann_colors <- list(Cluster=cluster_colours[names(cluster_colours) %in% as.character(seq(1,neuron_cutclusters))])
      }else{
        ann_colors <- NULL
      }
    }

    l <- pheatmap(tempMat,
             cellwidth = boxWidth,
             cellheight = boxHeight,
             cluster_rows = FALSE,
             cluster_cols = clusterCol,
             labels_col=labels_col,
             gaps_row = cumsum(gapsToPlace),
             scale="none",angle_col="0",
             colors,
             annotation_col = colAnnotation,
             breaks=breaks,
             legend_breaks = scaleBreaks,
             legend_labels = signif(scaleBreakValues,2),
             cutree_cols = neuron_cutclusters,
             ...)
    dev.off()
  }else{
    l <- pheatmap(tempMat,
             cellwidth = boxWidth,
             cellheight = boxHeight,
             cluster_rows = FALSE,
             cluster_cols = clusterCol,
             labels_col=labels_col,
             gaps_row = cumsum(gapsToPlace),
             scale="none",angle_col="0",
             colors,
             breaks=breaks,
             legend_breaks = scaleBreaks,
             legend_labels = signif(scaleBreakValues,2),
             annotation_col = colAnnotation,
             cutree_cols = neuron_cutclusters,
             ...)
  }
  # arrows(1.4, -1, 1.4, 0, xpd = TRUE)
  if(!is.na(neuron_cutclusters)){
    myMembers <- as.data.frame(cutree(l$tree_col,neuron_cutclusters))
    colnames(myMembers) <- "Cluster"
    myMembers$Cluster <- factor(myMembers$Cluster)
    ForN <- myMembers %>% rownames_to_column %>% mutate(Movie=gsub("_Neuron.*","",rowname)) %>% group_by(Cluster,Movie) %>% summarise(NofMovie=table(Movie)) %>% as.data.frame
  }else{
    myMembers <- NULL
    ForN <- NULL
  }
  
  return(list(CllusterInfo=list(tree_row=l$tree_row,tree_col=l$tree_col,kmeans=l$kmeans,gtable=l$gtable,CutMembers=myMembers,MoviesPerCluster=ForN),Values=tempMat,OriginalValues=tempMatOriginal,Min=min(tempMat),Max=max(tempMat)))
  
  # Heatmap(myHeatmaps,cluster_rows = FALSE,cluster_columns = FALSE,row_split=rep(c("Pump1", "Pump2","Pump3"), 9), row_gap = unit(5, "mm"))
}

joinTwoMovies <- function(movieA,movieB){
  newMovie <- list()
  newMovie$SummarisedAcrossMovies <- rbind(movieA$SummarisedAcrossMovies,movieB$SummarisedAcrossMovies)
  newMovie$allConcatenated <- rbind(movieA$allConcatenated,movieB$allConcatenated)
  newMovie <- c(newMovie,movieA[!names(movieA) %in% c("SummarisedAcrossMovies","allConcatenated")],movieB[!names(movieB) %in% c("SummarisedAcrossMovies","allConcatenated")])
  return(newMovie)
}


joinMovies <- function(x,...){
  myMoviesToJoin <- list(x,...)
  do.call(joinTwoMovies,myMoviesToJoin)
}



eventsToPumps <- function(input){
  # input <- read.delim("~/Downloads/19.9.10.F2_biopen003_events_ligand_info.csv",sep=",")
  if(any(colnames(input) %in% "Events")){
    myEvents <- unique(input$Events)
    myEvents <- as.vector(myEvents[myEvents != ""])
    names(myEvents) <- myEvents %>%  gsub(" :.*","",.)
    myEvents <- myEvents[order(names(myEvents))]
    evAll <- c(which(grepl("User",input$Events)),length(input$Events))
    eventList <- list()
    for(i in 1:length(myEvents)){
      temp <- rep(0,length(input$Events))
      evOI <- grep(myEvents[i],input$Events)
      for(l in 1:length(evOI)){
        evOII <- evOI[l]
        endOfevOI <- evAll[which(evAll == evOII)+1]
        temp[seq(evOII,endOfevOI)] <- 5
      }
      eventList[[i]] <- temp
    }
    newPumps <- do.call(cbind,eventList)
    colnames(newPumps) <- names(myEvents)
    fullMat <- matrix(nrow=nrow(newPumps),ncol=5,dimnames=list(NULL,paste0("PUMP",17:21)))
    for(i in 1:ncol(fullMat)){
      if(any(colnames(newPumps) %in% paste0("User ",i))){
        fullMat[,i] <- newPumps[,colnames(newPumps) %in% paste0("User ",i)]
      }else{
        fullMat[,i] <- rep(0,nrow(fullMat))
      }
    }  
    # newPumps <- newPumps[,!grepl("OFF",myEvents)]
    # colnames(newPumps) <- paste0("PUMP",
    #                              seq(sum(grepl("^PUMP",colnames(input)))+1,sum(grepl("^PUMP",colnames(input)))+ncol(newPumps)))
    input <- cbind(input,fullMat)
    input <- input[!grepl("User",input$Events),]
    input[,-grep("Event",colnames(input))]
  }
  return(input)
}


plotAsPCA <- function(heatOut){
  graphics.off()
  
  require(gridExtra)
  values <- heatOut$Values
  myPcs <- prcomp(t(values))$x %>% as.data.frame %>% rownames_to_column(var = "Neuron")
  if(is.null(heatOut$CllusterInfo$CutMembers)){
    p1_2 <- ggplot(myPcs,aes(x=PC1,y=PC2))+geom_point(size=4)+theme_bw()
    p1_3 <- ggplot(myPcs,aes(x=PC1,y=PC3))+geom_point(size=4)+theme_bw()
    p2_3 <- ggplot(myPcs,aes(x=PC2,y=PC3))+geom_point(size=4)+theme_bw()
  }else{
    memberShip <- heatOut$CllusterInfo$CutMembers %>%  as.data.frame %>% rownames_to_column(var = "Neuron")
    myPcs <- merge(memberShip,myPcs,by="Neuron")
    p1_2 <- ggplot(myPcs,aes(x=PC1,y=PC2,colour=Cluster))+geom_point(size=4)+theme_bw()
    p1_3 <- ggplot(myPcs,aes(x=PC1,y=PC3,colour=Cluster))+geom_point(size=4)+theme_bw()
    p2_3 <- ggplot(myPcs,aes(x=PC2,y=PC3,colour=Cluster))+geom_point(size=4)+theme_bw()
  }
  p <- grid.arrange(p1_2,p1_3,p2_3)
  p
  return(list(p,list(p1_2,p1_3,p2_3)))
}

plotClusters <- function(heatOut,filterTreatment=NULL,filterCluster=NULL,rawValues=TRUE,colourMap=NULL){
  require(tidyverse)
  require(ggplot2)
  require(gridExtra)
  graphics.off()
  
  if(rawValues){
    values <- heatOut$OriginalValues %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }else{
    values <- heatOut$Values %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }
  memberShip <- heatOut$CllusterInfo$CutMembers %>%  as.data.frame %>% rownames_to_column(var = "Neuron")
  TableToPlot <- merge(memberShip,gather(values,key = "Neuron",value= "Response",-Treatment),by="Neuron")
  TableToPlot$Cluster <- factor(TableToPlot$Cluster)
  if(!is.null(filterTreatment)){
    TableToPlot <- TableToPlot %>% dplyr::filter(Treatment %in% filterTreatment)
    TableToPlot$Treatment <- factor(as.vector(TableToPlot$Treatment),levels = filterTreatment)
  }
  if(!is.null(filterCluster)){
    TableToPlot <- TableToPlot %>% dplyr::filter(Cluster %in% filterCluster)
    TableToPlot$Cluster <- factor(as.vector(TableToPlot$Cluster),levels = filterCluster)
  }
  p <- ggplot(TableToPlot,aes(y=Response,x=Cluster,fill=Cluster))+geom_boxplot()+facet_grid(~Treatment)+theme_bw()
  p2 <- ggplot(TableToPlot,aes(y=Response,x=Treatment,fill=Treatment))+geom_boxplot()+facet_grid(~Cluster)+theme_bw()
  if(!is.null(colourMap)){
    p <- p+scale_fill_manual(values=colourMap)
    p2 <- p2+scale_fill_manual(values=colourMap)
    
  }
  require(gridExtra)
  pA <- grid.arrange(p,p2)
  pA
  return(list(pA,list(p,p2)))
}

plotScatterClusters <- function(heatOut,filterTreatment=NULL,filterCluster=NULL,rawValues=TRUE,plotByCluster=TRUE,colourMap=NULL){
  # heatOut
  # filterTreatment=NULL
  # filterCluster=NULL
  # rawValues=TRUE
  # plotByCluster=TRUE
  require(tidyverse)
  require(ggplot2)
  require(gridExtra)
  graphics.off()

  if(rawValues){
    values <- heatOut$OriginalValues %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }else{
    values <- heatOut$Values %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }
  
  MaxAndMin <- quantile(as.numeric(as.matrix(values[,-1])),probs = c(0,1))
  
  memberShip <- heatOut$CllusterInfo$CutMembers %>%  as.data.frame %>% rownames_to_column(var = "Neuron")
  TableToPlot <- merge(memberShip,gather(values,key = "Neuron",value= "Response",-Treatment),by="Neuron")
  TableToPlot$Cluster <- factor(TableToPlot$Cluster)
  TableToPlot$Treatment <- factor(TableToPlot$Treatment)
  if(!is.null(filterTreatment)){
    TableToPlot <- TableToPlot %>% dplyr::filter(Treatment %in% filterTreatment)
    TableToPlot$Treatment <- factor(as.vector(TableToPlot$Treatment),levels = filterTreatment)
  }
  if(!is.null(filterCluster)){
    TableToPlot <- TableToPlot %>% dplyr::filter(Cluster %in% filterCluster)
    TableToPlot$Cluster <- factor(as.vector(TableToPlot$Cluster),levels = filterCluster)
  }
  if(plotByCluster){
  TestConditions <- levels(TableToPlot$Treatment)
  dd <- list()
  for(i in 1:length(TestConditions)){

    myCond <- TestConditions[!TestConditions %in% TestConditions[i]]
    myOtherCond <- TestConditions[i]
    
    tempo <- TableToPlot %>% spread(Treatment,Response) %>% gather_("Treatment","Response",gather_cols=as.vector(myCond))
    if(is.null(filterTreatment)){
      filterTreatment2  <- unique(as.vector(myCond))
    }else{
      filterTreatment2  <-  filterTreatment[filterTreatment %in% unique(as.vector(myCond))]
    }
    tempo$Treatment <-   factor(tempo$Treatment,levels = filterTreatment2)
    dd[[i]] <-  tempo %>% 
    ggplot(aes_string(y=as.vector(myOtherCond),x="Response",colour="Treatment"))+geom_point(size=2,alpha=0.8)+theme_bw()+
      ylim(MaxAndMin[1],MaxAndMin[2])+xlim(MaxAndMin[1],MaxAndMin[2])
    if(length(unique(TableToPlot$Cluster)) > 1){
      dd[[i]] <-  dd[[i]] +facet_grid(~Cluster)
    }
    
  }
  }
  if(!plotByCluster){
    TestConditions <- levels(TableToPlot$Cluster)
    dd <- list()
    for(i in 1:length(TestConditions)){
      
      myCond <- TestConditions[!TestConditions %in% TestConditions[i]]
      myOtherCond <- TestConditions[i]
      
      tempo <- TableToPlot %>% spread(Treatment,Response) %>% gather_("Treatment","Response",gather_cols=as.vector(myCond))
      if(is.null(filterCluster)){
        filterTreatment2  <- unique(as.vector(myCond))
      }else{
        filterTreatment2  <-  filterCluster[filterCluster %in% unique(as.vector(myCond))]
      }
      tempo$Cluster <-   factor(tempo$Cluster,levels = filterTreatment2)
      dd[[i]] <-  tempo %>% 
        ggplot(aes_string(y=as.vector(myOtherCond),x="Response",colour="Treatment"))+geom_point(size=2,alpha=0.8)+theme_bw()+
        ylim(MaxAndMin[1],MaxAndMin[2])+xlim(MaxAndMin[1],MaxAndMin[2])
      if(length(unique(TableToPlot$Treatment)) > 1){
        dd[[i]] <-  dd[[i]] +facet_grid(~Treatment)
      }
    }  
  }
  if(!is.null(colourMap)){
    dd <- lapply(dd,function(x,colourMap)x+scale_colour_manual(values=colourMap),colourMap=colourMap)
  }
  require(gridExtra)
  MP <- marrangeGrob(dd,nrow=length(TestConditions),ncol=1)
  return(list(MP,dd))
}



expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

testClusters <- function(heatOut,clusterGroups=NULL,rawValues=FALSE){
  require(perm)
  
  
  # clusterGroups=list(
  #   list(
  #     Cluster=1,
  #     groups=list(A="Test",C="Test4",B="Test2")
  #   ),
  #   list(
  #     Cluster=2,
  #     groups=list(A="Test",B="Test2",c="Test4")
  #   )
  # )
  
  # heatOut <- swsw
  
  if(is.null(clusterGroups)){
    clusterGroups <- lapply(levels(swsw$CllusterInfo$CutMembers$Cluster),function(x){
      l=rownames(swsw$Values)
      names(l) <- l
      list(groups=as.list(l),
           Cluster=x)
      
    }
    )
    
  }
  # heatOut <- swsw
  # clusterGroups=list(list(Cluster=1,groups=list(A=c("Test","Test4"),B="Test2")))
  require(tidyverse)
  require(ggplot2)
  require(gridExtra)
  require(dunn.test)
  graphics.off()
  
  if(rawValues){
    values <- heatOut$OriginalValues %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }else{
    values <- heatOut$Values %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }
  memberShip <- heatOut$CllusterInfo$CutMembers %>%  as.data.frame %>% rownames_to_column(var = "Neuron")
  TableToPlot <- merge(memberShip,gather(values,key = "Neuron",value= "Response",-Treatment),by="Neuron")
  TableToPlot$Cluster <- factor(TableToPlot$Cluster)
  dugEphresh <- vector("list",length(clusterGroups))
  
  for(i in 1:length(clusterGroups)){
    tempEW <- clusterGroups[[i]]
    tempClusterTest <- TableToPlot[TableToPlot$Cluster %in% tempEW$Cluster,]
    valuesToTest <- vector("list",length(tempEW$groups))
    for(k in 1:length(valuesToTest)){
      valuesToTest[[k]] <- tempClusterTest[tempClusterTest$Treatment %in% tempEW$groups[[k]],"Response"]
    }
    names(valuesToTest) <-   names(tempEW$groups)
    temp <- expand.grid.unique(names(tempEW$groups),names(tempEW$groups)) %>% 
      as.data.frame %>%  
      transmute(GroupA=V1,GroupB=V2)
    res <- list()
    for(p in 1:nrow(temp)){
      
      GroupA_Name <- temp[p,] %>% pull(GroupA) %>% as.vector
      GroupA_Values <- valuesToTest[[GroupA_Name]]
      GroupA_mean <- mean(GroupA_Values)
      GroupA_median <- median(GroupA_Values)
      GroupA_sd <- sd(GroupA_Values)
      GroupA_Treatments=tempEW$groups[GroupA_Name] %>% unlist %>% paste(sep="",collapse=";")
      
      GroupB_Name <- temp[p,] %>% pull(GroupB) %>% as.vector
      GroupB_Values <- valuesToTest[[GroupB_Name]]
      GroupB_mean <- mean(GroupB_Values)
      GroupB_median <- median(GroupB_Values)
      GroupB_sd <- sd(GroupB_Values)
      GroupB_Treatments=tempEW$groups[GroupB_Name] %>% unlist %>% paste(sep="",collapse=";")
      
      wilcoxRes <- wilcox.test(GroupA_Values,GroupB_Values,paired = FALSE,alternative = "two.sided")
      ttestRes <- t.test(GroupA_Values,GroupB_Values,paired = FALSE,var.equal = FALSE,alternative = "two.sided")
      ttestRes_statistic <- ttestRes$statistic
      ttestRes_pvalue <- ttestRes$p.value
      wilcoxRes_statistic <- wilcoxRes$statistic
      wilcoxRes_pvalue <- wilcoxRes$p.value
      testPerm <- permTS(GroupA_Values,
                         GroupB_Values,method = "pclt",alternative = "two.sided")
      permutation_Pvalue <- testPerm$p.value
      
      res[[p]] <- data.frame( 
        Cluster=tempEW$Cluster %>% paste(sep="",collapse=";"),
        Comparison=paste0(GroupA_Name," - ",GroupB_Name),
        GroupA_Name,
        GroupA_Treatments,
        GroupB_Name,
        GroupB_Treatments,
        GroupA_mean,
        GroupA_median,
        GroupA_sd,
        GroupB_mean,
        GroupB_median,
        GroupB_sd,
        ttestRes_statistic,
        ttestRes_pvalue,
        wilcoxRes_statistic,
        wilcoxRes_pvalue,
        permutation_Pvalue)
    }
    dugEphreshTemp <- do.call(rbind,res)
    myFrame <- dunn.test(unname(unlist(valuesToTest)),rep(names(valuesToTest),lengths(valuesToTest)),table = FALSE)
    dugEphresh2Temp2 <- data.frame(Cluster=tempEW$Cluster %>% paste(sep="",collapse=";"),
                                   Comparison=myFrame$comparisons,
                                   Zscore=myFrame$Z,
                                   Pvalue=myFrame$P
    )
    dugEphresh[[i]] <- merge(dugEphreshTemp,dugEphresh2Temp2,all=TRUE,by="Comparison")
    
    
  }
  slikRic <- do.call(rbind,dugEphresh)
  return(slikRic)
  
}

testClusters_2 <- function(heatOut,nonLigands=c("mix","combo"),rawValues=FALSE,MuFor1SampleTest=0){
  require(perm)
  
  
  # clusterGroups=list(
  #   list(
  #     Cluster=1,
  #     groups=list(A="Test",C="Test4",B="Test2")
  #   ),
  #   list(
  #     Cluster=2,
  #     groups=list(A="Test",B="Test2",c="Test4")
  #   )
  # )
  
  # heatOut <- swsw
  
  clusterGroups <- lapply(levels(swsw$CllusterInfo$CutMembers$Cluster),function(x){
    l=rownames(swsw$Values)
    names(l) <- l
    list(groups=as.list(l),
         Cluster=x)
    
  }
  )
  
  
  # heatOut <- swsw
  # clusterGroups=list(list(Cluster=1,groups=list(A=c("Test","Test4"),B="Test2")))
  require(tidyverse)
  require(ggplot2)
  require(gridExtra)
  require(dunn.test)
  require(perm)
  graphics.off()
  
  if(rawValues){
    values <- heatOut$OriginalValues %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }else{
    values <- heatOut$Values %>% as.data.frame %>% rownames_to_column(var = "Treatment")
  }
  memberShip <- heatOut$CllusterInfo$CutMembers %>%  as.data.frame %>% rownames_to_column(var = "Neuron")
  TableToPlot <- merge(memberShip,gather(values,key = "Neuron",value= "Response",-Treatment),by="Neuron")
  TableToPlot$Cluster <- factor(TableToPlot$Cluster)
  dugEphresh <- vector("list",length(clusterGroups))
  
  for(i in 1:length(clusterGroups)){
    tempEW <- clusterGroups[[i]]
    tempClusterTest <- TableToPlot[TableToPlot$Cluster %in% tempEW$Cluster,]
    valuesToTest <- vector("list",length(tempEW$groups))
    for(k in 1:length(valuesToTest)){
      valuesToTest[[k]] <- tempClusterTest[tempClusterTest$Treatment %in% tempEW$groups[[k]],"Response"]
    }
    names(valuesToTest) <-   names(tempEW$groups)
    
    if(is.null(nonLigands)){
      nonLigands <- ""
    }
    
    valuesToTest <- valuesToTest[!names(valuesToTest) %in% nonLigands]
    
    res <- list()
    for(p in 1:length(valuesToTest)){
      
      GroupA_Name <- names(valuesToTest)[p]
      GroupA_Values <- valuesToTest[[GroupA_Name]]
      GroupA_mean <- mean(GroupA_Values)
      GroupA_median <- median(GroupA_Values)
      GroupA_sd <- sd(GroupA_Values)
      GroupA_Treatments=tempEW$groups[GroupA_Name] %>% unlist %>% paste(sep="",collapse=";")
      
      GroupB_Name <- paste(names(valuesToTest[!names(valuesToTest) %in% GroupA_Name]),collapse = ",")
      GroupB_Values <- valuesToTest[!names(valuesToTest) %in% GroupA_Name] %>% unlist %>% unname
      GroupB_mean <- mean(GroupB_Values)
      GroupB_median <- median(GroupB_Values)
      GroupB_sd <- sd(GroupB_Values)
      GroupB_Treatments=tempEW$groups[!tempEW$groups %in% GroupA_Name] %>% unlist %>% paste(sep="",collapse=";")
      
      wilcoxRes <- wilcox.test(GroupA_Values,GroupB_Values,paired = FALSE,alternative = "two.sided")
      wilcoxRes2 <- wilcox.test(GroupA_Values, mu = MuFor1SampleTest, alternative = "greater")
      ttestRes <- t.test(GroupA_Values,GroupB_Values,paired = FALSE,var.equal = FALSE,alternative = "two.sided")
      ttestRes_statistic <- ttestRes$statistic
      ttestRes_pvalue <- ttestRes$p.value
      wilcoxRes_statistic <- wilcoxRes$statistic
      wilcoxRes_pvalue <- wilcoxRes$p.value
      OneSample_WilcoxRes_statistic <- wilcoxRes2$statistic
      OneSample_WilcoxRes_pvalue <- wilcoxRes2$p.value
      
      testPerm <- permTS(GroupA_Values,
                         GroupB_Values,method = "pclt",alternative = "two.sided")
      permutation_Pvalue <- testPerm$p.value
      
      
      
      res[[p]] <- data.frame( 
        Cluster=tempEW$Cluster %>% paste(sep="",collapse=";"),
        Comparison=paste0(GroupA_Name," - ",GroupB_Name),
        GroupA_Name,
        GroupA_Treatments,
        GroupB_Name,
        GroupB_Treatments,
        GroupA_mean,
        GroupA_median,
        GroupA_sd,
        GroupB_mean,
        GroupB_median,
        GroupB_sd,
        ttestRes_statistic,
        ttestRes_pvalue,
        wilcoxRes_statistic,
        wilcoxRes_pvalue,
        OneSample_WilcoxRes_statistic,
        OneSample_WilcoxRes_pvalue,
        permutation_Pvalue)
    }
    
    dugEphresh[[i]] <-  do.call(rbind,res)
    
    
  }
  slikRic <- do.call(rbind,dugEphresh)
  return(slikRic)
  
}

###########

#########

