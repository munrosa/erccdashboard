loadERCCInfo <- function(expDat, erccmix = NULL, erccversion = NULL, userMixFile=NULL){
    # Get the ERCC Mix definition file provided by user and combine it with the 
    ## package ERCCDef file
    
    data(ERCC, envir = environment())
    
    stopifnot(!is.null(erccversion))
    
    if (erccversion == "ERCC1"){
      if(is.null(userMixFile)){
        MixDef <- ERCCMix1and2  
      }else{
        MixDef <- read.csv(userMixFile)
      }
      ERCCDef <- ERCCDef
    }
    if (erccversion == "ERCC2"){
      if(is.null(userMixFile)){
        MixDef <- ERCC2Mixes
      }else{
        MixDef <- read.csv(userMixFile)
      }
      ERCCDef <- ERCC2Def
    }
    
    
    names(MixDef)[1:2] <- c("Feature","Ratio") 
    
    # Sort by the feature column
    MixDefSort <- MixDef[do.call(order, MixDef[c(1)]), ]
    MixDef <- MixDefSort
    
    # Check that the MixDef has the same ERCCs as the ERCCDef data frame
    ERCCmatch <- ERCCDef[match(MixDef$Feature,ERCCDef$Feature),]
    
    # Combine the Mix definition and ERCCdef data frames
    idCols <- merge(ERCCmatch,MixDef, by = "Feature")
    idCols$Ratio <- as.factor(idCols$Ratio)
    #FCcode <- data.frame(Ratio=c("a","b","c","d"),#Default for Ambion pool
    #                         FC =  c(4,1,.667,.5))
    #legendLabels = c("4:1","1:1","1:1.5","1:2")#Default for Ambion pool
    
    if(erccmix == "Single"){
        # Fix the names and column identity, assigning the mix 1 concentration 
        # values to both columns
        idCols[c(6)] <- idCols[c(5)]   
        levels(idCols$Ratio) <- c("a1:1","b1:1","c1:1","d1:1")
    }
    
    names(idCols)[5:6] <- c("Conc1", "Conc2")
    
    FoldChanges <- round(idCols$Conc1/idCols$Conc2,digits=3)
    
    FC = rep(1:nlevels(idCols$Ratio), 0)
    for (i in 1:nlevels(idCols$Ratio)){
        FC[i] <- unique(FoldChanges[which(idCols$Ratio == 
                                              levels(idCols$Ratio)[i])])
    }
    FCcode = data.frame(Ratio = levels(idCols$Ratio), FC = FC) 
    
    legendLabels = as.character(FCcode$Ratio)
    erccInfo = list(idColsSRM = idCols, MixDef = MixDef, FCcode = FCcode, 
                    legendLabels = legendLabels)
    expDat$erccInfo <- erccInfo
    if(!exists("erccmix")){
        stop(paste("Please define the erccmix variable as either:",
                   "\"RatioPair\", or \"Single\""))
    }
    return(expDat)
}