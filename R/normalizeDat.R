normalizeDat <- function(expDat){
    
    sampleInfo <- expDat$sampleInfo
    erccInfo <- expDat$erccInfo
    expressDat = expDat$Transcripts
    datType <- sampleInfo$datType
    repNormFactor <- sampleInfo$repNormFactor
    #print(sampleInfo$repNormFactor)
    normVec = TRUE
    if(sampleInfo$datType == "count"){
        if (is.null(sampleInfo$repNormFactor)){
            normVec = FALSE
        }
        # print(normVec)
    }
    if(sampleInfo$datType == "array"){
        normVec = FALSE
    }
    
    # Library size normalize the data
    if (sampleInfo$isNorm == FALSE){
        if(sampleInfo$datType == "array"){
            #cat("\nUsing median intensity for ERCC 1:1 controls\n")
            #cat("\nUse median intensity from each array for normalization\n")
            #cat("\nUsing total intensity to normalize each array\n")
            cat(paste("\nUsing 75th percentile (upper quartile) intensity",
                      "to normalize each array\n"))
            ### Subset the 1:1 ercc controls and calculate the median intensity 
            ## for each column, build a normFactor vector  
            TranscriptsAll = expressDat
            normFactor = apply(expressDat[-c(1)],MARGIN=2,
                               FUN=quantile,probs=0.75)
            #normFactor  = colSums(TranscriptsAll[-c(1)])
            #print(normFactor)
            
            libAdjust = sweep(expressDat[-c(1)],2,normFactor,"/")
            expressDat = cbind(expressDat[c(1)], libAdjust)
        }
        if(sampleInfo$datType == "count"){
            if (normVec == FALSE){
                cat(paste("\nrepNormFactor is NULL,\n",
                          "Using Default Upper Quartile Normalization Method",
                          " - 75th percentile\n"))
                TranscriptsAll = expressDat
                normFactor = apply(expressDat[-c(1)],MARGIN=2,
                                   FUN=quantile,probs=0.75)
                datCols = expressDat[-c(1)]
                normFactor = normFactor#/(10^6) #per million mapped reads
                #Library size normalize the data  
                libAdjust = sweep(datCols, 2, normFactor,"/")
                expressDat = cbind(expressDat[c(1)], libAdjust)
            }
            
            if (normVec == TRUE){
                cat(paste("\nUsing read depth normalization factors provided",
                          "in repNormFactor\n"))
                TranscriptsAll = expressDat
                normFactor = repNormFactor
                #print(normFactor)
                datCols = expressDat[-c(1)]
                normFactor = normFactor/(10^6) #per million mapped reads
                
                #Library size normalize the data  
                libAdjust = sweep(datCols, 2, normFactor,"/")
                expressDat = cbind(expressDat[c(1)], libAdjust)
            }
        }
        
        cat("\nnormVec:\n")
        cat(normFactor)
        mnLibeFactor = mean(normFactor)
    }
    else{
        normFactor = NULL
        mnLibeFactor = NULL
    }
    # get just ERCC data in the expression data frame
    expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]
    
    
    expDat$normERCCDat <- expressDat
    expDat$normFactor <- normFactor
    expDat$mnLibeFactor <- mnLibeFactor
    
    
    return(expDat)
    
}