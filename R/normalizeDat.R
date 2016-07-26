normalizeDat <- function(expDat){
    
    sampleInfo <- expDat$sampleInfo
    erccInfo <- expDat$erccInfo
    expressDat = expDat$Transcripts
    datType <- sampleInfo$datType
    repNormFactor <- sampleInfo$repNormFactor
    
    normVec = TRUE
    
    if (is.null(sampleInfo$repNormFactor)){
        normVec = FALSE
    }
    
    # Normalize the data
    if (sampleInfo$isNorm == FALSE){
      #  if(sampleInfo$datType == "array"){
            if(normVec == FALSE){
                cat(paste("\nrepNormFactor is NULL,\n",
                          "Using Default Upper Quartile Normalization Method",
                          " - 75th percentile\n"))
                
                normFactor = apply(expressDat[-c(1)],MARGIN=2,
                                   FUN=quantile,probs=0.75,na.rm=TRUE)
                
                libAdjust = sweep(expressDat[-c(1)],2,normFactor,"/")
                expressDat = cbind(expressDat[c(1)], libAdjust)
            }
            if(normVec == TRUE){
                cat(paste("\nUsing normalization factors provided by user",
                          "in repNormFactor\n"))

                normFactor = repNormFactor
                
                datCols = expressDat[-c(1)]
                
                #normalize the data  
                libAdjust = sweep(datCols, 2, normFactor,"/")
                expressDat = cbind(expressDat[c(1)], libAdjust)
            }
       # }
       # if(sampleInfo$datType == "count"){
#             if(normVec == FALSE){
#                 cat(paste("\nrepNormFactor is NULL,\n",
#                           "Using Default Upper Quartile Normalization Method",
#                           " - 75th percentile\n"))
#                 normFactor = apply(expressDat[-c(1)],MARGIN=2,
#                                    FUN=quantile,probs=0.75)
#                 datCols = expressDat[-c(1)]
#                 #normFactor = normFactor#/(10^6) #per million mapped reads
#                 #normalize the data  
#                 libAdjust = sweep(datCols, 2, normFactor,"/")
#                 expressDat = cbind(expressDat[c(1)], libAdjust)
#             }
#             if(normVec == TRUE){
#                 cat(paste("\nUsing normalization factors provided",
#                           "in repNormFactor by user\n"))
#                 normFactor = repNormFactor
#                 #print(normFactor)
#                 datCols = expressDat[-c(1)]
#               #  normFactor = normFactor/(10^6) #per million mapped reads
#                  
#                 #Library size normalize the data  
#                 libAdjust = sweep(datCols, 2, normFactor,"/")
#                 expressDat = cbind(expressDat[c(1)], libAdjust)
#             }
#         }
        
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