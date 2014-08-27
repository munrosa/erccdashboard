dashboardFile <- function(exDat, filenameRoot){
    
    ## Assign local variables
    sampleInfo <- exDat$sampleInfo
    
    ## Create filename for results files
    filenameUse <- paste(filenameRoot,sampleInfo$sample1Name,
                        sampleInfo$sample2Name,sep=".")  
    
    cat(paste("Filename root is:", filenameUse, "\n"))
    
    exDat$sampleInfo$filenameRoot <- filenameUse
    
    return(exDat)
}