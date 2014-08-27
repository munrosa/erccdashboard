prepERCCDat <- function(expDat){
    sampleInfo <- expDat$sampleInfo
    expressDat <- expDat$Transcripts
    idCols <- expDat$idCols
    # get just ERCC data in the expression data frame
    expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]
    
    # Length normalize the Expected ERCC concentrations
    kNorm = sampleInfo$kNorm
    if (kNorm %in% c("N","n","no","No")){
        lengthFactor = (idCols$Length) #/(1000)
        ERCCxlabelIndiv <- expression(paste("Log2 ERCC Spike ", 
                                            "Amount (attomol nt ", mu,
                                           "g"^"-1"," total RNA)", sep = ""))
        
        ERCCxlabelAve <- expression(paste("Log2 Average ERCC Spike ",
                                          "Amount (attomol nt ", mu,
                                          "g"^"-1"," total RNA)", sep = ""))
    }
    if(kNorm %in% c("Y","y","yes","Yes")){
        lengthFactor = 1
        ERCCxlabelIndiv = expression(paste("Log2 ERCC Spike ",
                                           "Amount (attomol nt ", mu,
                                           "g"^"-1"," total RNA)", sep = ""))
        ERCCxlabelAve = expression(paste("Log2 Average ERCC Spike ",
                                         "Amount (attomol nt ", mu,
                                         "g"^"-1"," total RNA)", sep = ""))
        
    }
    
    
    
    #If length normalization of the expected concentrations is desired (default)
    idCols$Conc1 = (idCols$Conc1*lengthFactor)
    idCols$Conc2 = (idCols$Conc2*lengthFactor)  
    
    ### Calculate the per ERCC amount spiked attomoles of nt / ug total RNA
    ### ERCCdilution = 1; spikeVol = 50; totalRNAmass = 2.5*10^(3)
    spikeFraction <- (sampleInfo$erccdilution*sampleInfo$spikeVol)/(sampleInfo$totalRNAmass)
    idCols$Conc1 <- idCols$Conc1*spikeFraction
    idCols$Conc2 <- idCols$Conc2*spikeFraction
    
    
    expDat$idCols <- idCols
    expDat$plotInfo$ERCCxlabelIndiv <- ERCCxlabelIndiv
    expDat$plotInfo$ERCCxlabelAve <- ERCCxlabelAve
    expDat$spikeFraction <- spikeFraction
    return(expDat)
}