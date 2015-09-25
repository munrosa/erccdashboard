#' Initialize the exDat list
#'
#' @param datType       type is "count" or "array", unnormalized data is  
#'                      expected (normalized data may be accepted in future
#'                      version of the package). Default is "count" (integer 
#'                      count data),"array" is unnormalized fluorescent 
#'                      intensities from microarray
#'                      fluorescent intensities (not log transformed or 
#'                      normalized)
#' @param isNorm        default is FALSE, if FALSE then the unnormalized
#'                      input data will be
#'                      normalized in erccdashboard analysis. If TRUE then
#'                      it is expected that the data is already normalized
#' @param exTable      data frame, the first column contains names of 
#'                      genes or transcripts (Feature) and the remaining columns
#'                      are counts for sample replicates spiked with ERCC 
#'                      controls
#' @param repNormFactor optional vector of normalization factors for each 
#'                      replicate, default value is NULL and 75th percentile
#'                      normalization will be applied to replicates
#' @param filenameRoot  string root name for output files
#' @param sample1Name   string name for sample 1 in the gene expression 
#'                      experiment
#' @param sample2Name   string name for sample 2 in the gene expression
#'                      experiment
#' @param erccmix     Name of ERCC mixture design, "RatioPair" is 
#'                      default, the other option is "Single"
#' @param erccdilution  unitless dilution factor used in dilution of the Ambion 
#'                      ERCC spike-in mixture solutions 
#' @param spikeVol      volume in microliters of diluted ERCC mix spiked into
#'                      the total RNA samples
#' @param totalRNAmass  mass in micrograms of total RNA spiked with diluted ERCC
#'                      mixtures 
#' @param choseFDR      False Discovery Rate for differential expression testing
#'                      , default is 0.05
#' @param ratioLim      Limits for ratio axis on MA plot, default is c(-4,4)
#' @param signalLim     Limits for signal axis on dynamic range plot, default 
#'                      is c(-14,14)
#' @param userMixFile   optional filename input, default is NULL, if ERCC 
#'                      control ratio mixtures other than the Ambion product
#'                      were used then a userMixFile can be used for the 
#'                      analysis
#' @examples
#' 
#' data(SEQC.Example)
#' 
#' exDat <- initDat(datType="count", isNorm = FALSE, exTable=MET.CTL.countDat, 
#'                  filenameRoot = "testRun",sample1Name = "MET",
#'                  sample2Name = "CTL", erccmix = "RatioPair", 
#'                  erccdilution = 1/100, spikeVol = 1, totalRNAmass = 0.500,
#'                  choseFDR = 0.1)
#' summary(exDat)                      
#'                                         
#' 
#' @export


initDat <- function(datType=NULL, isNorm=FALSE, exTable=NULL, 
                    repNormFactor=NULL, filenameRoot=NULL,
                    sample1Name=NULL, sample2Name=NULL, 
                    erccmix="RatioPair", erccdilution=1,
                    spikeVol=1, totalRNAmass=1,choseFDR=0.05,
                    ratioLim=c(-4,4), signalLim=c(-14,14), 
                    userMixFile=NULL){
    cat("\nInitializing the exDat list structure...\n")
    
    if (missing(datType))
        stop("Need to specify datType.")
    if (missing(isNorm))
        stop("Need to specify isNorm.")
    if (missing(exTable))
        stop("Need to specify exTable.")
    if(missing(filenameRoot))
        stop("Need to specify filenameRoot.")
    if(missing(sample1Name))
        stop("Need to specify sample1Name.")
    if(missing(sample2Name))
        stop("Need to specify sample2Name.")
    
    myYLimMA <- ratioLim
    myXLimMA <- signalLim
    xlimEffects <- c(-15,15)
    
    myYLim <- myXLimMA
    myXLim <- NULL
    
    exDat<-NULL
    
    cat(paste("choseFDR =",choseFDR,"\n"))
    
    if(missing(userMixFile)){
        userMixFile <- NULL
    }
    
    if(is.null(repNormFactor)){
        #repNormFactor <- NULL
        cat("repNormFactor is NULL \n")
    }
    if(isNorm == TRUE){
        cat(paste("\nisNorm is TRUE, input data will be considered",
                  "to be normalized\n"))
        getKNorm<- function(){
            cat(paste("\nIs the expression data length normalized",
                      "(e.g. FPKM or RPKM)?\n"))
            readline("Enter Y or N: ")
        }
        kNorm <- as.character(getKNorm())  
    }else{
        kNorm = "N"
    }
    
    ##############################
    
    sampleInfo = list(sample1Name = sample1Name,
                      sample2Name = sample2Name, choseFDR = choseFDR,
                      erccdilution = erccdilution, erccmix = erccmix,
                      spikeVol = spikeVol, totalRNAmass = totalRNAmass,
                      isNorm = isNorm, kNorm = kNorm, datType = datType)
    
    plotInfo = list(myXLimMA = myXLimMA, myYLimMA = myYLimMA,
                    myXLim = myXLim, myYLim = myYLim, xlimEffects = xlimEffects)
    
    
    exDat <- list(sampleInfo = sampleInfo,plotInfo = plotInfo)
    
    exDat <- dashboardFile(exDat=exDat,filenameRoot=filenameRoot)  
    
    ############################################################################
    # Run loadERCCInfo function to obtain ERCC information
    exDat <- loadERCCInfo(exDat, erccmix, userMixFile)
    
    ############################################################################ 
    # Add experimental data (exTable) to exDat structure
    exDat <- loadExpMeas(exDat, exTable, repNormFactor)
    
    ############################################################################
    # normalize the data
    #if(isNorm == FALSE){
    exDat <- normalizeDat(exDat)  
    #}
    
    ############################################################################
    # length normalize the ERCC concentrations
    
    exDat <- prepERCCDat(exDat)
        
    exDat <- plotAdjust(exDat)
    
    # check for sample 1 names to be in the first set of exTable sample columns
    sample1Name <- exDat$sampleInfo$sample1Name
    
    col2Sample <- gsub(pattern = "_1",replacement = "",
                       x = names(exDat$Transcripts)[2])
    
    if(sample1Name != col2Sample){
        stop(paste("\nThe first set of columns in exTable are for sample",
                   col2Sample,"and \nmust match sample1Name =",sample1Name,
                   "\n","Please check sample and exTable names\n"))
    }
    
    # check for missing 1:1 ERCCs 
    if(!(any(exDat$erccInfo$FCcode$FC == 1))){
        stop("1:1 controls are missing, check userMixFile")    
    }
    
    
    return(exDat)
    
}