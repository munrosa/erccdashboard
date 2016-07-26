#' Run default erccdashboard analysis of ERCC control ratio mixtures
#'
#' @param datType       type is "count" (RNA-Seq) or "array" (microarray),
#'                      "count" is unnormalized integer count data (normalized
#'                      RNA-Seq data will be accepted in an updated version of
#'                      the package), "array" can be normalized or unnormalized 
#'                      fluorescent intensities from a microarray experiment.
#' @param isNorm        default is FALSE, if FALSE then the unnormalized
#'                      input data will be
#'                      normalized in erccdashboard analysis. If TRUE then
#'                      it is expected that the data is already normalized
#' @param exTable       data frame, the first column contains names of 
#'                      genes or transcripts (Feature) and the remaining columns
#'                      are expression measures for sample replicates spiked 
#'                      with ERCC controls
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
#' @param ratioLim      Limits for ratio axis on MA plot, default is c(-4,4)
#' @param signalLim      Limits for ratio axis on MA plot, default is c(-14,14)
#' @param userMixFile   optional filename input, default is NULL, if ERCC 
#'                      control ratio mixtures other than the Ambion product
#'                      were used then a userMixFile can be used for the analysis
#'                                         
#' 
#' @export
#' @examples
#' 
#' \donttest{
#' data(SEQC.Example)
#'      
#' exDat <- runDashboard(datType = "count",isNorm = FALSE,
#'                      exTable = MET.CTL.countDat, 
#'                      filenameRoot = "COH.ILM",
#'                      sample1Name = "MET", sample2Name = "CTL", 
#'                      erccmix = "RatioPair", erccdilution = 1/100, 
#'                      spikeVol = 1, totalRNAmass = 0.500,choseFDR = 0.1)
#'                  
#' summary(exDat)
#' }
#' 

runDashboard <- function(datType=NULL, isNorm = FALSE, 
                         exTable=NULL, repNormFactor=NULL,
                         filenameRoot = NULL,
                         sample1Name = NULL,sample2Name = NULL, 
                         erccmix = "RatioPair", erccdilution = 1,
                         spikeVol = 1, totalRNAmass = 1,choseFDR = 0.05,
                         ratioLim=c(-4,4),signalLim=c(-14,14),
                         userMixFile=NULL){
    
    # Initialize exDat structure
    # Required for all subsequent functions
    exDat <- initDat(datType=datType, isNorm = isNorm, exTable=exTable, 
                     repNormFactor=repNormFactor, filenameRoot=filenameRoot,
                     sample1Name=sample1Name, sample2Name=sample2Name, 
                     erccmix=erccmix, erccdilution=erccdilution, 
                     spikeVol=spikeVol, totalRNAmass=totalRNAmass,
                     choseFDR=choseFDR,ratioLim = ratioLim,
                     signalLim = signalLim,userMixFile=userMixFile)
    
    # Estimate the difference in mRNA fraction of total RNA for the two samples
    # Required for all subsequent functions
    exDat <- est_r_m(exDat)
    
    # Evaluate the dynamic range of the experiment (Signal-Abundance plot)
    # Not required for subsequent functions
    exDat <- dynRangePlot(exDat)
    
    # Test for differential expression between samples
    # Required for all subsequent functions
    exDat <- geneExprTest(exDat)
    
    # Generate ROC curves and AUC statistics
    # Not Required for subsequent functions
    exDat <- erccROC(exDat)
    
    # Estimate LODR for ERCC controls
    # Required for subsequent functions
    exDat <- estLODR(exDat,kind = "ERCC", prob=0.9)
    
    ## Estimate LODR using Simulated data from endogenous transcripts
    ## Not required for subsequent functions
    #exDat = estLODR(exDat,kind = "Sim", prob=0.9)
    
    # Generate MA plot (Ratio vs. Average Signal) with ERCC controls below LODR 
    ## annotated also flags possible False Negatives on DE gene list based on 
    ## LODR threshold from DE gene list
    # Not required for subsequent functions
    exDat <- annotLODR(exDat)
    
    
    ### Saving plots and results
    # Convenience function to save 4 main figures to file
    saveERCCPlots(exDat,saveas = "pdf")
    
    
    # Save exDat to a RData file for later use
    cat("\nSaving exDat list to .RData file...")
    nam <- paste(exDat$sampleInfo$filenameRoot, "exDat",sep = ".")
    assign(nam,exDat)
    
    to.save <- ls()
    saveName <- paste0(exDat$sampleInfo$filenameRoot,".RData")
    save(list = to.save[grepl(pattern = nam,
                              x=to.save)],file = saveName)
    
    # End analysis and return exDat to global environ. / workspace
    cat("\nAnalysis completed.")
    return(exDat)
    dev.off()
    
}