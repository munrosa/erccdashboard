#' Annotate signal-abundance and ratio-abundance plots with LODR
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @examples
#' data(SEQC.Example)
#'  
#' exDat <- initDat(datType="array", isNorm=FALSE, 
#'                  exTable=UHRR.HBRR.arrayDat,
#'                  filenameRoot="testRun", sample1Name="UHRR",
#'                  sample2Name="HBRR", erccmix="RatioPair",
#'                  erccversion = "ERCC1", 
#'                  erccdilution = 1, spikeVol = 50, 
#'                  totalRNAmass = 2.5*10^(3), choseFDR=0.01)
#'                  
#' exDat <- est_r_m(exDat)
#'                   
#' exDat <- dynRangePlot(exDat)
#' 
#' exDat <- geneExprTest(exDat)
#' \donttest{
#' exDat <- estLODR(exDat, kind="ERCC", prob=0.9)
#' 
#' exDat <- annotLODR(exDat)
#' 
#' exDat$Figures$maPlot
#' }
#' @export
#' 
annotLODR <- function(exDat){
    
    ## Assign LODR results to object in exDat
    LODR.annot.ERCC <- printLODRres(exDat)
    
    exDat$Results$LODR.annot.ERCC <- LODR.annot.ERCC
    
    ## Produce MA plots
    exDat <- maSignal(exDat)
    
    return(exDat)
}