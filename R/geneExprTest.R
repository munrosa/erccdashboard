#' Prepare differential expression testing results for spike-in analysis
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @details
#' This function wraps the edgeR differential expression testing package for
#' datType = "count" or uses the limma package for differential expression 
#' testing if datType = "array". Alternatively, for count data only, if
#' correctly formatted DE test results are provided,
#' then geneExprTest will bypass DE testing (with reduced runtime).
#' 
#' @examples
#' 
#' data(SEQC.Example)
#' 
#' exDat <- initDat(datType="array", isNorm=FALSE, 
#'                  exTable=UHRR.HBRR.arrayDat,
#'                  filenameRoot="testRun", sample1Name="UHRR",
#'                  sample2Name="HBRR", erccmix="RatioPair", 
#'                  erccdilution = 1, spikeVol = 50, 
#'                  totalRNAmass = 2.5*10^(3), choseFDR=0.01)
#'                  
#' exDat <- est_r_m(exDat)
#'                   
#' exDat <- dynRangePlot(exDat)
#' 
#' exDat <- geneExprTest(exDat)
#' 
#' @export

geneExprTest <- function(exDat){
    datType <- exDat$sampleInfo$datType
    isNorm <- exDat$sampleInfo$isNorm
    choseFDR <- exDat$sampleInfo$choseFDR
    cnt <- exDat$Transcripts
    designMat <- exDat$designMat
    sampleInfo <- exDat$sampleInfo
    info <- designMat
    
    allpvalFile <- paste(sampleInfo$filenameRoot,"All.Pvals.csv",sep=".")
    pvalERCC <- paste(sampleInfo$filenameRoot, "ERCC Pvals.csv",sep=" ")
    
    if(datType == "array"){
        if(is.null(choseFDR)){
            getPThresh<- function(){
                cat("\nFDR is NULL, to continue with LODR estimation\n")
                readline("Enter the threshold P-value: ")
            }
            exDat$Results$p.thresh <- as.numeric(getPThresh())  
        }
        
        if (file.exists(allpvalFile) == TRUE){
            deRes <- read.csv(allpvalFile)
            #if (!("qvals" %in% names(deRes))){
            deRes$qvals <- qvalue(deRes$Pval)$qvalues
            # print(summary(deRes$qvals))
            #}
            if(any(deRes$qvals<choseFDR)){
                p.thresh<-max(deRes$Pval[deRes$qvals<choseFDR])
            }
            cat(paste("\n Will use existing differential expression test",
                      "results for analysis.\n",
                      "Delete", allpvalFile, "if you want to repeat",
                      "differential \nexpression testing\n"))
        }else{
            exDat <- testDEArray(exDat)
            cat("\nFinished DE testing\n")
            p.thresh <- exDat$Results$p.thresh  
        }
        
        
    }
    if(datType == "count"){
        
        # set initial p.thresh
        p.thresh<-.1
        # First 3 columns of allpvalFile must contain Feature, Pval, and qvals
        # Decide to reuse results or run testDE
        
        if (file.exists(allpvalFile) == TRUE){
            deRes <- read.csv(allpvalFile)
            #if (("qvals" %in% names(deRes)) == FALSE){
            deRes$qvals <- qvalue(deRes$Pval)$qvalues
            #}
            if(any(deRes$qvals<choseFDR)){
                p.thresh<-max(deRes$Pval[deRes$qvals<choseFDR])
            }
            cat(paste("\n Found differential expression test results, will use",
                      " \nexisting P-values and Q-values for analysis.\n",
                      "Delete", allpvalFile, "if you want to repeat",
                      "differential \n",
                      "expression testing or view dispersion plots\n"))
        }else{
            if (isNorm == TRUE){
                cat(paste0("\nedgeR DE Testing for RNA-Seq requires count",
                           " (integer) data.\n",
                           "To estimate AUC and LODR for normalized RNA-Seq",
                           "data\nthe file '",sampleInfo$filenameRoot,
                           ".All.Pvals.csv' is required with columns for\n",
                           "'Feature','MnSignal','Pval', and 'Fold'\n"))
                return(exDat)
            }else{
                cat("\nStarting differential expression tests\n")
                exDat <- suppressWarnings(testDECount(sampleInfo, exDat, 
                                                      cnt=cnt, 
                                                      info=info))  
                deRes <- read.csv(allpvalFile)
                deRes$qvals <- qvalue(deRes$Pval)$qvalues
                if(any(deRes$qvals<choseFDR)){
                    p.thresh<-max(deRes$Pval[deRes$qvals<choseFDR])
                } 
            }
        }    
        
        
    }
    
    if(is.null(exDat$Figures$dispPlot)){
        cat(paste("\nDE testing results supplied without companion dispersion\n",
                  "plot. Dispersion plot is unavailable to print.\n"))
    }  
    cat("\nThreshold P-value\n")
    cat(p.thresh,"\n")
    
    if (p.thresh > .1){
        cat(paste("Threshold P-value is high for the chosen FDR of ", 
                  as.character(choseFDR)))
        cat(paste("\nThe sample comparison indicates a large amount of \n",
                  "differential expression in the measured transcript \n",
                  "populations\n"))
    }
    exDat$Results$p.thresh <- p.thresh  
    
    
    return(exDat)
    
}
