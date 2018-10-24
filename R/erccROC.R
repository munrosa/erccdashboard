#' Produce Receiver Operator Characteristic (ROC) Curves and AUC statistics
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @examples
#' 
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
#' 
#' exDat <- erccROC(exDat)
#' 
#' exDat$Figures$rocPlot
#'
#' @export
#' 
erccROC <- function(exDat){
    
    filenameRoot <- exDat$sampleInfo$filenameRoot
    folds <- exDat$erccInfo$FCcode
    legendLabels <- exDat$erccInfo$legendLabels
    idCols <- exDat$erccInfo$idColsSRM
    
    idCols <- idCols[-which(is.na(idCols$Ratio)),]
    
    legendLabelsDiff <- legendLabels[-(which(folds$FC == 1))]
    
    plotInfo <- exDat$plotInfo
    colScale <- plotInfo$colScale
    fillScale <- plotInfo$fillScale
    
    
    # Read in the p.values from the file
    #if (is.null(pValDat)){
    erccPval <- file.exists(paste(filenameRoot,"ERCC","Pvals.csv"))
    allPval <- file.exists(paste0(filenameRoot, ".All.Pvals.csv"))
    if((erccPval == TRUE)&(allPval==TRUE)){
        pValDat <- read.csv(file=paste(filenameRoot,"ERCC","Pvals.csv"),
                           header=TRUE)  
    }
    if(allPval == TRUE){
        pValDat <- read.csv(file=paste0(filenameRoot, ".All.Pvals.csv"),
                            header=TRUE)
        pValDat <- pValDat[grep("ERCC-00",pValDat$Feature),]
    }else{
        cat(paste0("\n",filenameRoot," ERCC Pvals.csv file is missing."))
        cat("\nExiting ROC Curve analysis...\n")
        return(exDat)
    }
    
    cat("\nGenerating ROC curve and AUC statistics...\n")
    names(pValDat)[1]= "Feature"
    pValDat <- pValDat[-2]
    
    pValDat$Feature <- as.character(pValDat$Feature)
    
    # now build the ROCR prediction objects
    # Format of FCcode <- data.frame(Ratio=c("a","b","c","d"),
    #                                   FC=c(4, 1,.667,.5));
    # First match the fold changes in pValDat
    folds <- folds[folds$FC %in% pValDat$Fold,]
    FCcodeC <- folds[-c(which(folds$FC == 1)),]
    pool.pred <- NULL
    FPR <- NULL
    TPR <- NULL
    FoldChange <- NULL
    ROCdat <- NULL
    AUCdat <-NULL
    for (i in 1:nrow(FCcodeC)){
        
        pool.predMeas <- prediction(1- pValDat$Pval[ (pValDat$Fold == 
                                                         FCcodeC$FC[i]) |
                                                        (pValDat$Fold == 1 )], 
                                   pValDat$Fold[(pValDat$Fold == 
                                                     FCcodeC$FC[i]) | 
                                                    (pValDat$Fold == 1)],
                                   label.ordering=c(1, FCcodeC$FC[i]))
        
        pool.perf <- performance(pool.predMeas, "tpr","fpr")
        pool.auc <- performance(pool.predMeas, "auc")
    
        # now build the three vectors for plotting - TPR, FPR, and FoldChange
        AUC <- unlist(pool.auc@y.values)
        
        AUCdatnew <- data.frame(Ratio=legendLabelsDiff[i],
                               AUC=round(AUC, digits=3), 
                               Detected=
                                   length(pValDat$Fold[(pValDat$Fold == 
                                                            FCcodeC$FC[i])]), 
                               Spiked=
                                   length(idCols$Ratio[idCols$Ratio == 
                                                           FCcodeC$Ratio[i]]))
        AUCdat <- rbind(AUCdat, AUCdatnew)
        FPR <- c( unlist(pool.perf@x.values)) 
        TPR <- c( unlist(pool.perf@y.values))
        #print(TPR)
        #print(FPR)
        Ratio <- c(rep(as.character(FCcodeC$Ratio[i]), 
                      length(unlist(pool.perf@y.values))))
        ROCdatnew <- data.frame(FPR=FPR, TPR=TPR, Ratio=Ratio)
        ROCdat <- rbind(ROCdat, ROCdatnew)
    }
    
    AUCAnnot <- AUCdat
    cat("\nArea Under the Curve (AUC) Results:\n")
    print(AUCAnnot, quote = FALSE, row.names = FALSE)
    
    AUCdat$xval <- 0.7
    AUCdat$yval <- seq(to=0.25, from=0.1, length.out=nrow(FCcodeC))
    
    ROCplot <- ggplot(data=ROCdat, aes(x=FPR, y=TPR)) + 
        geom_path(size=2, aes(colour=Ratio), alpha=0.7) + 
        geom_point(size=5, aes(colour=Ratio), alpha=0.7) + 
        colScale + geom_abline(intercept=0, slope=1, linetype=2) +
        theme_bw() + 
        annotation_custom(grob=tableGrob(AUCAnnot, rows=NULL),
                          xmin=0.375, xmax=1.0, ymin=0, ymax=0.25) +
        
#         annotation_custom(grob=
#                               tableGrob(AUCAnnot, show.rownames=FALSE, 
#                                         equal.width=TRUE, 
#                                         equal.height=TRUE,
#                                         gpar.corefill=gpar(fill="grey85",
#                                                              col="white"), 
#                                         gpar.rowfill=gpar(fill="grey80",
#                                                             col="white"),
#                                         gpar.colfill=gpar(fill="grey80",
#                                                             col="white")),
#                           xmin=0.375, xmax=1.0, ymin=0, ymax=0.25) +
        theme(legend.position=c(0.75, 0.5))
    
    exDat$Figures$rocPlot <- ROCplot
    exDat$Results$AUCdat <- AUCAnnot
    return(exDat)
    
}