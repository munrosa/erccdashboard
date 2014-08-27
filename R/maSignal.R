#' Generate MA plots with or without annotation using LODR estimates 
#'
#' @param exDat      list, contains input data and stores analysis results
#' @param alphaPoint  numeric value, for alpha (transparency) for plotted points,
#'                    range is 0 - 1
#' @param r_mAdjust   default is TRUE, if FALSE then the r_m estimate will not
#'                    used to offset dashed lines for empirical ratios on figure
#' @param replicate   default is TRUE, if FALSE then error bars will not be
#'                    produced
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
#' # generate MA plot without LODR annotation
#' exDat <- maSignal(exDat)
#' 
#' exDat$Figures$maPlot
#' \donttest{
#' exDat <- estLODR(exDat, kind = "ERCC", prob = 0.9)
#'
#' # Include LODR annotation
#' exDat <- annotLODR(exDat)
#' 
#' exDat$Figures$maPlot
#' }
#' @export

# Plots with target ratios and R adjusted ratios, MA plots and Ratio Summaries
maSignal <-function(exDat, alphaPoint=0.8, r_mAdjust=TRUE, replicate=TRUE){
    
    #ReplicateName = "Rep"
    # Melt the data
    #exDat <- meltexDat(exDat, cnt = exDat$Transcripts, 
    #                    designMat = exDat$designMat)
    #countPair <- exDat$expressDat_l
    sampleInfo <- exDat$sampleInfo
    erccInfo <- exDat$erccInfo
    plotInfo <- exDat$plotInfo
    
    cnt <- exDat$Transcripts
    designMat <- exDat$designMat
    sampleInfo <- exDat$sampleInfo
    normFactor <- exDat$normFactor
    datNames <- colnames(designMat)[-1]
    sample1 <- exDat$sampleNames[1]
    sample2 <- exDat$sampleNames[2]
    
    if (is.null(exDat$idColsAdj)){
        idColsAdj <- exDat$idCols
    }else{
        idColsAdj <- exDat$idColsAdj
    }
    
    datCols = cnt[-c(1)]
    if (!is.null(normFactor)){
        libAdjust = sweep(datCols, 2, normFactor,"/")
        sampleLibeDataNorm = cbind(cnt[c(1)],libAdjust)
        myDataERCC = sampleLibeDataNorm
        #expressDat = myDataERCC[-c(1)]
    }else{
        myDataERCC = cnt
    }
    
    sampleNameList = c(sample1,sample2)
    
    dat = cbind(Feature = myDataERCC[c(1)], Ratio = "Endo",
                myDataERCC[-c(1)])
    
    dat$Feature <- as.factor(as.character(
        dat$Feature))
    dat$Ratio <- as.character(dat$Ratio)
    #idx1 <- suppressWarnings(which(as.character(dat$Feature) == 
    #                                 as.character(exDat$idColsAdj$Feature)))
    idxERCC <- grep("ERCC-",dat$Feature)
    for (i in 1:length(idxERCC)){
        dat$Ratio[i] <- as.character(idColsAdj$Ratio)[match(dat$Feature[i],
                                                            idColsAdj$Feature)]
    }
    
    #idx1 <- match( exDat$idColsAdj$Feature, dat$Feature,
    #               nomatch=FALSE)
    #dat$Ratio[idx1] <- as.character(exDat$idColsAdj[idx1,c(4)])
    dat$Ratio <- as.factor(dat$Ratio)
    
    
    
    myYLim = plotInfo$myYLimMA
    myXLimMA = plotInfo$myXLimMA
    xlimMArange <- myXLimMA[2]-myXLimMA[1]
    if (xlimMArange > (plotInfo$myYLim[2] - plotInfo$myYLim[1])){
        myXLimMA <- c(plotInfo$myYLim[1],plotInfo$myYLim[2])
    }
    
    filenameRoot = sampleInfo$filenameRoot
    sample1 = exDat$sampleNames[1]
    sample2 = exDat$sampleNames[2]
    
    
    #idCols = exDat$idColsAdj
    r_m.res = exDat$Results$r_m.res
    
    #cutoffs = LODR.annot.ERCC$cutoffs
    cutoffs = exDat$Results$LODR.annot.ERCC$AbundCutoffs
    
    FCcode = erccInfo$FCcode
    legendLabels=erccInfo$legendLabels
    
    
    spikeFraction = exDat$spikeFraction
    
    r_m.mn = r_m.res$r_m.mn
    r_m.mnse = r_m.res$r_m.mnse
    
    if(is.null(r_m.mn)) {
        r_mAdjust = FALSE
        
    }
    #theme_update(legend.justification=c(1,0), legend.position=c(1,0))
    
    
    colScale <- plotInfo$colScale
    fillScale <- plotInfo$fillScale
    
    
    if (length(cutoffs) > 0){
        cat("\nLODR estimates are available to code ratio-abundance plot\n")  
    }
    
    maStats <- function(x, c1, c2){
        c(mean(log2(x[c1])-log2(x[c2])),sd(log2(x[c1])-log2(x[c2])),log2(mean(x)))
    } 
    
    totCol <- ncol(dat[-c(1:2)])
    
    if(odd(totCol)) stop("Uneven number of replicates for the two sample types")
    
    maStatDat <- data.frame(t(apply(dat[-c(1:2)],1,maStats, 
                                    c1 = c(1:(totCol/2)),
                                    c2 = c(((totCol/2)+1):totCol))))
    colnames(maStatDat) <- c("M.Ave","M.SD","A")
    maDatAll <- cbind(dat, maStatDat)
    
    maDatAll <- maDatAll[which(is.finite(maDatAll$M.Ave)),]
    
    ### Now subset and continue with just ERCCs
    maData <- subset(maDatAll, Ratio != "Endo")
    maData$Ratio <- factor(as.character(maData$Ratio),levels=FCcode$Ratio)
    maData$Feature <- as.factor(as.character(maData$Feature))
    #countPair <- subset(exDat$expressDat_l, subset = Ratio != "Endo")
    #countPair$Ratio <- as.factor(as.character(countPair$Ratio))
    
    maData$Nominal = FCcode$FC[1]
    for (i in 2:nlevels(FCcode$Ratio)){
        maData$Nominal[which(maData$Ratio == FCcode$Ratio[i])] = FCcode$FC[i]
    }
    
    
    if(r_mAdjust == TRUE){
        maData$Empirical = maData$Nominal/exp(r_m.mn)  
    }else{
        maData$Empirical = maData$Nominal
    }
    
    
    
    cutERCCs = unique(maData$Feature[which(is.na(maData$M.SD))])
    
    if(length(cutERCCs) != 0){
        cat(paste("\nThese ERCCs were not included in the ratio-abundance",
                  "plot, \n",
                  "because not enough non-zero replicate measurements of",
                  "these \n",
                  "controls were obtained for both samples:\n"))
        
        cat(paste(as.character(cutERCCs),collapse="\n"))
        
        maData = maData[-(which(maData$Feature %in% cutERCCs)),]
        searchCut = paste(cutERCCs,collapse="|")
        maData = maData[-(grep(searchCut,maData$Feature)),]
    }
    
    #maDataAveSD$Feature = as.factor(as.character(maDataAveSD$Feature))  
    maData$Feature = as.factor(as.character(maData$Feature))
    #write.csv(maData,file = paste(filenameRoot,"maDataFinite.csv",sep = "."))
    
    #avexlabel = exDat$ERCCxlabelAve
    if(sampleInfo$datType == "count"){
        avexlabel = "Log2 Average of Normalized Counts"
        ymalabel = "Log2 Ratio of Normalized Counts"
    }
    if(sampleInfo$datType == "array"){
        avexlabel = "Log2 Average of Normalized Intensity"
        ymalabel = "Log2 Ratio of Normalized Intensity"
        # myXLimMA = c(min(maData$A)-1, max(maData$A)+1)
        # exDat$plotInfo$myXLimMA <- myXLimMA
    }
    xlabel = xlab(avexlabel)
    
    # Estimate SD for all ratio measurements
    # For the ERCCs that are plotted take the log ratio data and subtract the 
    ## log nominal ratio, this will normalize the data
    normRats = maData$M.Ave - log2(maData$Nominal)
    #print(normRats)
    
    # Take the sd of all of the normalized log ratio data to find a global SD 
    ##for the ERCCs at this site
    sdGlobal = sd(normRats)
    #cat("\nGlobal Ratio SD for this sample pair is: ")
    #cat(sdGlobal, "\n")
    
    ratVarDat <- maData
    #ratVarDat <- subset(ratVarDat, A > 0)
    #ratVarDat$A <- log2((2^(ratVarDat$A.Ave))/(spikeFraction))
    
    # set null values to appease R CMD Check
    Ratio <- A <- M.Ave <- M.SD <- LODR <- Nominal <- Empirical <- NULL
    if (length(cutoffs)>0){
        maData$LODR = "below"
        FCcodeC = FCcode[-c(which(FCcode$FC == 1)),]
        for (i in 1:length(cutoffs)){
            maData$LODR[which((maData$A > cutoffs[i])&
                                  (maData$Ratio ==
                                       FCcodeC$Ratio[i]))] = "above"
        }
        
        maData$LODR <- as.factor(maData$LODR)
        if (is.null(r_m.mn)){
            rm_dat = data.frame("Mean r_m" = as.character("'No Estimate'"),
                                "Weighted SE" = as.character(" "))
        }else{
            rm_dat = data.frame("Mean r_m" = signif(r_m.mn,digits = 4),
                                "Weighted SE" = signif(r_m.mnse,digits = 4))  
        }
        
        colnames(rm_dat) <- c("'Weighted Mean'",
                              "'(+/-) Weighted Standard Error'")
        rownames(rm_dat) <- expression(log(r[m]))
        
        # Plot ratio-signal data coding points below the LODR with open circles
        maPlot <- ggplot(maData, aes(x = A, y = M.Ave) ) + 
            geom_point(data = subset(maDatAll, maDatAll$Ratio == "Endo"),
                       aes(x = A, y = M.Ave),
                       colour = "grey80", alpha = 0.5) +
            geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                              colour = Ratio),size = 1,alpha = alphaPoint) + 
            geom_point(aes(colour = Ratio),size = 5, alpha = alphaPoint) +
            geom_point(data = subset(maData, (LODR == "below")),
                       colour = "white",size = 2.5) + 
            geom_hline(aes(yintercept = log2(Nominal), colour = Ratio), 
                       alpha = 0.7) +
            geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                       size = 1, linetype = "longdash") + 
            ylab(ymalabel) + xlabel + 
            coord_cartesian(xlim = myXLimMA, ylim = myYLim) + colScale + 
            annotation_custom(tableGrob(rm_dat,parse=TRUE, 
                                        gpar.corefill = gpar(fill = "grey85",
                                                             col = "white"), 
                                        gpar.rowfill = gpar(fill = "grey80",
                                                            col = "white"),
                                        gpar.colfill = gpar(fill = "grey80",
                                                            col = "white")), 
                              #xmin = quantile(maData$A,probs=0.25),
                              #xmax = max(maData$A),
                              ymin = (myYLim[2]) - 0.25*myYLim[2], 
                              ymax = myYLim[2]) + 
            scale_y_continuous(breaks = seq(myYLim[1],myYLim[2],1))+ theme_bw()+
            theme( legend.justification = c(1,0),legend.position=c(1,0)) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank())
        
    }else{
        # Plot ratio-signal data without LODR coding
        maPlot <- ggplot(maData, aes(x = A, y = M.Ave, colour = Ratio) ) + 
            geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD),
                          size = 1, alpha = alphaPoint) + 
            geom_point(size = 5, alpha = alphaPoint) + 
            geom_hline(aes(yintercept = log2(Nominal), colour = Ratio), 
                       alpha = 0.7) +
            geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                       size = 1, linetype = "longdash") + 
            ylab(ymalabel) + xlabel + colScale+ 
            coord_cartesian(xlim = myXLimMA, ylim = myYLim) + theme_bw()+
            theme( legend.justification = c(1,0),legend.position=c(1,0)) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank())
    }
    
    ratioVarPlot <- ggplot(maData) + geom_violin(aes(x = Ratio, 
                                                     y = M.SD, 
                                                     fill = Ratio), 
                                                 alpha = alphaPoint) +
        ylab("SD of Log2 Ratios") + colScale + fillScale + theme_bw()
    
    
    exDat$Figures$maPlot <- maPlot
    #exDat$Results$ratVarDat <- ratVarDat
    #exDat$Results$modRatVar <- stdevCoef
    exDat$Results$maDatAll <- maDatAll
    #exDat$Figures$ratioSDPlot <- sdRatioplotFit
    #exDat$Figures$ratioSDPlot <- ratioVarPlot
    
    return(exDat)
}