#' Estimate the mRNA fraction differences for the pair of samples using 
#' replicate data 
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @details
#' This is the first function to run after an exDat structure is initialized
#' using initDat, because it is needed for all additional analysis. An r_m of
#' 1 indicates that the two sample types under comparison have 
#' similar mRNA fractions of total RNA. The r_m estimate is used to adjusted 
#' the expected ERCC mixture ratios in this analysis and may indicate a need for
#' a different sample normalization approach.
#' 
#' @examples
#' data(SEQC.Example)
#' 
#' exDat <- initDat(datType="count", isNorm = FALSE, exTable=MET.CTL.countDat, 
#'                  filenameRoot = "testRun",sample1Name = "MET",
#'                  sample2Name = "CTL", erccmix = "RatioPair", 
#'                  erccdilution = 1/100, spikeVol = 1, totalRNAmass = 0.500,
#'                  choseFDR = 0.1)
#'
#' exDat <- est_r_m(exDat)
#' 
#' @export

est_r_m <- function(exDat){
    cat("\nCheck for sample mRNA fraction differences(r_m)...\n")
    
    cnt = exDat$Transcripts
    sampleInfo <- exDat$sampleInfo
    plotInfo <- exDat$plotInfo
    erccInfo <- exDat$erccInfo
    
    site <- sampleInfo$siteName
    avexlabel <- plotInfo$ERCCxlabelAve
    myXLim <- plotInfo$myXLim
    idCols <- exDat$idCols
    # requires that the sample1 columns are first in the table -> 
    # need to add a stopifnot statement for this 
    
    colScale <- plotInfo$colScale
    fillScale <- plotInfo$fillScale
    
    dat = unique(cnt)
    Features = make.names(dat$Feature,unique=TRUE)
    Features = gsub(".","-", Features, fixed = TRUE)
    rownames(dat)<-Features; dat<-as.matrix(dat[,-1])
    
    colnames(dat)<-paste(rep(c(exDat$sample1,exDat$sample2),
                             each=ncol(dat)/2),c(1:(ncol(dat)/2),
                                                 1:(ncol(dat)/2)),sep="")
    
    if (sampleInfo$isNorm == TRUE){
        cat("\nData is normalized, no r_m estimation\n")
        exDat$Results$r_m.res <- list(r_m.mn = NULL, r_m.mnse =NULL)
        exDat$Results$r_mDat <- NULL
        return(exDat)
    }  
    if (sampleInfo$datType == "array"){
        cat("\ndatType is array, using 1:1 ERCC controls for r_m estimate\n")
        # Pull out all ERCCs from dat except for the 1:1 ERCCs
        ERCCcut <- idCols$Feature[idCols$Ratio != "1:1"]
        dat <- dat[-(which(rownames(dat) %in% as.character(ERCCcut))),]
    }
    ## Get ERCC names
    ERCC<-rownames(dat[substr(rownames(dat),1,5)=="ERCC-",])
    
    ## Specify Sample (A or B)
    trt<-rep(1:2,each=ncol(dat)/2)
    design.list<-list(trt,rep(1,ncol(dat)))
    
    ## Compute offset (e.g. total counts, 75% quantile, TMM, etc) could modify
    # in future to enable other methods
    #log.offset<-log(colSums(dat))
    log.offset <- NULL
    if(sampleInfo$isNorm == TRUE){
        log.offset = NULL
    }else{
        log.offset = log(exDat$normFactor)
    }
    
    ######################################################
    #### Estimate r_m from each ERCC using NegBin GLM ####
    ######################################################
    
    ERCC.FC = idCols[c(1,4)];rownames(ERCC.FC)<-ERCC.FC[,1]
    ERCC.FC$NomRatio <- 1
    for (i in 1:nlevels(erccInfo$FCcode$Ratio)){
        ERCC.FC$NomRatio[which(ERCC.FC$Ratio == as.character(
            erccInfo$FCcode$Ratio[i]))] = erccInfo$FCcode$FC[i]  
    }
    ERCC.Ratio = ERCC.FC[c(1,2)]
    
    ERCC.FC = ERCC.FC[-c(2)]
    
    #library(MASS)
    
    r_m<-NULL
    cat("\nNumber of ERCC Controls Used in r_m estimate\n")
    cat(length((1:nrow(dat))[substr(rownames(dat),1,5)=="ERCC-"]),"\n")
    
    for( i in (1:nrow(dat))[substr(rownames(dat),1,5)=="ERCC-"]){ 
        
        #Obtain estimated log fold-change and standard error for each ERCC
        
        r_m<-rbind(r_m,summary(suppressWarnings(glm.nb(dat[i,]~as.factor(trt) +
                                                           offset(log.offset)
        )))$coefficients[2,1:2])
    }
    
    #### Collect Per ERCC r_m estimates and se in data frame
    colnames(r_m)<-c("r_m.hat","r_m.se")
    rownames(r_m)<-rownames(dat)[substr(rownames(dat),1,5)=="ERCC-"]
    r_m <- data.frame(r_m)
    
    #### Add nominal log fold change to results    
    r_m$nominal<- - log(ERCC.FC[rownames(r_m),2])
    
    #### Make 95% CIs based on t-distribution
    dft <- ncol(dat) - 2  
    quant<-qt(.975, dft )
    
    #### Plot log r_m estimates and corresponding 95% CIs
    r_m$Ratio = ERCC.Ratio$Ratio[match(row.names(r_m), ERCC.Ratio$Feature)]
    # Subset data frame to remove missing ERCC controls
    r_m <- subset(r_m, is.finite(Ratio))
    
    #### Get global site r_m estimate and sd
    #  # weighted mean of the individual means
    
    #  # standard deviation, sigma of the weighted mean
    
    
    ### Try new version SM code added 20140410
    w_vec <- 1/(r_m[,2]^2)
    r_m.mn <- sum((r_m[,1]-r_m[,3])*w_vec)/sum(w_vec)
    
    
    len_w_vec <- length(w_vec[w_vec > 0])
    #print(len_w_vec)
    r_m.mnse <- sqrt((sum(w_vec*(((r_m[,1]-r_m[,3])-
                                      r_m.mn)^2)))/(((len_w_vec-1)*
                                                         sum(w_vec))/
                                                        (len_w_vec)))
    
    # add Feature column to left side of r_m matrix
    r_m = cbind(as.character(row.names(r_m)), r_m)
    names(r_m)[1]= "Feature" 
    
    r_m.mnlog = exp(r_m.mn)
    
    idCols$Conc2 = r_m.mnlog*idCols$Conc2
    idColsAve = idCols
    
    idColsAve = idCols[match(r_m$Feature,idCols$Feature),]
    
    r_m$AveConc = log2((idColsAve$Conc1+ idColsAve$Conc2)/2)
    
    r_m = r_m[order(r_m$AveConc),]
    
    
    r_m$ymax = r_m$r_m.hat - r_m$nominal + (quant)*r_m$r_m.se
    r_m$ymin = r_m$r_m.hat - r_m$nominal - (quant)*r_m$r_m.se
    
    textDat = subset(r_m, (r_m.mn<(r_m.hat - nominal - 
                                       (quant*r_m.se))) | 
                         (r_m.mn>(r_m.hat - nominal + (quant*r_m.se))))
    
    cat("\nOutlier ERCCs for GLM r_m Estimate:\n")
    if (length(textDat$Feature) > 0){
        outlierERCCs <- as.character(textDat$Feature)
        for (j in seq(from=1,to=length(outlierERCCs),by=5)){
            k = j+4
            if (k > length(outlierERCCs)) k = length(outlierERCCs)
            cat(outlierERCCs[j:k])
            cat("\n")
        } 
        
    }else{
        cat("None","\n")
    }
    
    
    cat(paste("\nGLM log(r_m) estimate:\n"))
    cat(r_m.mn,"\n")
    
    cat("\nGLM log(r_m) estimate weighted s.e.:\n")
    cat(r_m.mnse,"\n")
    
    
    #cat(paste(site,"\nGLM r_m estimate:\n"))
    #cat(exp(r_m.mn),"\n")
    
    #cat("\nGLM r_m weighted s.e.\n")
    #cat(exp(r_m.mnse),"\n")
    
    # appease CMD check
    
    Ratio <- r_m.hat <- nominal <- r_m.se <- AveConc <- ymin <- ymax <- NULL
    Feature <- NULL
    
    if (nrow(textDat)>1){
        plotSiter_m = ggplot(r_m, aes(x = AveConc, y = r_m.hat - nominal, 
                                      colour =  Ratio)) + 
            geom_point(size = 6, alpha = 0.7) + 
            geom_errorbar(aes(ymin = ymin,ymax = ymax),alpha = 0.7) + 
            coord_cartesian(xlim = myXLim, ylim = c(-1.5,1.5)) + 
            xlab(avexlabel) + 
            ylab(expression(log(r[m]))) + 
            geom_text(data = textDat, aes(x = AveConc, y = (r_m.hat - nominal),
                                          label = gsub("ERCC-00","",Feature)),
                      colour="black", size=6,show_guide=FALSE,angle=90) + 
            geom_hline(yintercept = r_m.mn) + colScale + theme_bw()
        theme(legend.justification=c(1,0), legend.position=c(1,0))  
    }else{
        plotSiter_m = ggplot(r_m, aes(x = AveConc, y = r_m.hat - nominal, 
                                      colour =  Ratio)) + 
            geom_point(size = 6, alpha = 0.7) + 
            geom_errorbar(aes(ymin = ymin,ymax = ymax),alpha = 0.7) + 
            coord_cartesian(xlim=myXLim, ylim=c(-1.5,1.5)) + xlab(avexlabel) + 
            ylab(expression(log(r[m]))) + geom_hline(yintercept = r_m.mn) + 
            colScale + theme_bw() + 
            theme(legend.justification=c(1,0), legend.position=c(1,0))
    }
    
    exDat$idColsAdj <- idCols
    exDat$Results$r_m.res <- list(r_m.mn = r_m.mn, r_m.mnse = r_m.mnse)
    
    exDat$Figures$r_mPlot <- plotSiter_m
    exDat$Results$r_mDat <- r_m  
    return(exDat)
}