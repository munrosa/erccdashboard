testDEArray <- function(exDat){
    #library(limma)
    y <- exDat$Transcripts
    choseFDR <- exDat$sampleInfo$choseFDR
    if(is.null(choseFDR)){
        p.thresh <- exDat$Results$p.thresh
    }
    
    #normy <- sweep(y[-c(1)], 2, exDat$normFactor/1000,"/")
    #y <- cbind(y[c(1)],normy)
    erccInfo <- exDat$erccInfo
    sampleInfo <- exDat$sampleInfo
    row.names(y) <- make.names(y$Feature, unique=TRUE)
    y <- y[-c(1)]
    
    yERCC <- y[grep("ERCC.00",row.names(y)),]
    yAll <- y[-grep("ERCC.00",row.names(y)),]
    
    # adjust for r_m before hypothesis testing
    if(!is.null(exDat$Results$r_m.res$r_m.mn)){
        adj <- exp(exDat$Results$r_m.res$r_m.mn)
        
        yERCC[c(1:(ncol(yERCC)/2))] <- yERCC[c(1:(ncol(yERCC)/2))]*adj
    }
    
    y <- rbind(yERCC, yAll)
    
    #y <- yERCC
    if(is.null(exDat$normFactor)&(exDat$sampleInfo$isNorm==TRUE)){
        cat("\nisNorm is TRUE, array data is already normalized\n")
        ynorm <- y
    }else{
        ynorm <- sweep(y, 2, exDat$normFactor, "/")
        
    }  
    
    ylog <- log2(ynorm)
    
    if(odd(ncol(ynorm))) stop(paste("\nUneven number of replicates for the",
                                    "two sample types\n"))
    
    design <- cbind(Grp1=1,Grp1vs2=c(rep(x=1,times=ncol(ynorm)/2), 
                                     rep(x=0,times=ncol(ynorm)/2)))
    
    fit <- lmFit(ylog,design)
    
    fit <- eBayes(fit)
    
    res <- topTable(fit,sort.by="none",number = dim(ylog)[1],coef = 2)
    
    ### generate qvals
    if(!is.null(choseFDR)){
        
        pval <- res$P.Value
        
        res$qvals <- qvalue(pval)$qvalues
        
        if(any(res$qvals<choseFDR)){
            p.thresh<-max(res$P.Value[res$qvals<choseFDR])
        }
    }
    res$Feature <- row.names(res)
    
    erccFC <- data.frame(Feature = erccInfo$idColsSRM$Feature, 
                         FC = round(erccInfo$idColsSRM$Conc1/
                                        erccInfo$idColsSRM$Conc2,digits=3))
    
    yMns <- data.frame(Feature = row.names(y), MnSignal = rowMeans(y))
    mergedRes <- merge(res,yMns,by="Feature")
    
    ERCCres<- mergedRes[grep("ERCC-",mergedRes$Feature),]
    Endores <- mergedRes[-grep("ERCC-",mergedRes$Feature),]
    
    ercc.pval.res <- data.frame( Feature = ERCCres$Feature,
                                 MnSignal = ERCCres$MnSignal,
                                 Pval = ERCCres$P.Value,
                                 Fold = erccFC$FC[match(ERCCres$Feature, 
                                                        erccFC$Feature,
                                                        nomatch=0)])
    
    write.csv(ercc.pval.res, paste(sampleInfo$filenameRoot, "ERCC Pvals.csv"),
              row.names = FALSE)
    
    
    endo.pval.res <- data.frame( Feature = Endores$Feature, 
                                 MnSignal = Endores$MnSignal, 
                                 Pval = Endores$P.Value,
                                 Fold = NA)
    
    all.pval.res <- rbind(ercc.pval.res, endo.pval.res)
    
    write.csv(all.pval.res, paste0(sampleInfo$filenameRoot, ".All.Pvals.csv"),
              row.names = FALSE)
    exDat$Results$limma.res <- all.pval.res
    exDat$Results$ERCC.pval <- ercc.pval.res
    exDat$Results$p.thresh <- p.thresh
    return(exDat)
}