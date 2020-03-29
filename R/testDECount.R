testDECount <- function(sampleInfo, exDat, cnt = cnt, info = info){
    # Pull info from exDat
    erccInfo <- exDat$erccInfo
    plotInfo <- exDat$plotInfo
    filenameRoot = sampleInfo$filenameRoot
    legendLabels = sampleInfo$legendLabels
    FCcode = erccInfo$FCcode
    
    idCols = exDat$idColsAdj
    r_m.mn = exDat$Results$r_m.res$r_m.mn 
    repNormFactor = sampleInfo$repNormFactor 
    normFactor <- exDat$normFactor
    sample1 = exDat$sampleNames[1]
    sample2 = exDat$sampleNames[2]
    
    colScale <- plotInfo$colScale
    fillScale <- plotInfo$fillScale
    
    ## Organize the count table
    cnt = unique(cnt)
    Features = make.names(cnt$Feature,unique=TRUE)
    Features = gsub(".","-", Features, fixed = TRUE)
    rownames(cnt)<-Features
    cnt<-as.matrix(cnt[,-1])
    
    if(odd(ncol(cnt))) stop(paste("Uneven number of replicates for",
                                  "the two sample types"))
    
    colnames(cnt)<-paste(rep(c(sample1,sample2),
                             each=ncol(cnt)/2),
                         c(1:(ncol(cnt)/2),1:(ncol(cnt)/2)),sep="")
    
    ## Get ERCC names
    ERCC<-rownames(cnt[substr(rownames(cnt),1,5)=="ERCC-",])
    
    ## Specify Sample (A or B)
    trt<-rep(1:2,each=ncol(cnt)/2)
    
    ## Get design.list for edgeR
    design.list<-list(trt,rep(1,ncol(cnt)))
    
    ## Compute offset (e.g. total counts, 75% quantile, TMM, etc)
    if(is.null(normFactor)){
        log.offset <- c(rep(0,ncol(cnt)))
    }else{
        log.offset<- log(normFactor)
    }
    
    cat("\nShow log.offset\n")
    cat(log.offset,"\n")
    
    ERCC.FC = idCols[c(1,4)];rownames(ERCC.FC)<-ERCC.FC[,1]
    
    ERCC.FC$NumRatio <- NA
    FCcode$Ratio <- as.factor(FCcode$Ratio)
    for (i in 1:nlevels(FCcode$Ratio)){
      print(i)
        ERCC.FC$NumRatio[which(ERCC.FC$Ratio == 
                                   FCcode$Ratio[i])]=FCcode$FC[i]  
    }
    ERCC.Ratio = ERCC.FC[c(1,2)]
    ERCC.FC = ERCC.FC[-c(2)]
    
    group <- as.factor(trt)
    d <- DGEList(counts=cnt,group=group)
    # use log.offset for the library size
    d$samples$lib.size <- exp(log.offset) 
    
    #Dispersion trend
    design <- model.matrix(~group)
    d1 <- estimateGLMCommonDisp(d,design,verbose=TRUE)
    d1 <- estimateGLMTrendedDisp(d1,design)
    d1 <- estimateGLMTagwiseDisp(d1,design)
    
    #######################################################
    ###  Simulate ERCC data from negative binomial fit ####
    ###  (to be used for sim-based LODR)               ####
    #######################################################
    # Define function simcnt.lodr
    simcnt.lodr<-function(cnt,disp,trt,fold,log.offset){
        #### 'cnt' is the matrix of counts for endogenous genes
        #### 'disp' is the central trend fitted to the estimated dispersions 
        ####  (vs. expression)
        #### 'trt' is a vector specifying treatments for each column in 'cnt'
        #### 'fold' is the desired fold change (trt 2/trt 1)                  
        #### 'log.offset' is used to account for differences in library size
        
        ### mimick every 797th gene (roughly), when sorted by total count
        sim.ind<-round(nrow(cnt)*c(1:49)/50)  
        norm.cnt<-t(t(cnt)/exp(log.offset))
        sim.mn<-matrix(sort(rowMeans(norm.cnt))[sim.ind],
                       nrow=length(sim.ind), ncol=ncol(cnt))
        sim.mn<-sim.mn*2/(1+fold)
        sim.mn[,trt==2]<-sim.mn[,trt==2]*fold
        sim.mn<-t(t(sim.mn)*exp(log.offset))
        
        sim.disp<-disp[order(rowMeans(norm.cnt))][sim.ind]
        size<-matrix(1/sim.disp,length(sim.disp),ncol(cnt))
        simcnt<-matrix(rnbinom(length(sim.disp)*ncol(cnt),
                               mu=sim.mn,size=size),length(sim.disp),ncol(cnt))
        rownames(simcnt)<-paste("Sim",fold,"Fold",1:length(sim.ind),sep="")
        return(simcnt)
    }
    
    #### Simulate data for each of the fold changes used in the ERCCs
    simcnt<-NULL
    for(fold in unique(ERCC.FC[!is.na(ERCC.FC[,2]),2])){
        simcnt<-rbind(simcnt,simcnt.lodr(cnt,disp=d1$trended.dispersion,
                                         trt=trt,fold=fold,
                                         log.offset=log.offset))
    }
    
    ##### Analyze combination of observed and simulated data with edgeR
    group <- as.factor(trt)
    d2 <- DGEList(counts=rbind(cnt,simcnt),group=group)
    d2$samples$lib.size <- exp(log.offset)
    
    design <- model.matrix(~group)
    d2 <- estimateGLMCommonDisp(d2,design,verbose=TRUE)
    d2 <- estimateGLMTrendedDisp(d2,design)
    d2 <- estimateGLMTagwiseDisp(d2,design)
    
    NBdisp<-d2$tagwise.dispersion
    names(NBdisp)<-rownames(rbind(cnt,simcnt))
    NBdisptrend<-d2$trended.dispersion
    names(NBdisptrend)<-rownames(rbind(cnt,simcnt))
  
    use.fit <- glmQLFit(d2, design)
    
    qlf.res <- glmQLFTest(use.fit)
    
    use.res <- qlf.res
    
    ### Collect results for simulated data, to be passed along to LODR function
    Feature <- row.names(simcnt)
    MnSignal <- as.numeric(rowMeans(simcnt))
    Pval <- use.res$table$PValue[-(1:nrow(cnt))]
    LogPval <- log10(use.res$table$PValue)[-(1:nrow(cnt))]
    F.stat <- use.res$table$F[-(1:nrow(cnt))]
    Fold <- rep(unique(ERCC.FC[!is.na(ERCC.FC[,2]),2]),each=49)
    sim.pval.res <- data.frame(Feature, MnSignal, Pval, LogPval, F.stat, Fold)
    colnames(sim.pval.res)<-c("Feature","MnSignal","Pval","LogPval","F.stat",
                               "Fold")
    rownames(sim.pval.res)<-rownames(simcnt)
    
    write.csv(sim.pval.res[-c(4,5)],file=paste(filenameRoot,"Sim Pvals.csv"),
    row.names = FALSE)
    
    ## remove results for simulated data (using indexing with 1:nrow(cnt))
    use.fit2<-use.fit
    use.fit2$coefficients <-use.fit2$coefficients[1:nrow(cnt),]
    use.fit2$fitted.values <- use.fit2$fitted.values[1:nrow(cnt),]
    use.fit2$deviance <- use.fit2$deviance[1:nrow(cnt)]
    use.fit2$counts <- use.fit2$counts[1:nrow(cnt),]
    use.fit2$unshrunk.coefficients <- use.fit2$unshrunk.coefficients[1:nrow(cnt),]
    use.fit2$df.residual <- use.fit2$df.residual[1:nrow(cnt)]
    use.fit2$offset <- use.fit2$offset[1:nrow(cnt),]
    use.fit2$dispersion <- use.fit2$dispersion[1:nrow(cnt)]
    use.fit2$AveLogCPM <- use.fit2$AveLogCPM[1:nrow(cnt)]
    use.fit2$df.residual.zeros <- use.fit2$df.residual.zeros[1:nrow(cnt)]
    use.fit2$var.post <- use.fit2$var.post[1:nrow(cnt)]
    use.fit2$var.prior <- use.fit2$var.prior[1:nrow(cnt)]
    
    use.res2 <- glmQLFTest(use.fit2)
    
    ###################################
    #### Examine results for ERCCs ####
    ###################################
    pvals<-use.res2$table$PValue
    names(pvals)<-rownames(cnt)
    ERCC.pvals<-pvals[ERCC]
    
    ## Reanalyze ERCC transcripts using adjusted offsets to center fold 
    ## change estimates
    adj <- r_m.mn #### Use r_m estimated from NegBin GLM
    use.fit.adj <- glmQLFit(cnt[ERCC,], design, dispersion = NBdisptrend[ERCC],
                            offset = log.offset - rep(c(adj,0),each=ncol(cnt)/2))
    est.FC.adj<-use.fit.adj$coefficients[ERCC,2]
    
    use.fit3<-use.fit2
    # Substitute ERCC centered data into full use.fit3 structure
    use.fit3$coefficients[ERCC,] <- use.fit.adj$coefficients[ERCC,]
    use.fit3$unshrunk.coefficients[ERCC,] <- use.fit.adj$unshrunk.coefficients[ERCC,]
    use.fit3$var.prior[ERCC] <- use.fit.adj$var.prior[ERCC]
    use.fit3$var.post[ERCC] <- use.fit.adj$var.post[ERCC]
    
    # deal with CompressedMatrix format to add the new ercc offsets...
    expandedfit3 <- expandAsMatrix(use.fit3$offset)
    expandedfit3[1:length(ERCC),] <- expandAsMatrix(use.fit.adj$offset)
    recompress <- makeCompressedMatrix(expandedfit3)
    use.fit3$offset <- recompress
    
    use.res.adj<-glmQLFTest(use.fit3)
     
    # collect the results for plotting and writing a table
    Feature <- row.names(use.res.adj$table)
    MnSignal <- as.numeric(rowMeans(cnt))
    #replace with edgeR results
    Pval <- use.res.adj$table$PValue
    LogPval <- log10(use.res.adj$table$PValue)
    Qval <- qvalue(Pval)$qvalues
    F.stat <- use.res.adj$table$F
    
    names(Qval)<-names(F.stat)<-names(LogPval)<-names(Pval)<-rownames(cnt)
    
    Log2FC <- use.res.adj$table$logFC
    allDE.res <- data.frame(Feature = names(pvals),
                               MnSignal = rowMeans(cnt), 
                               Fold = 
                                   c(ERCC.FC[ERCC,2], rep(x=NA, length.out=
                                                              (length(pvals) -
                                                                   (length
                                                                    (ERCC.FC[
                                                                        ERCC,2]
                                                                    ))))), 
                               Log2Rat=Log2FC, Pval=Pval,
                               qvals=Qval, log.pvals=LogPval, F.stat=F.stat)
    
    write.csv(allDE.res[c(1,2,5,3)],
              file=paste0(filenameRoot, ".All.Pvals.csv"), row.names = FALSE)
    
    ERCC.Pval.adj<-Pval[ERCC]
    ERCC.F.stat.adj<-F.stat[ERCC]
    ERCC.LogPval.adj<-LogPval[ERCC]
    
    ### Collect results for ERCCs; to be passed along to LODR function
    ERCC.pval.res<-data.frame(row.names(cnt[ERCC,]),rowMeans(cnt[ERCC,]),
                              ERCC.Pval.adj,ERCC.LogPval.adj,ERCC.F.stat.adj,
                              ERCC.FC[ERCC,2])
    colnames(ERCC.pval.res)<-c("Feature","MnSignal","Pval","LogPval","F.stat","Fold")
    #print(str(pval.res))
    row.names(ERCC.pval.res) <- NULL
    write.csv(ERCC.pval.res[-c(4,5)],file=paste(filenameRoot,"ERCC Pvals.csv"),
              row.names = FALSE)
    cat("Finished DE testing")
    
    exDat$Results$allDE.res <- allDE.res
    exDat$Results$ERCCpvals <- ERCC.pval.res
    
    ### Using QLDisp code from edgeR to create similar ggplot with use.fit
      glmfit <- use.fit
      xlab = "Average Log2 CPM (Counts per million)"
      ylab = "Squeezed QL Dispersion Estimates (Quarter-Root Mean Deviance)"
      
      A <- glmfit$AveLogCPM
    
      if (is.null(A)) 
        A <- aveLogCPM(glmfit)
      s2 <- glmfit$deviance/glmfit$df.residual.zeros
      if (is.null(glmfit$var.post)) {
        stop("need to run glmQLFit before getting QL dispersion estimates")
      }
      
      squeezedDisp <- data.frame(A = A, Dispersion = sqrt(sqrt(glmfit$var.post)))
      dispERCC <- squeezedDisp[ERCC,]
      dispERCC$Ratio <- ERCC.Ratio[ERCC,2]
      
      if (length(glmfit$var.prior) == 1L) {
        trendPlot <- geom_abline(yintercept = sqrt(sqrt(glmfit$var.prior)))
      } else {
        o <- order(A)
        dispTrend <- data.frame(Atrend = A[o], Dispersion = sqrt(sqrt(glmfit$var.prior[o])) )
        trendPlot <-  geom_line(data = dispTrend, aes(Atrend, Dispersion))
      }
      quasiDispPlot <- ggplot() + geom_point(data = squeezedDisp, aes(x = A, y = Dispersion),
                                             colour = "grey80", size = 5,
                                             alpha = 0.6) +
        geom_point(data = dispERCC, aes(x = A, y = Dispersion,colour = Ratio), 
                   size = 5, alpha = 0.6) + xlab(xlab) + ylab(ylab) +
        trendPlot +
        #scale_x_log10() + 
        colScale + theme_bw() +
        theme(legend.justification=c(1,1), legend.position=c(0.9,0.9)) 
  
    exDat$Figures$dispPlot <- quasiDispPlot
    exDat$Results$simcnt <- simcnt 
    cat("\nFinished examining dispersions\n")
    
    return(exDat)
    ### end Edit Sarah Munro 20140216
}