testDECount<- function(sampleInfo, exDat, cnt = cnt, info = info){
    #library(QuasiSeq)
    #library(DESeq)
    #library(edgeR)
    erccInfo <- exDat$erccInfo
    plotInfo <- exDat$plotInfo
    filenameRoot = sampleInfo$filenameRoot
    legendLabels = sampleInfo$legendLabels
    FCcode = erccInfo$FCcode
    #totalSeqReads = sampleInfo$totalSeqReads
    
    
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
        #stop(cat("\nlibe size normalization is missing\n"))
        
        #log.offset<-log(colSums(cnt))
        #cat("\nUsing Mapped Reads\n")
        #cat(colSums(cnt),"\n")  
    }else{
        #log.offset<- log(repNormFactor)
        log.offset<- log(normFactor)
        #cat("\nUsing repNormFactor\n")
        #cat(repNormFactor,"\n")
        
    }
    
    cat("\nShow log.offset\n")
    cat(log.offset,"\n")
    
    ERCC.FC = idCols[c(1,4)];rownames(ERCC.FC)<-ERCC.FC[,1]
    
    ERCC.FC$NumRatio <- 1
    for (i in 1:nlevels(FCcode$Ratio)){
        ERCC.FC$NumRatio[which(ERCC.FC$Ratio == 
                                   as.character(FCcode$Ratio[i]))]=FCcode$FC[i]  
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
    
    ### Plot estimated dispersions
    xdat = ydat = xsim = ysim = xERCC = yERCC = NULL
    dispcnt = data.frame(xdat = rowMeans(rbind(cnt,simcnt)), ydat = NBdisp)
    dispSimcnt = data.frame(xsim = rowMeans(simcnt), 
                            ysim=NBdisp[-c(1:nrow(cnt))])
    dispERCCcnt = data.frame(xERCC = rowMeans(cnt)[ERCC],yERCC = NBdisp[ERCC])
    
    dispPlot = ggplot(dispcnt) + geom_point(aes(x = xdat, y = ydat)) + 
        scale_x_log10() + xlab("Average Count") + ylab("Estimated Dispersion") +
        geom_point(data = dispSimcnt,aes(x = xsim, y = ysim),colour = 3) + 
        geom_point(data = dispERCCcnt, aes(x = xERCC, y = yERCC),colour = 2) + 
        theme(legend.justification=c(1,1), legend.position=c(1,1))
    #print(dispPlot) 
    
    ### Run QuasiSeq
    use.fit<-QL.fit(rbind(cnt,simcnt), design.list,
                    log.offset=log.offset, 
                    NBdisp=NBdisptrend)# using trended dispersion from edgeR
    
    use.res<-QL.results(use.fit, Plot = FALSE)
    
    ### Collect results for simulated data; to be passed along to LODR function
    #   sim.pval.res<-cbind(row.names(simcnt), rowMeans(simcnt),
    #                       use.res$P.values$QLSpline[-(1:nrow(cnt))],
    #                       use.res$log.P.values$QLSpline[-(1:nrow(cnt))],
    #                       use.res$F.stat$QLSpline[-(1:nrow(cnt))],
    #                       rep(unique(ERCC.FC[!is.na(ERCC.FC[,2]),2]),each=49),
    #                       rep(use.res$d0[2],nrow(simcnt)))
    Feature <- row.names(simcnt)
    MnSignal <- as.numeric(rowMeans(simcnt))
    Pval <- use.res$P.values$QLSpline[-(1:nrow(cnt))]
    LogPval <- use.res$log.P.values$QLSpline[-(1:nrow(cnt))]
    F.stat <- use.res$F.stat$QLSpline[-(1:nrow(cnt))]
    Fold <- rep(unique(ERCC.FC[!is.na(ERCC.FC[,2]),2]),each=49)
    Den.df <- rep(use.res$d0[2],nrow(simcnt))
    sim.pval.res <- data.frame(Feature, MnSignal, Pval, LogPval, F.stat, Fold, 
                               Den.df)
    #colnames(sim.pval.res)<-c("Feature","MnSignal","Pval","LogPval","F.stat",
    #                           "Fold", "Den.df")
    #rownames(sim.pval.res)<-rownames(simcnt)
    
    write.csv(sim.pval.res[-c(4,5,7)],file=paste(filenameRoot,"Sim Pvals.csv"),
              row.names = FALSE)
    
    ## remove results for simulated data
    use.fit2<-use.fit
    use.fit2$LRT<-matrix(use.fit2$LRT[1:nrow(cnt),],nrow(cnt),1)
    use.fit2$phi.hat.dev<-use.fit2$phi.hat.dev[1:nrow(cnt)]
    use.fit2$phi.hat.pearson<-use.fit2$phi.hat.pearson[1:nrow(cnt)]
    use.fit2$mn.cnt<-use.fit2$mn.cnt[1:nrow(cnt)]
    use.fit2$NB.disp<-use.fit2$NB.disp[1:nrow(cnt)]
    use.fit2$fitted.values<-use.fit2$fitted.values[1:nrow(cnt),]
    use.fit2$coefficients<-use.fit2$coefficients[1:nrow(cnt),]
    
    use.res2<-QL.results(use.fit2, Plot = FALSE)
    
    ###################################
    #### Examine results for ERCCs ####
    ###################################
    pvals<-use.res2$P.values$QLSpline
    names(pvals)<-rownames(cnt)
    ERCC.pvals<-pvals[ERCC]
    
    rownames(use.fit2$coefficients)<-rownames(cnt)
    est.FC<-use.fit2$coefficients[ERCC,2]
    
    ## Reanalyze ERCC transcripts using adjusted offsets to center fold 
    ## change estimates
    adj <- r_m.mn #### Use r_m estimated from NegBin GLM
    use.fit.adj<-QL.fit(cnt[ERCC,],design.list,log.offset=
                            log.offset-rep(c(adj,0),each=ncol(cnt)/2),
                        NBdisp=NBdisptrend[ERCC])
    
    rownames(use.fit.adj$coefficients)<-ERCC;
    est.FC.adj<-use.fit.adj$coefficients[ERCC,2]
    
    use.fit3<-use.fit2
    rownames(use.fit3$LRT)<-rownames(cnt)
    # substitute the r_m adjusted LRT for the ERCCs
    use.fit3$LRT[ERCC,]<-use.fit.adj$LRT 
    
    use.res.adj<-QL.results(use.fit3, Plot = FALSE)
    
    pvals<-use.res.adj$P.values$QLSpline
    log.pvals<-use.res.adj$log.P.values$QLSpline
    F.stat<-use.res.adj$F.stat$QLSpline
    names(F.stat)<-names(log.pvals)<-names(pvals)<-rownames(cnt)
    
    qvals<-use.res.adj$Q.values$QLSpline
    logRatio <- function(x, c1, c2){
        log2(x[c1])-log2(x[c2])
    }
    totCol <- ncol(cnt)
    
    if(odd(totCol)) stop("Uneven number of replicates for the two sample types")
    
    ratioDat <- data.frame(t(apply(cnt,1,logRatio, 
                                   c1=c(1:(totCol/2)),
                                   c2=c(((totCol/2)+1):totCol))))
    
    colnames(ratioDat)<- "Log2Rat"
    quasiSeq.res <- data.frame(Feature = names(pvals),
                               MnSignal = rowMeans(cnt), 
                               Fold = 
                                   c(ERCC.FC[ERCC,2], rep(x=NA, length.out=
                                                              (length(pvals) -
                                                                   (length
                                                                    (ERCC.FC[
                                                                        ERCC,2]
                                                                    ))))), 
                               Log2Rat=ratioDat$Log2Rat, Pval=pvals,
                               qvals=qvals, log.pvals=log.pvals, F.stat=F.stat, 
                               den.df=rep(use.res.adj$d0[2], length(pvals)))
    
    write.csv(quasiSeq.res[c(1,2,5,3)],
              file=paste0(filenameRoot, ".All.Pvals.csv"), row.names = FALSE)
    
    ERCC.pvals.adj<-pvals[ERCC]
    ERCC.F.stat.adj<-F.stat[ERCC]
    ERCC.log.pvals.adj<-log.pvals[ERCC]
    
    ### Collect results for ERCCs; to be passed along to LODR function
    pval.res<-data.frame(row.names(cnt[ERCC,]),rowMeans(cnt[ERCC,]),
                         ERCC.pvals.adj,ERCC.log.pvals.adj,ERCC.F.stat.adj,
                         ERCC.FC[ERCC,2],rep(use.res.adj$d0[2],
                                             length(ERCC.pvals.adj)))
    colnames(pval.res)<-c("Feature","MnSignal","Pval","LogPval","F.stat","Fold",
                          "Den.df")
    #print(str(pval.res))
    row.names(pval.res) <- NULL
    write.csv(pval.res[-c(4,5,7)],file=paste(filenameRoot,"ERCC Pvals.csv"),
              row.names = FALSE)
    print("Finished DE testing")
    
    exDat$Results$quasiSeq.res <- quasiSeq.res
    exDat$Results$ERCCpvals <- pval.res
    
    #### Code from QuasiSeq
    #if (!Dispersion %in% c("Deviance", "Pearson")) 
    #  stop("Unidentified Dispersion: Dispersion must be either 'Deviance' 
    # or 'Pearson'.")
    fit <- use.fit
    spline.df <- NULL
    LRT <- fit$LRT
    phi.hat <- fit$phi.hat.dev
    
    mn.cnt <- fit$mn.cnt
    den.df <- fit$den.df
    num.df <- fit$num.df
    Model = fit$Model
    
    phi.hat<-use.fit$phi.hat.dev[1:nrow(cnt)]
    
    mn.cnt<-rowMeans(cnt)
    den.df=length(trt)-length(unique(trt))
    
    
    #if (Dispersion == "Pearson") 
    #  phi.hat <- fit$phi.hat.pearson
    if (length(num.df) == 1) 
        num.df <- rep(num.df, ncol(LRT))
    shrink.phi <- function(phi.hat, den.df) {
        phi.hat[phi.hat <= 0] <- min(phi.hat[phi.hat > 0])
        z <- log(phi.hat)
        z[z == Inf] <- max(z[z != Inf])
        z[z == -Inf] <- min(z[z != -Inf])
        mnz <- mean(z)
        d0arg <- var(z) - trigamma((den.df)/2)
        if (d0arg > 0) {
            dif <- function(x, y) abs(trigamma(x) - y)
            inverse.trigamma <- function(y) {
                optimize(dif, interval = c(0,10000),y = y)$minimum
            }
            d0 <- 2 * inverse.trigamma(d0arg)
            phi0 <- exp(mnz - digamma((den.df)/2) + digamma(d0/2) - 
                            log(d0/(den.df)))
            phi.shrink <- ((den.df) * phi.hat + d0 * phi0)/(den.df + 
                                                                d0)
        }
        else {
            phi.shrink <- rep(exp(mnz), length(z))
            d0 <- Inf
            phi0 <- exp(mnz)
        }
        return(list(phi.shrink = phi.shrink, d0 = d0, phi0 = phi0))
    }
    phi.hat[phi.hat < 0] <- min(phi.hat[phi.hat > 0])
    phi.hat2 <- phi.hat
    if (Model == "Poisson") 
        phi.hat2[phi.hat < 1] <- 1
    shrink <- shrink.phi(phi.hat, den.df)
    phi.shrink <- shrink[[1]]
    est.d0 <- shrink[[2]]
    if (Model == "Poisson") 
        phi.shrink[phi.shrink < 1] <- 1
    y <- log(phi.hat)
    y[y == -Inf] <- min(y[y != -Inf])
    y[y == Inf] <- max(y[y != Inf])
    spline.fit <- if (is.null(spline.df)) 
        smooth.spline(x = log(mn.cnt), y = y)
    else smooth.spline(x = log(mn.cnt), y = y, df = spline.df)
    spline.pred <- predict(spline.fit, x = log(mn.cnt))$y
    fit.method <- "spline"
    y2 <- phi.hat/exp(spline.pred)
    shrink <- shrink.phi(y2, den.df)
    D0 <- shrink[[2]]
    phi0 <- shrink[[3]]
    print(paste("Spline scaling factor:", phi0))
    
    mean = mn.cnt
    names(y)<-rownames(cnt)
    
    dispcnt = data.frame(mean = mean,y = y)
    
    xsort = NULL
    ysort = NULL
    dispcntSort = data.frame(xsort = sort(mean), 
                             ysort = spline.pred[order(mn.cnt)]) # original
    
    mean <- mean[ERCC]
    y <- y[ERCC]
    Ratio <- NULL
    dispERCC = data.frame(mean = mean[ERCC],y = y[ERCC], 
                          Ratio=
                              ERCC.Ratio$Ratio[match(ERCC,ERCC.Ratio$Feature)])
    
    dispERCC$Ratio <- as.factor(dispERCC$Ratio) 
    
    
    quasiDispPlot = ggplot() + geom_point(data = dispcnt, aes(x = mean, y = y),
                                          colour = "grey80", size = 5,
                                          alpha = 0.6) +
        geom_point(data = dispERCC, aes(x = mean, y = y,colour = Ratio), 
                   size = 5, alpha = 0.6) + xlab("Mean Counts") + 
        ylab("Log Dispersion Estimates from QuasiSeq)") + 
        stat_smooth(data = dispcntSort,aes(x = xsort, y = ysort),
                    colour = "black") + 
        scale_x_log10() + colScale + theme_bw() +
        theme(legend.justification=c(1,1), legend.position=c(1,1)) 
    
    exDat$Figures$dispPlot <- quasiDispPlot
    exDat$Results$simcnt <- simcnt 
    cat("\nFinished examining dispersions\n")
    
    return(exDat)
    ### end Edit Sarah Munro 20140216
}