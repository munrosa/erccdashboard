#' Save erccdashboard plots to a pdf file
#'
#' @param exDat         list, contains input data and stores analysis results
#' @param plotsPerPg    string, if "main" then the 4 main plots are printed  
#'                      to one page, if "single" then a single plot is printed 
#'                      per page from the plotlist argument
#' @param saveas    Choose file format from "pdf", "jpeg" or "png"
#' @param outName   Choose output file name, default will be fileName from exDat
#' @param res       Choose the file resolution                             
#' @param plotlist  list, contains plots to print
#' 
#' @description
#' The function savePlots will save selected figures to a pdf file. The default 
#' is the 4 manuscript figures to a single page (plotsPerPg = "manuscript"). 
#' If plotsPerPg = "single" then each plot is placed on an 
#' individual page. If plotlist is not defined (plotlist = NULL) or 
#' if plotlist = exDat$Figures then all plots in exDat$Figures are printed to a PDF file.
#' 
#' @examples
#' \donttest{
#' data(SEQC.Example)
#'  
#' exDat <- initDat(datType="count", isNorm=FALSE, exTable=MET.CTL.countDat, 
#'                  filenameRoot="testRun", sample1Name="MET",
#'                  sample2Name="CTL", erccmix="RatioPair", 
#'                  erccdilution=1/100, spikeVol=1, totalRNAmass=0.500,
#'                  choseFDR=0.1)
#'                  
#' exDat <- est_r_m(exDat)
#'                   
#' exDat <- dynRangePlot(exDat)
#'
#' exDat <- geneExprTest(exDat)
#' 
#' exDat <- erccROC(exDat)
#' 
#' exDat <- estLODR(exDat, kind="ERCC", prob=0.9)
#' 
#' exDat <- annotLODR(exDat)
#' 
#' #to print 4 main plots to a single page pdf file
#' saveERCCPlots(exDat, plotsPerPg = "manuscript",saveas = "pdf")
#' 
#' #to print 4 plots to a jpeg file
#' saveERCCPlots(exDat, plotsPerPg = "manuscript",saveas = "jpeg")
#' 
#' # or to create a multiple page pdf of all plots produced
#' saveERCCPlots(exDat, plotsPerPg = "single", plotlist = exDat$Figures)
#' 
#' # or to create a multiple page pdf of just 2 plots
#' saveERCCPlots(exDat, plotsPerPg = "single", 
#'               plotlist = list(exDat$Figures$lodrPlot, exDat$Figures$maPlot))
#' 
#' }
#' 
#' 
#' @export

saveERCCPlots<-function(exDat, plotsPerPg = "main", saveas = "pdf", outName,
                        plotlist, res){
    # Options are either the default of printing the plots as shown in 
    # the erccdashboard 2014 Nature Communications publication 
    # (plotsPerPg = "main" and plotlist is NULL) or 
    # to print plots one per page choose plotsPerPg = "single" and provide any
    # list of plots as the plotlist arguement
    
    if(missing(plotsPerPg)){
        plotsPerPg <- "main"
        cat(paste0("\nSaving main dashboard plots to ",saveas," file...\n"))
    }
    if(missing(saveas)){
        saveas <- "pdf"
        cat("Printing to PDF...")
    }
    if(missing(outName)){
        outName <- exDat$sampleInfo$filenameRoot
    }
    if(missing(res)) res = 200
    if (plotsPerPg == "main"){
        cols = 2
        nFigs = 4
        pwidth = 7*cols
        pheight = 7*nFigs/cols
        switch(saveas,
               png = png(filename = paste(outName,"png",sep="."),
                         width=pwidth,height = pheight,
                         units = "in", res = res),
               jpeg = jpeg(filename = paste(outName,"jpeg",sep="."),
                           width=pwidth,height = pheight, 
                           units = "in", res = res),
               pdf = pdf(file = paste(outName,"pdf",sep="."),
                         width=pwidth,height = pheight),
               stop(cat("\"saveas\" = ",saveas,
                        ", it must be \"pdf\", \"png\", or \"jpeg\"")))
        multiplot(exDat$Figures$dynRangePlot, exDat$Figures$rocPlot,
                  exDat$Figures$maPlot, exDat$Figures$lodrERCCPlot, cols=cols)
        dev.off()
    }
    if (plotsPerPg == "single"){
        pwidth = 7
        pheight = 7
        if (missing(plotlist)){
            plotlist = exDat$Figures
        }
        if(is.list(plotlist)){           
            nPlotsSum <- dim(summary(plotlist))[1]
            nPlotsLength <- length(plotlist)
            if(nPlotsSum == nPlotsLength){
                if(saveas != "pdf") cat("Ignoring \"saveas\" = \"",saveas,"\"\n")
                cat("Printing plots to multiple pages in one PDF file.\n")
                pdf(file=paste(outName,"pdf",sep="."),
                    onefile=TRUE,width=pwidth,
                    height=pheight)
                print(plotlist)
                dev.off()
            }
        }else{
            stop(cat("\"plotlist\" is not a list! ")) 
        }    
    } 
}