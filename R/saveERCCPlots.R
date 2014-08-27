#' Save erccdashboard plots to a pdf file
#'
#' @param exDat     list, contains input data and stores analysis results
#' @param plotsPerPg string, if "manuscript" then the 4 main plots are printed  
#'                   to one page, if "single" then each plot is printed to page 
#'                   in the pdf file
#' @param saveas     Choose file format from "pdf", "jpeg" or "png"
#' @param res        Choose the file resolution                             
#' @param plotlist   list, contains plots to print
#' 
#' @description
#' The function savePlots will save selected figures to a pdf file. The default 
#' is the 4 manuscript figures to a single page (plotsPerPg = "manuscript"). 
#' If plotsPerPg = "single" then each plot is placed on an 
#' individual page in one pdf file. If plotlist is not defined (plotlist = NULL)
#'  then all plots in exDat$Figures are printed to the file.
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
#' exDat <- estLODR(exDat, kind="ERCC", prob=0.9)
#' 
#' exDat <- annotLODR(exDat)
#' 
#' #to print 4 plots from manuscript to a single page pdf file
#' saveERCCPlots(exDat, plotsPerPg = "manuscript")
#' 
#' # or to create a multiple page pdf of all plots produced
#' saveERCCPlots(exDat, plotsPerPg = "single", plotlist = exDat$Figures)
#' }
#' 
#' 
#' @export

saveERCCPlots<-function(exDat,plotsPerPg = "manuscript", saveas = "pdf", 
                        res = 180, plotlist = NULL){
    # Options are either the default of printing the plots as shown in 
    ## publication (plotsPerPg = "manuscript" and plotlist is NULL) or 
    # to print plots one per page choose (plotsPerPg = "single" and provide any
    # list of plots as the plotlist arguement
    
    
    # Open PDF file to write results
    outName <- exDat$sampleInfo$filenameRoot 
    #   if (plotsPerPg == "manuscript"){
    #     cols = 2
    #     pwidth = 7*cols
    #     pheight = 7*6/cols
    #     pdf(filename = paste(outName,"pdf",sep="."), 
    #         width=pwidth,height = pheight)
    #     
    #     multiplot(exDat$Figures$rocPlot,exDat$Figures$dynRangePlot, 
    #               exDat$Figures$lodrERCCPlot,exDat$Figures$rangeResidPlot, 
    #               exDat$Figures$dispPlot,exDat$Figures$maPlot,cols=2)
    #     dev.off()
    #   } 
    
    cat(paste0("\nSaving main dashboard plots to ",saveas," file...\n"))
    if (plotsPerPg == "manuscript"){
        cols = 2
        nFigs = 4
        pwidth = 7*cols
        pheight = 7*nFigs/cols
        if (saveas == "png"){
            png(filename = paste(outName,"png",sep="."),
                width=pwidth,height = pheight, units = "in", res = res)
            multiplot(exDat$Figures$dynRangePlot, exDat$Figures$rocPlot,
                      exDat$Figures$maPlot, exDat$Figures$lodrERCCPlot, cols=cols)
            dev.off()
        }
        if (saveas == "pdf"){
            pdf(file = paste(outName,"pdf",sep="."),
                width=pwidth,height = pheight)
            multiplot(exDat$Figures$dynRangePlot, exDat$Figures$rocPlot,
                      exDat$Figures$maPlot, exDat$Figures$lodrERCCPlot, cols=cols)
            dev.off()
        }
        if (saveas == "jpeg"){
            jpeg(filename = paste(outName,"jpeg",sep="."),
                width=pwidth,height = pheight, units = "in", res = res)
            multiplot(exDat$Figures$dynRangePlot, exDat$Figures$rocPlot,
                      exDat$Figures$maPlot, exDat$Figures$lodrERCCPlot, cols=cols)
            dev.off()
        }
        else{
            stop(cat("Choose 'pdf', 'png', or 'jpeg' file formats"))
        }
    }
    if (plotsPerPg == "single"){
        if (is.null(plotlist)){
            plotlist = exDat$Figures
        } 
        if (saveas == "png"){
            png(filename = paste(outName,"png",sep="."),width=7,
                height = 7, units = "in", res =res)
            print(plotlist)
            dev.off()
        }
        if (saveas == "pdf"){
            pdf(file = paste(outName,"pdf",sep="."),width=7,
                height = 7)
            print(plotlist)
            dev.off()    
        }
        if (saveas == "jpeg"){
            jpeg(filename = paste(outName,"jpeg",sep="."),width=7,
                height = 7, units = "in", res =res)
            print(plotlist)
            dev.off()
        }
        else{
            stop(cat("Choose 'pdf', 'png', or 'jpeg' file formats"))
        }
        
    } 
}
