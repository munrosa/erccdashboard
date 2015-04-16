#' Multiplot function from R cookbook
#'
#' @param ...       comma separated list of plots to include
#' @param plotlist  list, naming plots to include
#' @param cols      number of columns for grid of plots
#' 
#' @description
#' # Function is originally from R Cookbook
#' # http://wiki.stdout.org/rcookbook/Graphs/Multiple%20
#' #          graphs%20on%20one%20page%20(ggplot2)/
#' 
#' @examples
#' data(SEQC.Example)
#' 
#' exDat <- runDashboard(datType="array", isNorm=FALSE, 
#'                  exTable=UHRR.HBRR.arrayDat,
#'                  filenameRoot="testRun", sample1Name="UHRR",
#'                  sample2Name="HBRR", erccmix="RatioPair", 
#'                  erccdilution = 1, spikeVol = 50, 
#'                  totalRNAmass = 2.5*10^(3), choseFDR=0.01)
#'                 
#' # print 2 plots to page with 2 columns                
#' multiplot(exDat$Figures$lodrERCCPlot,exDat$Figures$rocPlot,cols = 2)
#'
#' # print 2 plots to page with 1 column
#' multiplot(exDat$Figures$lodrERCCPlot,exDat$Figures$rocPlot,cols = 1)
#' 
#' @export
#' 
multiplot <- function(..., plotlist=NULL, cols) {
    #require(grid)
    # function from R cookbook: 
    ##http://wiki.stdout.org/rcookbook/Graphs/Multiple%20
    ##          graphs%20on%20one%20page%20(ggplot2)/
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # Make the panel
    plotCols = cols # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed
    
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
    
}