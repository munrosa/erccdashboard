getDesignMat <- function(expressionData, 
                         factorList=c("Study", "Platform", "Site", "Sample",
                                      "Library", "Lane", "Barcode", "FlowCell",
                                      "AlignMethod"), patternSplit = '_'){
  designMatrix <- data.frame(countSet=
                                 as.character(names(expressionData)[-c(1)]))
  factorSet <- colsplit(designMatrix$countSet,pattern=patternSplit,
                       names=factorList)
  factorSet <- data.frame(factorSet, stringsAsFactors= TRUE)
  designMatrix <- cbind(designMatrix, factorSet)
  designMatrix <- as.data.frame(lapply(designMatrix,as.factor))
  
  return(designMatrix)
}
