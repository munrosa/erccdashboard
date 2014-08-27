loadExpMeas<- function(exDat, exTable, repNormFactor){
  
  designFactors <- c("Sample","Rep")
  datType <- exDat$sampleInfo$datType
  isNorm <- exDat$sampleInfo$isNorm
  # Check for exTable input errors
  for (i in 2:length(colnames(exTable))){
      if (str_count(colnames(exTable[c(i)]),"_") != 1){
        stop("Check exTable column names for use of underscore (_)")
      }
  }
  if (anyDuplicated(colnames(exTable)) > 0){
    stop("Column names must be unique sample IDs")  
  }
  sampleInfo = exDat$sampleInfo
  # Pull idCols out of exDat
  idCols <- exDat$erccInfo$idCols 

  if (missing(repNormFactor)){
    repNormFactor <- NULL
  }
  #}else{
    exDat$sampleInfo$repNormFactor = repNormFactor
  #}
  
  Transcripts = exTable
  
  # Import data based on analysis type from SEQC main project
    
   # force names to be ERCC- and first column name to Feature
   names(Transcripts)[1] = "Feature"
   
   #Transcripts$Feature = gsub(".","-",Transcripts$Feature)
   Transcripts$Feature = gsub(":","",Transcripts$Feature)
   row.names(Transcripts) <- gsub("[[:punct:]]", "_", row.names(Transcripts))
   Transcripts$Feature <- gsub("[[:punct:]]", "_", Transcripts$Feature)
  Transcripts$Feature = gsub("ERCC_","ERCC-",Transcripts$Feature)
  #print(head(Transcripts))
   # get data frames with just the ERCCs and just the human genes
   TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
   TranscriptsHumanOnly = Transcripts[-c(grep("ERCC-0", Transcripts$Feature)),]
   
   # Remove ERCCs in the definition file that are not in the count data file
   idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]
   
   # Remove ERCCs without a Ratio
   idCols = idCols[which(is.finite(idCols$Ratio)),]
   
   # Remove ERCCs from count data and idCols that are absent from the experiment
   TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,
                                                  TranscriptsERCCOnly$Feature),]
   Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)
   
   #############################################################################
   sample1 <- sampleInfo$sample1Name 
   sample2 <- sampleInfo$sample2Name

  designMat <- getDesignMat(expressionData = Transcripts,
                               factorList = designFactors,
                               patternSplit = '_')
  
  ### Filter the transcripts
  if ((datType == "count")&(isNorm == FALSE)){
    lengthinit <- dim(Transcripts)[1]
    idxsample <- which((rowMeans(Transcripts[-c(1)])>1)&(rowSums(
      Transcripts[-c(1)]!=0)>2))
    
    Transcripts <- Transcripts[idxsample,]
    
    Transcripts$Feature <- as.factor(as.character(Transcripts$Feature))
    
    measERCCs <- Transcripts$Feature[grep("ERCC-0", Transcripts$Feature)]
    
    insuffDat <- setdiff(idCols$Feature, measERCCs)
    
    cat(paste("\nTranscripts were removed with a mean count < 1 or more than 2",
              "\nreplicates with 0 counts."))
    cat(paste("\nOriginal data contained ",lengthinit,"transcripts.",
              "\nAfter filtering ",length(idxsample),"transcripts remain for ",
              "analysis."))
    cat(paste("\nA total of",length(insuffDat),"out of",length(idCols$Feature),
              "\nERCC controls were filtered from the data set"))
    if(length(insuffDat > 0)){
      cat("\nThe excluded ERCCs are:\n")
      for (j in seq(from=1,to=length(insuffDat),by=5)){
        k = j+4
        if (k > length(insuffDat)) k = length(insuffDat)
        cat(insuffDat[j:k])
        cat("\n")
      }  
    }
     
  }
  
  
  # write Transcript csv file to directory
  #write.csv(Transcripts, paste(sampleInfo$filenameRoot,"Transcripts.csv",
  ##        sep="."), row.names = F)
  # collect everything to add to exDat
  exDat = append(exDat, list(Transcripts = Transcripts,
                               designMat = designMat,
                               sampleNames = c(sample1,sample2),
                               idCols = idCols))
  return(exDat)

}