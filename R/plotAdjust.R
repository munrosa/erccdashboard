plotAdjust <- function(expDat){
  # this function sets the general plot style parameters
  theme_set(theme_bw(base_size=12))
  sampleInfo <- expDat$sampleInfo
  erccInfo <- expDat$erccInfo
  
  if((sampleInfo$erccmix =="RatioPair")|
      (sampleInfo$erccmix =="Single")){
    #Create a custom color scale
    # yellow, green, blue, red
    myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
    
    names(myColors) <- erccInfo$FCcode$Ratio
    colScale <- scale_colour_manual(name = "Ratio",values = myColors)
     #                               labels = sampleInfo$legendLabels)
    fillScale <- scale_fill_manual(name = "Ratio", values = myColors) 
      #                             labels = sampleInfo$legendLabels)
  }else{
    colScale <- scale_colour_hue(name = "Ratio")
    fillScale <- scale_fill_hue(name = "Ratio")
  }
    
  expDat$plotInfo$colScale <- colScale
  expDat$plotInfo$fillScale <- fillScale
  return(expDat)
  
}