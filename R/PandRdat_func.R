
#' Generation of PandR data for scenario evaluation
#'
#' @param realdata 
#' @param randdata 
#' @param metricList 
#' @param steps 
#' @param specScenario 
#' @param save.data 
#'
#' @return
#' @export
#'
#' @examples
PandRdat <- function(realdata, randdata, metricList = metricList, steps = 0.001,
                     specScenario = "AllFragments", save.data = TRUE){
  
  PandRList <- list()
  
  if(is.null(metricList)){
    
    pandr <- data.frame(matrix(data = NA, nrow = 1/steps, ncol = 3))
    colnames(pandr) <- c("Threshold", "False pos", "True pos")
    
  #  metric <- metricList[k]
    real.npair <- length(realdata)
    rand.npair <- length(randdata)
    
    
    for(i in 1:(1/steps+1)){
      
      thres <- (i-1)*steps
      pandr[i,1] <- thres
      pandr[i,2] <- length(which(randdata > thres)) / rand.npair
      pandr[i,3] <- length(which(realdata > thres)) / real.npair
      
    }
    
    pandr[,4] <- 1 - pandr[,"True pos"]
    pandr[,5] <- 1 - pandr[,"False pos"]
    colnames(pandr) <- c("Threshold", "False pos", "True pos", "False neg", "True neg")
    
    
    #if(save.data == TRUE){
   # pandr$metric <- c(rep(metric, 1/steps+1))
    PandRList[[1]] <- pandr
    
  }else{
    
    for(k in 1:length(metricList)){
      
      pandr <- data.frame(matrix(data = NA, nrow = 1/steps, ncol = 3))
      colnames(pandr) <- c("Threshold", "False pos", "True pos")
      
      metric <- metricList[k]
      real.npair <- length(realdata[,paste(metric)])
      rand.npair <- length(randdata[,paste(metric)])
      
      
      for(i in 1:(1/steps+1)){
        
        thres <- (i-1)*steps
        pandr[i,1] <- thres
        pandr[i,2] <- length(which(randdata[,paste(metric)] > thres)) / rand.npair
        pandr[i,3] <- length(which(realdata[,paste(metric)] > thres)) / real.npair
        
      }
      
      pandr[,4] <- 1 - pandr[,"True pos"]
      pandr[,5] <- 1 - pandr[,"False pos"]
      colnames(pandr) <- c("Threshold", "False pos", "True pos", "False neg", "True neg")
      
      
      #if(save.data == TRUE){
      pandr$metric <- c(rep(metric, 1/steps+1))
      PandRList[[k]] <- pandr
      #   name(PandRList[k]) <- metric
      #   return(PandRList)
      #  }
      
    }
  }
  
  return(PandRList)
  
}

