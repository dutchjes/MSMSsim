
#' PandR function for looking at FPR vs. TPR.
#'
#' @param realdata 
#' @param randdata 
#' @param steps 
#' @param metricList 
#' @param specScenario 
#' @param PandRplot 
#' @param legend.lab 
#' @param save.data 
#' @param auc 
#'
#' @return Large wrapper for all PandR functions. Has largely been replaced by ROC_func
#' @export
#'
#' @examples
PandR <- function(realdata, randdata, steps = 0.001, metricList = metricList,
                  specScenario = "AllFragments", PandRplot = TRUE, legend.lab = metricList, save.data = TRUE, auc = TRUE){
  
  #library(pracma)
  
  #   if(metricList %in% colnames(realdata) == FALSE){
  #     print("Metrics not in realdata!")
  #     stop
  #       }
  # 
  #   if(metricList %in% colnames(randdata) == FALSE){
  #     print("Metrics not in realdata!")
  #     stop
  #   }
  
  
  # pandr <- data.frame(matrix(data = NA, nrow = 1/steps, ncol = 3))
  #  colnames(pandr) <- c("Threshold", "False pos", "True pos")
  
  if(save.data == TRUE){
    PandRList <- list()
  }
  
  if(auc == TRUE){
    aucs <- data.frame(matrix(data = NA, ncol = 2, nrow = length(metricList)))
  }
  
  for(k in 1:length(metricList)){
    
    pandr <- data.frame(matrix(data = NA, nrow = 1/steps, ncol = 3))
    colnames(pandr) <- c("Threshold", "False pos", "True pos")
    
    metric <- metricList[k]
    real.npair <- length(na.omit(realdata[,paste(metric)]))
    rand.npair <- length(na.omit(realdata[,paste(metric)]))
    
    
    for(i in 1:(1/steps+1)){
      
      thres <- (i-1)*steps
      pandr[i,1] <- thres
      pandr[i,2] <- length(which(randdata[,paste(metric)] > thres)) / real.npair
      pandr[i,3] <- length(which(realdata[,paste(metric)] > thres)) / rand.npair
      
    }
    
    if(save.data == TRUE){
      pandr$metric <- c(rep(metric, 1/steps+1))
      PandRList[[k]] <- pandr
      #   name(PandRList[k]) <- metric
      #   return(PandRList)
    }
    
    if(auc == TRUE){
      aucs[k,1] <- metric
      aucs[k,2] <- as.numeric(trapz(pandr[,3], pandr[,2]))
    }
    
    if(PandRplot == TRUE){
      
      metcol <- rainbow(length(metricList))
      
      if(k == 1){
        png(filename = paste(specScenario, "_PandR.png", sep = ""))
        plot.new()
        box()
        title(main = paste(specScenario, "ROC Plot", sep = " "))
        axis(side = 1)
        axis(side = 2)
        mtext(side = 1, "False positive rate (FPR)", line = 2.5)
        mtext(side = 2, "True positive rate (TPR)", line = 2.5)
        
        
        points(pandr[,2:3], pch = 19, col = metcol[k])
      }else{
        points(pandr[,2:3], pch = 19, col = metcol[k])
        
      }
      
      if(k == length(metricList)){
        # abline(v = 0.05, lty = 2)
        #  abline(h = 0.80, lty = 2)
        grid(col = "darkgray")
        legend("bottomright", legend = legend.lab, col = metcol, pch = 19, cex = 1.2, #title = "Combinations",
               ncol = ifelse(length(metricList)>8,2,1))
        
        dev.off()
      }
      
    }
    
  }
  
  if(save.data == TRUE){
    output <- do.call(rbind, PandRList)
    write.csv(output, file = paste("PandRList_", specScenario, ".csv", sep = ""))
    print(paste("Saved PandRList to ", getwd(), sep = ""))
  }
  
  if(auc == TRUE){
    write.csv(aucs, file = paste("AUCs_", specScenario, ".csv", sep = ""))
    print(paste("Saved AUC table to ", getwd(), sep = ""))
  }
  
  # return(PandRList)
  # return(aucs)
}
