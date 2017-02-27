

#' Function to get the Similarity Score Threshold for desired FPR
#'
#' @param fptp output from the PandRdat function. A list with each of the scenarios to evaluate. Will also accept a data frame 
#' @param minFPR the desired minimum false positive rate
#'
#' @return a vector equal to the length of fptp input, with the minimum similarity score needed for each scenario to reach the desired false positive rate
#' @export
#'
#' @examples
#' 
getCutoff <- function(fptp, minFPR = 0){
  
  if(is.data.frame(fptp)==TRUE){
    fptpdat <- fptp
    fptpdat[,2:5] <- fptpdat[,2:5] * 100
    fptpdat$TP <- fptpdat$`True pos`
    fptpdat$FN <- fptpdat$TP + fptpdat$`False neg`
    
    fptpdat$FP <- fptpdat$`False pos`
    fptpdat$TN <- fptpdat$FP + fptpdat$`True neg`
    cutoff <- min(which(fptpdat$`False pos` <= minFPR))
    result <- fptpdat[cutoff, "Threshold"]
  }else{
    
    result <- rep(NA, length(fptp))
    
    for(i in 1:length(fptp)){
      
      
      fptpdat <- fptp[[i]]
      fptpdat[,2:5] <- fptpdat[,2:5] * 100
      fptpdat$TP <- fptpdat$`True pos`
      fptpdat$FN <- fptpdat$TP + fptpdat$`False neg`
      
      fptpdat$FP <- fptpdat$`False pos`
      fptpdat$TN <- fptpdat$FP + fptpdat$`True neg`
      cutoff <- min(which(fptpdat$`False pos` <= minFPR))
      result[i] <- fptpdat[cutoff, "Threshold"]
    }
    
  }
  
  return(result)
  
}
