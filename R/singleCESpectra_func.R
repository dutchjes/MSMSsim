## 

#' spectra from single CE function
#'
#' @param all.frag 
#' @param energy collision energy desired
#' @param dmz 
#'
#' @return
#' @export
#'
#' @examples
#' 
singleCESpectra <- function(all.frag, energy = "15", dmz = 0.001){
  
  frags <- subset(all.frag, all.frag$CollisionEnergy == energy)
  IDs <- as.numeric(unique(frags$CompoundID))
  IDs <- IDs[order(IDs, decreasing = FALSE)]
  singleCE.spec <- list();
  
  for(i in 1:length(IDs)){
    
    where <- which(frags$CompoundID == IDs[i])
    singleCE.spec[[i]] <- as.data.frame(frags[where,c("FragmentMZ", "FragmentIntensity")])
    colnames(singleCE.spec[[i]]) <- c("mz", "int")

    singleCE.spec[[i]] <- filterMerge(singleCE.spec[[i]], dmz = dmz)
  }
  
  names(singleCE.spec) <- IDs
  return(singleCE.spec)
}