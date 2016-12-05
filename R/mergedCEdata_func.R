#' Merged CE data for comparison
#'
#' @param all.frag 
#' @param pairs 
#' @param dmz 
#' @param intensity.type 
#'
#' @return list of two lists, first with merged spectra of parent, second with merged spectra of TP. New function mergedSpectra 
#' is better because doesn't require input of the pairs.
#' @export
#'
#' @examples
mergedCEdata <- function(all.frag, pairs, dmz = 0.001, intensity.type){
  
  if(intensity.type != "absolute" & intensity.type != "relative"){
    stop("Which type of intensity: absolute or relative?")
  }
  
  TP.ID <- as.integer(pairs$TP_ID)
  par.ID <- as.integer(pairs$par_ID)
  pair.ID <- as.integer(pairs$pair_ID)
  
  
  sdt <- list()
  sdp <- list()
  sdtCE <- list()
  sdpCE <- list()
  
  for(i in 1:length(TP.ID)){
    
    tpid <- as.integer(paste(TP.ID[i]))
    where1 <- which(all.frag$CompoundID == tpid)
    
    ptid <- as.integer(paste(par.ID[i]))
    where2 <- which(all.frag$CompoundID == ptid)
    
    if(length(where1) == 0){
      
      if(length(where2) == 0){
        #   pairs[i,"result"] <- paste("neither measured")
        
        sdt[[i]] <- as.data.frame(c(0,0))
        sdp[[i]] <- as.data.frame(c(0,0))
        sdtCE[[i]] <- c(0)
        sdpCE[[i]] <- c(0)
        next
      }
      
      sdt[[i]] <- as.data.frame(c(0,0))
      sdtCE[[i]] <- c(0)
      
      if(intensity.type == "absolute"){
          find <- "FragmentIntensity"
      }
      
      if(intensity.type == "relative"){
        find <- "FragmentRelativeIntensity"
      }
      
      sdp[[i]] <- as.data.frame(all.frag[where2,c("FragmentMZ", paste(find), "CollisionEnergy")])
      sdpCE[[i]] <- as.integer(unique(sdp[[i]][,"CollisionEnergy"]))
      sdp[[i]] <- sdp[[i]][,c("FragmentMZ", paste(find))]
      colnames(sdp[[i]]) <- c("mz", "int")
      # sdp[[i]][,"int"] <- sdp[[i]][,"int"]/max(sdp[[i]][,"int"])
      
      sdp[[i]] <- filterMerge(sdp[[i]], dmz = dmz)
      next
      
    }
    
    if(length(where2) == 0){
      sdt[[i]] <- as.data.frame(all.frag[where1,c("FragmentMZ", paste(find), "CollisionEnergy")])
      sdtCE[[i]] <- as.integer(unique(sdt[[i]][,"CollisionEnergy"]))
      sdt[[i]] <- sdt[[i]][,c("FragmentMZ", paste(find))]
      colnames(sdt[[i]]) <- c("mz", "int")
      sdt[[i]] <- filterMerge(sdt[[i]], dmz = dmz)
      
      sdp[[i]] <- as.data.frame(c(0,0))
      sdpCE[[i]] <- c(0)
      next
    }
    
    ## making merged table for the TP
    sdt[[i]] <- as.data.frame(all.frag[where1,c("FragmentMZ", paste(find), "CollisionEnergy")])
    sdtCE[[i]] <- as.integer(unique(sdt[[i]][,"CollisionEnergy"]))
    sdt[[i]] <- sdt[[i]][,c("FragmentMZ", paste(find))]
    colnames(sdt[[i]]) <- c("mz", "int")
    sdt[[i]] <- filterMerge(sdt[[i]], dmz = dmz)
    
    ## making the merged table for the parent
    sdp[[i]] <- as.data.frame(all.frag[where2,c("FragmentMZ", paste(find), "CollisionEnergy")])
    sdpCE[[i]] <- as.integer(unique(sdp[[i]][,"CollisionEnergy"]))
    sdp[[i]] <- sdp[[i]][,c("FragmentMZ", paste(find))]
    colnames(sdp[[i]]) <- c("mz", "int")
    sdp[[i]] <- filterMerge(sdp[[i]], dmz = dmz)
    
    
  }
  
  return(list(sdp = sdp, sdt = sdt))
  
}
