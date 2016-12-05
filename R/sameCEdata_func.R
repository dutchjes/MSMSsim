
### 

#' Producing the list for comparing MS/MS spectra from same CEs for the pairs. Largely replaced by singleCESpectra_func
#'
#' @param all.frag 
#' @param pairs 
#'
#' @return Returns list of three lists. First is the summary of each collision energy, second is the data of parents at all 
#' CEs, third is the data of TPs at all CEs.
#' @export
#'
#' @examples
#' 
sameCEdata <- function(all.frag, pairs){


TP.ID <- as.integer(pairs$TP_ID)
par.ID <- as.integer(pairs$par_ID)
pair.ID <- as.integer(pairs$pair_ID)
y <- list()
z <- list()
nowp <- list(0)
nowt <- list(0)

nCE <- dim(table(all.frag$CollisionEnergy))
sameCE <- as.list(rep("", nCE))
names(sameCE) <- unique(all.frag$CollisionEnergy)


for(j in 1:nCE){
  
  opCE <- names(table(all.frag$CollisionEnergy))[j]
  sameCE[[j]] <- data.frame(matrix(NA, ncol = 1, nrow = nrow(pairs)))
  
  #colnames(sameCE) <- names(table(all.frag$CollisionEnergy))
  rownames(sameCE[[j]]) <- rownames(pairs)
  
  sameCE[[j]]$result <- c()
  sameCE[[j]]$parMZ <- c()
  sameCE[[j]]$parIon <- c()
  sameCE[[j]]$parOrigin <- c()
  sameCE[[j]]$parRes <- c()
  sameCE[[j]]$parCE <- c()
  sameCE[[j]]$TPMZ <- c()
  sameCE[[j]]$TPIon <- c()
  sameCE[[j]]$TPOrigin <- c()
  sameCE[[j]]$TPRes <- c()
  sameCE[[j]]$TPCE <- c()
  sameCE[[j]]$pairID <- pair.ID
  
  for(i in 1:length(TP.ID)){
    
    
    tpid <- as.integer(paste(TP.ID[i]))
    wheret <- which(all.frag$CompoundID == tpid)
    if(length(wheret) == 0){
      y[[i]] <- as.data.frame(c(0,0))
      z[[i]] <- as.data.frame(c(0,0))
      sameCE[[j]][i,"result"] <- paste("No TP")
      next
    }
    
    tpmass <- as.numeric(paste(all.frag[wheret[1], "PrecursorMZ"]))
    sameCE[[j]][i, "TPMZ"] <- tpmass
    
    ptid <- as.integer(paste(par.ID[i]))
    wherep <- which(all.frag$CompoundID == ptid)
    if(length(wherep) == 0){
      y[[i]] <- as.data.frame(c(0,0))
      z[[i]] <- as.data.frame(c(0,0))
      sameCE[[j]][i,"result"] <- paste("No par")
      next
    }
    
    ptmass <- as.numeric(paste(all.frag[wherep[1], "PrecursorMZ"]))
    sameCE[[j]][i, "parMZ"] <- ptmass
    
    tpion <- paste(all.frag[wheret[1], "IonMode"])
    sameCE[[j]][i,"parIon"] <- tpion
    ption <- paste(all.frag[wherep[1], "IonMode"])
    sameCE[[j]][i,"TPIon"] <- ption
    
    if(ption != tpion){
      sameCE[[j]][i,"result"] <- paste("Different ionization modes!!")
      y[[i]] <- as.data.frame(c(0,0))
      z[[i]] <- as.data.frame(c(0,0))
      next
    }
    
    
    ### Extraction of TP fragments
    
    tpdat <- subset(all.frag, all.frag$CompoundID == tpid)
    tpdat <- subset(tpdat, tpdat$CollisionEnergy == opCE)
    
    if(nrow(tpdat) == 0){
      y[[i]] <- as.data.frame(c(0,0))
      z[[i]] <- as.data.frame(c(0,0))
      sameCE[[j]][i, "result"] <- paste("No TP at this CE")
      next
    }
    
    ptdat <- subset(all.frag, all.frag$CompoundID == ptid)
    ptdat <- subset(ptdat, ptdat$CollisionEnergy == opCE)
    ptdat <- subset(ptdat, ptdat$IonMode == tpion)
    
    if(nrow(ptdat) == 0){
      y[[i]] <- as.data.frame(c(0,0))
      z[[i]] <- as.data.frame(c(0,0))
      sameCE[[j]][i, "result"] <- paste("No parent at this CE")
      next
    }
    
    
    res <- unique(tpdat$Resolution)
    
    if(length(res) > 1){
      
      take <- max(res)
      tpdat <- subset(tpdat, tpdat$Resolution == take)
      sameCE[[j]][i,"TPRes"] <- take
    }else{
      sameCE[[j]][i,"TPRes"] <- res
      
    }
    
    if(mean(table(tpdat$FragmentRank)) != 1 ){
      
      take <- unique(tpdat$origin)[1]
      tpdat <- subset(tpdat, tpdat$origin == take)
      sameCE[[j]][i,"TPOrigin"] <- take
      # break
      
    }else{
      sameCE[[j]][i,"TPOrigin"] <- unique(tpdat$origin)
    }
    
    tp <- as.data.frame(cbind(tpdat$FragmentMZ, tpdat$FragmentIntensity))
    y[[i]] <- as.data.frame(tp[order(tp$V1), ])
    y[[i]] <- unique(y[[i]])
    
    
    ### Extraction of parent fragments
    
    
    res <- unique(ptdat$Resolution)
    
    if(length(res) > 1){
      
      take <- max(res)
      ptdat <- subset(ptdat, ptdat$Resolution == take)
      sameCE[[j]][i,"parRes"] <- take
    }else{
      sameCE[[j]][i,"parRes"] <- res
    }
    
    if(mean(table(ptdat$FragmentRank)) != 1 ){
      
      take <- unique(ptdat$origin)[1]
      ptdat <- subset(ptdat, ptdat$origin == take)
      sameCE[[j]][i,"parOrigin"] <- take
      # break
      
    }else{
      sameCE[[j]][i,"parOrigin"] <- unique(ptdat$origin)
    }
    
    pt <- as.data.frame(cbind(ptdat$FragmentMZ, ptdat$FragmentIntensity))
    z[[i]] <- as.data.frame(pt[order(pt$V1), ])
    z[[i]] <- unique(z[[i]])
    
    sameCE[[j]][i,"result"] <- paste("both measured")
    
  }
  
  nowp[[j]] <- z
  nowt[[j]] <- y

}

rm(pt, ptdat, tp, tpdat, opCE, ptid, ption, ptmass, res, take, tpid, tpion, tpmass, wherep, wheret)

return(list(sameCE = sameCE, nowt = nowt, nowp = nowp))

}