
#' Analysis of Single Collion Energies
#'
#' @param pairs Data frame containing the pairs which are to be compared. First column should be "pair_ID", second column "TP_ID",
#'  and the third column "par_ID". 
#' @param nowp List which contains all the MSMS spectra for the parent compounds. Output from sameCEdata function
#' @param nowt List which contains all the MSMS spectra for the TP compounds. Output from sameCEdata function
#' @param sameCE List which contains the measurement information for each of the collision energies. Output from sameCEdata function
#' @param wMD Should the spectra of the TP be shifted by the mass difference of the pair? Default is FALSE
#' @param y m/z weighting factor for simTableWeighting(). Default is 0
#' @param z intensity weighting factor for simTableWeighting(). Default is 1
#' @param t m/z tolerance for merging factors for simTable(). Default is 0.005
#' @param b relative intensity cutoff for relCutoff(). Default is 0.05.
#' @param replace.na value to replace NAs with. Default is 0

#' @param sink Should the data be output into a .txt file? Default in FALSE
#'
#' @return List with the first entry the summary data frame with scores of all pairs at each collision energy. Second entry is a list, with
#' an entry for each of the collision energies analyzed with the measurement information.
#' 
#' @export
#'
#' @examples

sameCEanalysis <- function(pairs, nowp, nowt, sameCE, wMD = FALSE ,y=1, z=1, t=0.005, b=0.5, replace.na = 0, sink = FALSE){
  
summarySim <- data.frame(c(pairs$pair_ID))
par.ID <- pairs$par_ID
TP.ID <- pairs$TP_ID
nCE <- length(names(sameCE))
engs <- names(sameCE)

for(j in 1:nCE){
  
  opCE <- engs[j]
  goodpairs_sameCE <- which(sameCE[[j]][,"result"] == "both measured")
  sameCE[[j]][,"sim"] <- c()
  name <- paste("SimilarityResults_CE", opCE, ".txt", sep = "")
  sameCE[[j]][,"MD"] <- sameCE[[j]][,"TPMZ"] - sameCE[[j]][,"parMZ"]
  
  #exist <- c(1:nrow(sameCE))
  #good.pairs_sameCE <- setdiff(exist, good.pairs_sameCE)
  
  
  if(sink == TRUE){
    sink(name)
  }
  
  for(i in 1:length(goodpairs_sameCE)){
    
    ivalue <- as.integer(paste(goodpairs_sameCE[i]))
    
    par.name <- par.ID[ivalue]
    TP.name <- TP.ID[ivalue]
    
    if(sink==TRUE){
      print(paste(par.name, TP.name, "Pair", sameCE[[j]][ivalue, "pairID"], sep = " "))
    }
    
    parMSMS <- nowp[[j]][[ivalue]]
    tpMSMS <- nowt[[j]][[ivalue]]
    
    if(wMD == TRUE){
      MD <- as.numeric(sameCE[[j]][ivalue,"MD"])
      
      tpMSMS$V1 <- tpMSMS$V1 - MD
      
    }
    
    MSMSsim <- OrgMSSim(spec.top = parMSMS, spec.bottom = tpMSMS, y=y, z=z, t = t, b = b, replace.na = replace.na
                      )
    sameCE[[j]][ivalue,"sim"] <- MSMSsim$score
 
    if(sink==TRUE){
      print(MSMSsim$result)
    }
  }
  
  if(sink==TRUE){
    sink()
  }
  
  summarySim[,j+1] <- sameCE[[j]][,"sim"]
  
  
}

colnames(summarySim) <- c("pair_ID", engs)
return(list(summarySim = summarySim, sameCE = sameCE))

}
