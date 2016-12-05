
## Plotting of breakdown curves for each compound
## Need to load a all.frag table generated from MassBank

#' Title
#'
#' @param all.frag 
#' @param plot default is TRUE
#' @param start passed to rainbow() function
#' @param end passed to rainbow() function
#'
#' @return Just a graphing function
#' @export
#'
#' @examples
breakdown <- function(all.frag, plot = TRUE, start = 0, end = 1){
  

for(i in 1:length(unique(all.frag$CompoundID))){
  
  take <- unique(all.frag$CompoundID)[i]
  dat <- subset(all.frag, all.frag$CompoundID == take);
  dat <- subset(dat, dat$FragmentationMode == "HCD")
  mode <- as.data.frame(table(dat$IonMode));
  ionmode <- mode[which(mode$Freq == max(mode$Freq)),1];
  datnow <- subset(dat, dat$IonMode == ionmode);
  
  res <- as.data.frame(table(dat$Resolution));
  resol <- res[which(res$Freq == max(res$Freq)),1];
  datnow <- subset(datnow, datnow$Resolution == resol)
  
  #ces <- unique(datnow$CollisionEnergy)
  #colo <- rainbow(length(ces))
  #datnow <- subset(datnow, datnow$CollisionEnergy == ces[1])
  
  frag <- unique(datnow$FragmentAnnotatedFormula)
  frags <- rainbow(length(frag), start = start, end = end)
  datnow <- subset(datnow, datnow$FragmentAnnotatedFormula == frag[1])
  fragmax <- as.data.frame(c(rep(NA, length(frag))))
  
  if(plot == TRUE){
  png(filename = paste("FragmentsatCEsforCompound", take, ".png", sep = ""));
 
  plot.new()   
  par(fig = c(0,0.85,0,1), new = TRUE)
  
  plot(datnow$CollisionEnergy, log10(datnow$FragmentIntensity), xlim = c(0,100), ylim = c(log10(min(dat$FragmentIntensity)),
                                                                                          log10(max(dat$FragmentIntensity))+1),
       pch = 19, col = frags[1], type = "b", main = paste("Compound", take, "IonMode", ionmode, "Resol", resol, sep = " "),
       xlab = "Collision Energy (NCE)", ylab = "Log10(Intensity)", xaxt = "n");
  axis(side = 1, at = c(0,15,30,45,60,75,90,100))
    
  fragmax[1,1] <- max(datnow$FragmentIntensity)
  
  # legend("topright", legend = frag, pch = 19, col = colo, cex = 0.5);
  
  #   png(filename = paste("FragmentsatCEsforCompound", take, ".png", sep = ""));
  #   plot(datnow$FragmentMZ, datnow$FragmentIntensity, xlim = c(50,max(dat$FragmentMZ)+50), ylim = c(0,max(dat$FragmentIntensity)+1000),
  #        pch = 19, col = colo[1], main = paste("Compound", take, "IonMode", ionmode, "Resol", resol, sep = " "));
  #   legend("topright", legend = ces, pch = 19, col = colo);
  #   
  for(j in 2:length(frag)){
    datnow <- subset(dat, dat$IonMode == ionmode);
    datnow <- subset(datnow, datnow$Resolution == resol)
    datnow <- subset(datnow, datnow$FragmentAnnotatedFormula == frag[j])
    
    points(datnow$CollisionEnergy, log10(datnow$FragmentIntensity), pch = 19, col = frags[j], type = "b")
    fragmax[j,1] <- max(datnow$FragmentIntensity)
  }
  
#   
  fragmax[,2] <- c(rep(1, length(fragmax[,1])))
  par(fig = c(0.7,1,0,1), new = TRUE)
  plot(log10(fragmax[,1]) ~ fragmax[,2], pch = 19, col = frags, axes = FALSE, xlim = c(0.95,1.05), ylim = c(log10(min(dat$FragmentIntensity)),
                                                                                                          log10(max(dat$FragmentIntensity))+1),
       xlab = "merged spectrum", ylab = "")
  box()

  
  dev.off()
  
    }
  }
}