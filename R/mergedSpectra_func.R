## merged spectra function

mergedSpectra <- function(all.frag, intensity.type = "absolute", dmz = 0.001){
  
  if(intensity.type != "absolute" & intensity.type != "relative"){
    stop("Which type of intensity: absolute or relative?")
  }
  
  if(intensity.type == "absolute"){
    find <- "FragmentIntensity"
  }
  
  if(intensity.type == "relative"){
    find <- "FragmentRelativeIntensity"
  }
  
  IDs <- as.numeric(unique(all.frag$CompoundID))
  IDs <- IDs[order(IDs, decreasing = FALSE)]
  merged.spec <- list();
  
  for(i in 1:length(IDs)){
    
    where <- which(all.frag$CompoundID == IDs[i])
    merged.spec[[i]] <- as.data.frame(all.frag[where,c("FragmentMZ", paste(find))]) ## add in here if you want extra column from all.frag
    colnames(merged.spec[[i]]) <- c("mz", "int")

    merged.spec[[i]] <- filterMerge(merged.spec[[i]], dmz = dmz)
  }
  
  names(merged.spec) <- IDs
  return(merged.spec)
}