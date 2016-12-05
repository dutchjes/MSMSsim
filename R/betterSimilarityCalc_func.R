

### Merged Fragment Table for Comparison
## t is the mz tolerance for merging fragments
## replace.na is the value for undetected fragments

#' Builds the merged fragment table 
#'
#' @param spec.top MSMS spectrum 1
#' @param spec.bottom MSMS spectrum 2
#' @param t m/z tolerance for merging fragments. Default is 0.005
#' @param replace.na value to replace NAs with. Default is 0
#'
#' @return data frame with fragment m/z values, intensity in MSMS spectrum 1, and intensity in MSMS spectrum 2
#' @export
#'
#' @examples
#' 
simTable <- function(spec.top, spec.bottom, t = 0.005, replace.na = 0){
  
  spec.top <- data.frame(mz = spec.top[, 1], intensity = spec.top[,2] ) 

  spec.bottom <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[,2] )

  for (i in 1:nrow(spec.bottom)) 
    spec.top[, 1][spec.bottom[, 1][i] >= spec.top[,1] - t & spec.bottom[, 1][i] <= spec.top[, 1] + t] <- spec.bottom[,1][i]
  
  alignment <- merge(spec.top, spec.bottom, by = 1, all = TRUE)
  colnames(alignment) <- c("mz","intensity.top", "intensity.bottom")
#  if (length(unique(alignment[, 1])) != length(alignment[, 1])) 
#    warning("the m/z tolerance is set too high")
  alignment[is.na(alignment)] <- replace.na
  return(alignment)
}

### Relative cutoff
#' Relative intensity cutoff
#'
#' @param alignment merged fragment table, for example output from simTable() 
#' @param b relative intensity cutoff, default is 0.05. all intensities below this value will be replaced with 0. 
#'
#' @return data frame with fragment m/z values, intensity in MSMS spectrum 1, and intensity in MSMS spectrum 2
#' @export
#'
#' @examples
#' 
relCutoff <- function(alignment, b = 0.05, replace.na = 0){
  
  alignment <- data.frame(mz = alignment[,1], intensity.top = alignment[,2], intensity.bottom = alignment[,3])
  
  alignment$int.top.normalized <- alignment$intensity.top / max(alignment$intensity.top) * 100
  alignment$int.bottom.normalized <- alignment$intensity.bottom / max(alignment$intensity.bottom) * 100
  
  below <- which(alignment$int.top.normalized <= b)
  alignment[below,"intensity.top"] <- replace.na
  
  below <- which(alignment$int.bottom.normalized <= b)
  alignment[below,"intensity.bottom"] <- replace.na
  
  alignment <- alignment[,c(1:3)]
  
  match <-  as.vector(rep("NA", nrow(alignment)))
  for(i in 1:nrow(alignment)){
    if(alignment[i, "intensity.top"] != 0 && alignment[i, "intensity.bottom"] != 0){
      match[i] <- as.integer(1)
    }else{
      match[i] <- 0
    }
  }
  
  alignment$match <- as.numeric(match)
  
  return(alignment)
}

### Weighting of the mz and intensity
## y is for m/z weighting
## z is for intensity weighting

#' Weighting of m/z and intensity
#'
#' @param alignment merged fragment table, for example output from simTable() or relCutoff()
#' @param y weighting factor for m/z. Default is 0
#' @param z weight factor for intensity. Default is 1
#' @details In order to compare spectra in MSsim function, vectors need to be generated for each spectra.
#' These vectors can be generated from only the intensities (default; y = 0, z = 1) or using different weighting factors. For example,
#' the dot product similarity scoring from the NIST database uses y = 3, z = 0.6 and from MassBank y = 2, z = 0.5.
#' Here the input intensity values for each spectrum are replaced by the weighted vector. The weighted vector is calculated as:
#' 
#' (m/z)^y * (intensity)^z
#' 
#' @return data frame with fragment m/z values, weighted vector for MSMS spectrum 1, and weighted vector of MSMS spectrum 2
#' @export
#'
#' @examples
#' 
simTableWeighting <- function(alignment, y = 0, z = 1){

  alignment <- data.frame(mz = alignment[,1], intensity.top = alignment[,2], intensity.bottom = alignment[,3])
  alignment$top.weighted <- alignment[,1]^y  * alignment[,2]^z
  alignment$bottom.weighted <- alignment[,1]^y * alignment[,3]^z
  
  alignment <- alignment[,c(1,4,5)]
  
  return(alignment)
} 


### Similarity score calculation

#' Similarity score calculation
#'
#' @param alignment merged fragment table, for example output from simTable(), relCutoff(), or simTableWeighting(). 
#' The two vectors which are being compared should be in columns 2 and 3
#'
#' @details correlation coeffiecient is calculated as (u %*% v) / ( sqrt(sum(u*u)) * sqrt(sum(v*v)) ) where u in the first input vector and
#' v is the second input vector and %*% is the dot product
#' @return correlation coeffiecient of the input vectors
#' @export
#'
#' @examples
MSsim <- function(alignment){
  
  u <- alignment[,2]
  v <- alignment[,3]
  score <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))
  return(score)
  
}


### Wrapper function

#' Calculates similarity score of two MSMS spectra using modified cosine
#'
#' @param spec.top MSMS spectrum 1
#' @param spec.bottom MSMSpectrum 2
#' @param y m/z weighting factor for simTableWeighting(). Default is 0
#' @param z intensity weighting factor for simTableWeighting(). Default is 1
#' @param t m/z tolerance for merging factors for simTable(). Default is 0.005
#' @param replace.na value to replace NAs with. Default is 0
#' @param b relative intensity cutoff for relCutoff(). Default is 0.05.
#'
#' @return List with the correlation coefficient (ie the similarity score) and data frame showing output from each step (ie merging, intensity
#' cutoff, and weighting)
#' @export
#'
#' @examples
#' 
OrgMSSim <- function(spec.top, spec.bottom, y = 0, z = 1, t = 0.005, replace.na = 0,  b = 0.05, min.frag = 0){
  
  result <- simTable(spec.top = spec.top, spec.bottom = spec.bottom, t = t, replace.na = replace.na)
  CUresult <- relCutoff(result, b = b)
  weighted <- simTableWeighting(CUresult, y = y, z = z)
  score <- MSsim(weighted)
  
  if(sum(CUresult$match) < min.frag){
    score <- 0
  }
  
  result <- cbind(weighted, CUresult[,c(2,3,4)], result[,c(2,3)])
  
  return(list(score=score, result=result))
}



#' Wrapper function
#' 
#' Calculates similarity score of two MSMS spectra using count of matched fragments
#'
#' @param spec.top MSMS spectrum 1
#' @param spec.bottom MSMSpectrum 2
#' @param y m/z weighting factor for simTableWeighting(). Default is 0
#' @param z intensity weighting factor for simTableWeighting(). Default is 1
#' @param t m/z tolerance for merging factors for simTable(). Default is 0.005
#' @param replace.na value to replace NAs with. Default is 0
#' @param b relative intensity cutoff for relCutoff(). Default is 0.05.
#' @param min.frag the minimum number of matched fragments required. If number of matched fragments is less than this value, 
#' score is automatically 0. Default is 0
#'
#' @return List with the the similarity score (sum of number of matched fragments) and data frame showing output from each step
#'  (i.e., merging and intensity cutoff). First column is m/z, then the top and bottom weighted values used for the similarity
#'  score calculation. Fourth and fifth columns show the result after applying the relative intensity cutoff, followed by
#'  if fragments are matching. Final two column show the original input intensities prior to any processing.
#' @export
#'
#' @examples
#' 
fragCountSim <- function(spec.top, spec.bottom, y = 0, z = 1, t = 0.005, replace.na = 0,  b = 0.05, min.frag = 0){
  
  result <- simTable(spec.top = spec.top, spec.bottom = spec.bottom, t = t, replace.na = replace.na)
  CUresult <- relCutoff(result, b = b)
 # weighted <- simTableWeighting(CUresult, y = y, z = z)
 # score <- MSsim(weighted)
  
  score <- sum(CUresult$match)
  
  if(sum(CUresult$match) < min.frag){
    score <- 0
  }
  
 
  result <- cbind(CUresult, result[,c(2,3)])
  
  return(list(score=score, result=result))
}


