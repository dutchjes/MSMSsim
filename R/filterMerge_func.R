#' Merging of spectra
#'
#' @param peaks Data frame with the peaks which are to be merged. Column 1 should be mz, column 2 should be absolute intensity
#' @param dmz Delta m/z requested for merging, in Da. Default is 0.001
#' @param int Intensity cutoff for merging. Currently not applied
#'
#' @return data frame with merged peaks. Column 1 is mz, column 2 is absolute intensity, and a new column 3 is the frequency of the fragment. Data frame is ordered by frequency
#' @export
#'
#' @examples
#' 
filterMerge <- function (peaks, dmz = 0.001, int = 1e1) 
{
  cutoff_int_limit <- int
  cutoff_mz_limit <- dmz
  peaks <- data.frame(mz = peaks[,1], int = peaks[,2])
  peaks_o <- peaks[order(peaks$int, decreasing = TRUE), ]
  n <- 1
  peaks_o$count <- 0
  while (n <= nrow(peaks_o)) {
    
    windowpeaks <- which( (peaks_o$mz > peaks_o[n, "mz"] - cutoff_mz_limit) & (peaks_o$mz < peaks_o[n, "mz"] + cutoff_mz_limit))
    if(length(windowpeaks) > 1) {
      peaks_o <- peaks_o[-windowpeaks[-1],]
    }
    peaks_o[n, "count"] <- length(windowpeaks)
    n <- n + 1
  }
  return(peaks_o[order(peaks_o$mz), ])
}
