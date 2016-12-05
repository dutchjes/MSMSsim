wMD <- function(data1, data2, add = "wMD" ) {
  for(i in 1:length(colnames(data2))){
    colnames(data2)[i] <- paste(colnames(data2)[i], add, sep = "")
  }
  data1 <- cbind(data1, data2)
  return(data1)
  
}