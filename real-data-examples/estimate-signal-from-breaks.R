
estimate_signal_from_breaks <- function(data, cpt_locs)
{
  fhat <- data * 0
  
  brks <- c(0, cpt_locs, length(data))
  
  for(ii in 1:(length(cpt_locs) + 1))
  {
    int <- (brks[ii] + 1):brks[ii + 1]
    fhat[int] <- mean(data[int]) 
  }
  
  return(fhat)
}
