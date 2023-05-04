
library(nsp)

library(lpSolve)

library(mosum)

library(lpSolve)

library(stepR)

library(strucchange)

library(ChangePointInference)

library(stats)


## MOSUM wrapper
##

multiscale_mosum_ints <- function(xx, alpha, min_width)
{
  bandwidths <- bandwidths.default(length(xx))
  
  mlp <- multiscale.localPrune(xx, G = bandwidths)
  
  ci <- confint(mlp, level = alpha)
  
  return(list(intervals = ci$CI[,4:5], threshold = NULL))
}


uniscale_mosum_ints <- function(xx, alpha, min_width)
{
  mlp <- multiscale.localPrune(xx, G = min_width)
  
  ci <- confint(mlp, level = alpha)
  
  return(list(intervals = ci$CI[,4:5], threshold = NULL))
}


## SMUCE wrapper
##

smuce_ints <- function(xx, alpha)
{
  smuce_fit <- smuceR(xx, family = "gauss", confband = TRUE, alpha = alpha)
  
  smuce_ints <- jumpint(smuce_fit)[-1,1:2]
  
  return(list(intervals = smuce_ints, threshold = NULL))
}


## BP wrapper 
##

Bai_Perron_ints <- function(xx, alpha = 0.1, degree = 0, min_width) 
  {
  
  if (degree == 0) cpt_est <- breakpoints(xx ~ 1, h = max(min_width,degree+2), breaks = NULL)
  
  if (degree > 0) cpt_est <- breakpoints(xx ~ poly(seq_along(xx),degree), h = max(min_width,degree+2), breaks = NULL)
  
  
  if (is.na(cpt_est$breakpoints[1])) return(list(intervals = matrix(0,0,2), threshold = NULL))
  
  
  NN <- length(cpt_est$breakpoints)
  
  BP_ints <- try(confint(cpt_est, level = 1 - alpha / NN), silent = TRUE)
  
  if (class(BP_ints) == "try-error")
  {
    BP_ints <- matrix(0,0,2)
  } else {
    BP_ints <- BP_ints$confint[,c(1,3),drop=FALSE] 
  }
  
  return(list(intervals = BP_ints, threshold = NULL))
}
