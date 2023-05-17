
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



multiscale_mosum_ints <- function(xx, alpha = 0.1, min_width = sqrt(length(xx)))
{
  min_width <- floor(min_width)
  
  bandwidths <- bandwidths.default(length(xx), G.min = min_width)
  
  mlp <- multiscale.localPrune(xx, G = bandwidths)
  
  ci <- confint(mlp, level = alpha)
  
  return(list(intervals = ci$CI[,4:5], threshold = NULL))
}


uniscale_mosum_ints <- function(xx, alpha = 0.1, min_width = sqrt(length(xx)))
{
  min_width <- floor(min_width)
  
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


dep_smuce_ints <- function(xx, alpha = 0.1, min_width = sqrt(length(xx)))
{
  est_sd <- sdrobnormNonparam(x = xx, param = log(min_width) / log(length(xx)))
  
  print(est_sd)
  
  dep_smuce_est <- stepFitNonparam(xx, sd = 1, alpha = alpha)
  # 
  # dep_smuce_est
  
  dep_smuce_ints <- cbind(dep_smuce_est$leftEnd[-1], dep_smuce_est$rightEnd[-length(dep_smuce_est$rightEnd)])

  return(list(intervals = dep_smuce_ints, threshold = NULL))
}


## BP wrapper 
##

Bai_Perron_ints <- function(xx, alpha = 0.1, degree = 0, min_width = sqrt(length(xx))) 
  {
  
  min_width <- floor(min_width)
  
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



if (file.exists("wiener_holder_sim"))
{

    load("wiener_holder_sim")

} else {
  
  wh <- sim_max_holder(100, 500, .03)
  
  save(wh, file = "wiener_holder_sim")
}



## NSP wrapper
##

nsp_selfnorm_ar <- function(xx, alpha = 0.1, degree = 0, ord = 1, M = 1000)
{
  nn <- length(xx)

  x.c <- matrix(0, nn, degree + 1 + ord)

  if (degree > 0) x.c[,2:(degree+1)] <- poly(nn, degree, raw = TRUE)

  for (kk in 1:ord) x.c[(1+kk):nn, degree + 1 + kk] <- xx[1:(nn-kk)]

  x.c <- x.c[(ord+1):nn,]

  xx <- xx[(ord+1):nn]

  thresh.val <- as.numeric(stats::quantile(wh, 1-alpha))

    nsp_selfnorm(y = xx, x = x.c, M = M, lambda = thresh.val)
}

