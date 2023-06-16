
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

h_smuce_ints <- function(xx, alpha)
{
  smuce_fit <- stepFit(xx, family = "hsmuce", confband = TRUE, alpha = alpha)
  
  smuce_ints <- jumpint(smuce_fit)[-1,1:2]
  
  return(list(intervals = smuce_ints, threshold = NULL))
}

sdrobnormNonparam <- function (x, param, supressWarningNA = FALSE, supressWarningResultNA = FALSE) 
{
  if (!is.logical(supressWarningNA) || length(supressWarningNA) != 1 || is.na(supressWarningNA)) 
    {
    stop("supressWarningNA must be a single logical (not NA)")
    }
  
  if (!is.logical(supressWarningResultNA) || length(supressWarningResultNA) != 1 || is.na(supressWarningResultNA)) 
    {
    stop("supressWarningResultNA must be a single logical (not NA)")
    }
  
  if (any(is.na(x)))
    {
    x <- x[!is.na(x)]
    if (!supressWarningNA) {
      warning("the data vector 'x' contains NAs")
    }
    }
  
  n = length(x)
  kn = floor(n^(param))
  mn = floor(n/kn)
  A = rep(0, mn)
  for (i in (0:(mn - 1))) {
    for (j in (1:kn)) {
      A[i + 1] = A[i + 1] + (1/kn) * x[j + i * kn]
    }
  }
  
  B = 0
  
  for (i in 1:(mn - 1)) B = B + (abs(A[i + 1] - A[i]))^2
  
  return(sqrt(kn)/(sqrt(2 * (mn - 1))) * sqrt(B))
}


dep_smuce_ints <- function(xx, alpha = 0.1, min_width = (length(xx))^(1/3))
{
  est_sd <- sdrobnormNonparam(x = xx, param = log(min_width) / log(length(xx)))
  
  smuce_fit <- smuceR(xx, family = "gauss", confband = TRUE, alpha = alpha, param = est_sd, lengths = floor(min_width):length(xx))
  
  smuce_ints <- jumpint(smuce_fit)[-1,1:2]
  
  return(list(intervals = smuce_ints, threshold = NULL))
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

