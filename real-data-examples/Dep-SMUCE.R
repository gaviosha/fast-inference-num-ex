
library(stepR)

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



dep_smuce <- function(xx, alpha = 0.1, min_width = (length(xx))^(1/3))
{
  est_sd <- sdrobnormNonparam(x = xx, param = log(min_width) / log(length(xx)))
  
  smuce_fit <- smuceR(xx, family = "gauss", confband = TRUE, alpha = alpha, param = est_sd, lengths = floor(min_width):length(xx))
  
  return(smuce_fit)
}