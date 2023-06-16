#--------------
#
# Setup
#
#--------------

library(doParallel)

library(foreach)

library(xtable)


methods_names <- c("DIF1-MAD", "DIF2-SD", "DIF2-TAVC")

## Simulation parameters
##

set.seed(100)

KK <- 500 # number of replications

alpha <- 0.1 # desired coverage

sample_sizes <- c(100, 500, 1000, 2000)

polynomial_degrees <- 0:2

degrees_names <- paste(rep("degree", length(polynomial_degrees)), polynomial_degrees)



##
##

mk_col_names <- function(sample_sizes, n_methods)
{
  tmp <- paste(rep("n =",length(sample_sizes)), sample_sizes)
  
  out <- c()
  
  for (ii in seq_along(tmp))
  {
    out <- c(
      out, 
      rep("",floor((n_methods-1)/2)),
      tmp[ii],
      rep("",ceiling((n_methods-1)/2))
      )
  }
  
  return(out)
}

## Noise types
##

noise_types <- list(
  N1 = function(nn, sd) rnorm(nn, sd = sd),
  N2 = function(nn, sd) sd * sqrt(3/5) * rt(nn, df = 5), 
  N3 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, sd = sd / sqrt(1-0.5**2)), 
  N4 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) sd * sqrt(3/5) * rt(n, df = 5))  
)


## Start cluster
##

cl <- makeCluster(8)

registerDoParallel(cl)



#--------------
#
# Run 
#
#--------------


for (ii in seq_along(noise_types))
{
  out <- matrix(0,length(methods_names)*length(sample_sizes), length(polynomial_degrees))
  
  for (jj in seq_along(sample_sizes))
  {
    for (kk in seq_along(polynomial_degrees))
    {
      out_tmp <- foreach(ll = 1:KK, .combine = "cbind") %dopar%
        {
          library(ChangePointInference)
          
          progress_file <- file("progress.txt")
          
          writeLines(
            paste(
              "running itteration ", ll, 
              "; noise =", names(noise_types)[ii], 
              "; sample size = ", sample_sizes[jj],
              "; poly deg = ", polynomial_degrees[kk]
              ),
            progress_file
            )
          
          xx <- noise_types[[ii]](sample_sizes[jj],1)
          c(
            nrow(diffInf(xx, degree = polynomial_degrees[kk], alpha = alpha, gaussian_noise = TRUE, independent_noise = TRUE, min_scale = sample_sizes[jj]^{(1/3)})$intervals) == 0, 
            nrow(diffInf(xx, degree = polynomial_degrees[kk], alpha = alpha, gaussian_noise = FALSE, independent_noise = TRUE, min_scale = sample_sizes[jj]^{(1/3)})$intervals) == 0, 
            nrow(diffInf(xx, degree = polynomial_degrees[kk], alpha = alpha, gaussian_noise = FALSE, independent_noise = FALSE, min_scale = sample_sizes[jj]^{(1/3)})$intervals) == 0
          )
        }
      
      out[((jj-1)*length(methods_names)+1):(jj*length(methods_names)),kk] <- apply(out_tmp,1,mean)
    }
  }
  
  ## save 
  
  save(out, file = paste("../RData/", names(noise_types)[ii], "-large-coverage-1", sep = ""))
  
  ## post-process
  
  out <- data.frame(round(out,2))
  
  out <- cbind(
    mk_col_names(sample_sizes, length(methods_names)),
    rep(methods_names, length(sample_sizes)),
    out
  )
  
  colnames(out) <- c("","",degrees_names)
  
  ## make table
  
  print(
    xtable(out, align = "|l|l|l|c|c|c|"), 
    file = paste("../tables/", names(noise_types)[ii], "-large-coverage-1.tex", sep = ""),
    include.rownames = FALSE, 
    floating = FALSE
  )
}


stopCluster(cl)

rm(cl)
