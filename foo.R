
## Required packages
##

library(ChangePointInference)

source("R/other-methods-wrappers.R")

library(doParallel)

library(foreach)



nn <- 500 # signal length 

KK <- 10 # number of replications

alpha <- 0.1 # desired coverage

cl <- makeCluster(8)

registerDoParallel(cl)



#-------------------
# N1
#-------------------


## degree = 0 
##

foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    source("R/other-methods-wrappers.R")
    
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N1, degree = 0"),progress_file)
    
    xx <- rnorm(nn)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = sqrt(nn))$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly(xx, alpha = alpha)$intervals) == 0 
    
    out[5] <- 0 
    
    out[6] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    out[7] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    
    out[8] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
    
    out 
  } -> vv


## degree = 1 
##

foreach(ii=1:KK, .combine = "rbind") %dopar%
  {
    
  }

## degree = 2 
##

foreach(ii=1:KK, .combine = "rbind") %dopar%
  {
    
  }


stopCluster(cl)


file_conn <- file("progress.txt")
  
writeLines(paste("On itteration: ", 10),file_conn)



stopCluster(cl)
