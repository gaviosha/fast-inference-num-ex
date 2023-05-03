
## Required packages
##

source("methods-wrappers.R")

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

N1_deg_0 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    source("methods-wrappers.R")
    
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N1, degree = 0"),progress_file)
    
    xx <- rnorm(nn)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out[6] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    out[7] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    
    out[8] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
    
    out 
  }

N1_deg_0 <- apply(N1_deg_0, 1, mean)


## degree = 1 
##

N1_deg_1 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N1, degree = 1"),progress_file)
    
    xx <- rnorm(nn)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly(xx, alpha = alpha, deg = 1)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out
  }

N1_deg_1 <- apply(N1_deg_1, 1, mean)


## degree = 2 
##

N1_deg_2 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N1, degree = 2"),progress_file)
    
    xx <- rnorm(nn)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly(xx, alpha = alpha, deg = 2)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out
  }

N1_deg_2 <- apply(N1_deg_2, 1, mean)


## Save outputs
##

N1_coverage <- cbind(N1_deg_0, N1_deg_1, N1_deg_2)

save(N1_coverage, file = "N1-coverage")



#-------------------
# N2
#-------------------


## degree = 0 
##

N2_deg_0 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    source("methods-wrappers.R")
    
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N2, degree = 0"),progress_file)
    
    xx <- sqrt(3/5) * rt(nn, df = 5)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out[6] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    out[7] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = 20)$intervals) == 0
    
    out[8] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
    
    out 
  }

N2_deg_0 <- apply(N2_deg_0, 1, mean)


## degree = 1 
##

N2_deg_1 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N2, degree = 1"),progress_file)
    
    xx <- sqrt(3/5) * rt(nn, df = 5)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 1)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out
  }

N2_deg_1 <- apply(N2_deg_1, 1, mean)


## degree = 2 
##

N2_deg_2 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N2, degree = 2"),progress_file)
    
    xx <- sqrt(3/5) * rt(nn, df = 5)
    
    out <- numeric(8)
    
    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = 20)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = 20)$intervals) == 0
    
    out[4] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 2)$intervals) == 0 
    
    out[5] <- 0 # Bain and Perron TO-DO
    
    out
  }

N2_deg_2 <- apply(N2_deg_2, 1, mean)


## Save outputs
##

N2_coverage <- cbind(N2_deg_0, N2_deg_1, N2_deg_2)

save(N2_coverage, file = "N2-coverage")

stopCluster(cl)