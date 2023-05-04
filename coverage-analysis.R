
## Required packages
##

source("methods-wrappers.R")

library(doParallel)

library(foreach)

library(xtable)

diff_names <- c("DIF1-MAD","DIF1-TAVC","DIF2-MAD","DIF2-TAVC")

other_names <- c("B&P","MOSUM (uniscale)", "MOSUM (multiscale)", "SMUCE")

degrees <- c("degree 0", "degree 1", "degree 2")

## Simulation params 
##

set.seed(42)

nn <- 750 # signal length 

KK <- 100 # number of replications

alpha <- 0.1 # desired coverage

min_width <- 20 # minimum segment length

cl <- makeCluster(11)

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

    out <- numeric(11)

    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 0, min_width = min_width)$intervals) == 0

    out[9] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    out[10] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0

    out[11] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0

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

    out <- numeric(11)

    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 1)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 1, min_width = min_width)$intervals) == 0

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

    out <- numeric(9)

    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 2)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 2, min_width = min_width)$intervals) == 0

    out
  }

N1_deg_2 <- apply(N1_deg_2, 1, mean)


## Save outputs
##

N1_coverage <- data.frame(cbind(N1_deg_0, N1_deg_1, N1_deg_2))

save(N1_coverage, file = "RData/N1-coverage")


rownames(N1_coverage) <- c(diff_names, "NSP", other_names)

colnames(N1_coverage) <- degrees

N1_coverage[9:11,2:3] <- "-"

print(
  xtable(N1_coverage, align = "|l|c|c|c|"),
  file = "tables/N1-coverage.tex",
  floating = FALSE
)

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

    out <- numeric(11)

    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 0, min_width = min_width)$intervals) == 0

    out[9] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    out[10] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0

    out[11] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0

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

    out <- numeric(11)

    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 1)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 1, min_width = min_width)$intervals) == 0

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

    out <- numeric(11)

    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 2)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 2, min_width = min_width)$intervals) == 0

    out
  }

N2_deg_2 <- apply(N2_deg_2, 1, mean)


## Save outputs
##

N2_coverage <- data.frame(cbind(N2_deg_0, N2_deg_1, N2_deg_2))

save(N2_coverage, file = "RData/N2-coverage")


rownames(N2_coverage) <- c(diff_names, "NSP-SN", other_names)

colnames(N2_coverage) <- degrees

N2_coverage[9:11,2:3] <- "-"

print(
  xtable(N2_coverage, align = "|l|c|c|c|"),
  file = "tables/N2-coverage.tex",
  floating = FALSE
)


#-------------------
# N3
#-------------------


## degree = 0 
##

N3_deg_0 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    source("methods-wrappers.R")
    
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N3, degree = 0"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn)
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    
    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 0, min_width = min_width)$intervals) == 0
    
    out[9] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    out[10] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    
    out[11] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
    
    out 
  }

N3_deg_0 <- apply(N3_deg_0, 1, mean)


## degree = 1 
##

N3_deg_1 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N3, degree = 1"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn)
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0

    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 1)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 1, min_width = min_width)$intervals) == 0
    
    out
  }

N3_deg_1 <- apply(N3_deg_1, 1, mean)


## degree = 2 
##

N3_deg_2 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N3, degree = 2"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn)
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    
    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 2)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 2)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 2, min_width = min_width)$intervals) == 0
    
    out
  }

N3_deg_2 <- apply(N3_deg_2, 1, mean)


## Save outputs
##

N3_coverage <- data.frame(cbind(N3_deg_0, N3_deg_1, N3_deg_2))

save(N3_coverage, file = "RData/N3-coverage")


rownames(N3_coverage) <- c(diff_names, "NSP-AR", other_names)

colnames(N3_coverage) <- degrees

N3_coverage[9:11,2:3] <- "-"

print(
  xtable(N3_coverage, align = "|l|c|c|c|"), 
  file = "tables/N3-coverage.tex",
  floating = FALSE
  )


#-------------------
# N4
#-------------------


## degree = 0 
##

N4_deg_0 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    source("methods-wrappers.R")
    
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N4, degree = 0"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) rt(n, df = 5))
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 0, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 0, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    
    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 0, min_width = min_width)$intervals) == 0
    
    out[9] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    out[10] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
    
    out[11] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
    
    out 
  }

N4_deg_0 <- apply(N4_deg_0, 1, mean)


## degree = 1 
##

N4_deg_1 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N4, degree = 1"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) rt(n, df = 5))
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 1, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 1, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    
    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 1)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 1, min_width = min_width)$intervals) == 0
    
    out
  }

N4_deg_1 <- apply(N4_deg_1, 1, mean)


## degree = 2 
##

N4_deg_2 <- foreach(ii=1:KK, .combine = "cbind") %dopar%
  {
    progress_file <- file("progress.txt")
    
    writeLines(paste("running itteration ",ii, "; noise = N4, degree = 2"),progress_file)
    
    xx <- arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) rt(n, df = 5))
    
    out <- numeric(11)
    
    out[1] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = FALSE, alpha = alpha)$intervals) == 0 
    out[2] <- nrow(diffInf(xx, degree = 2, noise_type = "gaussian", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[3] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0 
    out[4] <- nrow(diffInf(xx, degree = 2, noise_type = "non_gaussian_dependent", dependent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
    
    out[5] <- nrow(nsp_poly(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[6] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = 0)$intervals) == 0
    out[7] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = 0)$intervals) == 0
    
    out[8] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = 0, min_width = min_width)$intervals) == 0
    
    out
  }

N4_deg_2 <- data.frame(apply(N4_deg_2, 1, mean))


## Save outputs
##

N4_coverage <- cbind(N4_deg_0, N4_deg_1, N4_deg_2)

save(N4_coverage, file = "RData/N4-coverage")


rownames(N4_coverage) <- c(diff_names,"NSP-AR" ,"NSP-SN", other_names)

colnames(N4_coverage) <- degrees

N4_coverage[9:11,2:3] <- "-"

print(
  xtable(N4_coverage, align = "|l|c|c|c|", digits = c(0,0,0,0)), 
  file = "tables/N4-coverage.tex",
  floating = FALSE
  )



## Stop cluster
##

stopCluster(cl)
