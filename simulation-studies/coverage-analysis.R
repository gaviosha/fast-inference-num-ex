
library(doParallel)

library(foreach)

library(xtable)

methods_names <- c("DIF1-MAD","DIF2-SD","DIF2-TAVC","NSP","NSP-SN","NSP-AR","B&P","MOSUM (uniscale)", "MOSUM (multiscale)", "SMUCE", "H-SMUCE", "Dep-SMUCE")

degrees_names <- c("degree 0", "degree 1", "degree 2")


#------------------------
#
# Simulation parameters
#
#------------------------

set.seed(42)

nn <- 750 # signal length 

KK <- 500 # number of replications

alpha <- 0.1 # desired coverage

min_width <- floor(sqrt(nn)/2) # minimum segment length

polynomial_degrees <- 0:2

noise_types <- list(
  N1 = function(nn, sd) rnorm(nn, sd = sd),
  N2 = function(nn, sd) sd * sqrt(3/5) * rt(nn, df = 5),
  N3 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, sd = sd / sqrt(1-0.5**2)),
  N4 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) sd * sqrt(3/5) * rt(n, df = 5))
)

cl <- makeCluster(10)

registerDoParallel(cl)



#--------------
#
# Run 
#
#--------------


for (ii in seq_along(noise_types))
{
  
  N_coverage <- matrix(0, 12, length(polynomial_degrees))
  
  for (jj in seq_along(polynomial_degrees))
  {
    N_deg_coverage <- foreach(kk = 1:KK, .combine = "cbind") %dopar%
      {
        source("methods-wrappers.R")
        
        progress_file <- file("progress.txt")
        
        writeLines(paste("running itteration ", kk, "; noise =", names(noise_types)[ii], "degree = ", polynomial_degrees[jj]),progress_file)
        
        xx <- noise_types[[ii]](nn,1)
        
        out <- numeric(12)
        
        out[1] <- nrow(diffInf(xx, degree = polynomial_degrees[jj], gaussian_noise = TRUE, independent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
        out[2] <- nrow(diffInf(xx, degree = polynomial_degrees[jj], gaussian_noise = FALSE, independent_noise = TRUE, alpha = alpha, min_scale = min_width)$intervals) == 0
        out[3] <- nrow(diffInf(xx, degree = polynomial_degrees[jj], gaussian_noise = FALSE, independent_noise = FALSE, alpha = alpha, min_scale = min_width)$intervals) == 0

        out[4] <- nrow(nsp_poly(xx, alpha = alpha, deg = polynomial_degrees[jj])$intervals) == 0
        out[5] <- nrow(nsp_poly_selfnorm(xx, alpha = alpha, deg = polynomial_degrees[jj])$intervals) == 0
        out[6] <- nrow(nsp_poly_ar(xx, alpha = alpha, deg = polynomial_degrees[jj])$intervals) == 0

        if (polynomial_degrees[jj] == 0)
        {
          out[7] <- nrow(Bai_Perron_ints(xx, alpha = alpha, degree = polynomial_degrees[jj], min_width = min_width)$intervals) == 0

          out[8] <- nrow(uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0
          out[9] <- nrow(multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals) == 0

          out[10] <- nrow(smuce_ints(xx, alpha = alpha)$intervals) == 0
          out[11] <- nrow(h_smuce_ints(xx, alpha)$intervals) == 0
          out[12] <- nrow(dep_smuce_ints(xx, alpha)$intervals) == 0
        }
        
        out
      }
    
    N_coverage[,jj] <- apply(N_deg_coverage, 1, mean)
  }
  
  N_coverage <- data.frame(N_coverage)
  
  save(N_coverage, file = paste("../RData/", names(noise_types)[ii], "-coverage", sep = ""))
  
  
  rownames(N_coverage) <- methods_names
  
  colnames(N_coverage) <- degrees_names
  
  N_coverage[8:12,2:3] <- "-"
  
  print(
    xtable(N_coverage, align = "|l|c|c|c|"),
    file = paste("../tables/", names(noise_types)[ii], "-coverage.tex", sep = ""),
    floating = FALSE
  )
  
  rm(N_coverage)
  
  rm(N_deg_coverage)
  
}


stopCluster(cl)

rm(cl)
