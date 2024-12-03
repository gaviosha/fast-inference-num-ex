
#--------------------------------------------
#
# Coverage under strong serial dependence
#
#---------------------------------------------

set.seed(42)

library(ChangePointInference)

library(xtable)


## Helper functions
##

get_ar1_tavc <- function(ww, phi)
{
  tavc <- sum(toeplitz(phi^(0:(ww-1)))) * (1/ww) * (1/(1-phi^2))
  
  return(sqrt(tavc))
}


## Simulation params
##

KK <- 100

nn <- 750

ww <- floor(sqrt(nn)/2.5)

pb <- txtProgressBar(0, KK, style = 3)

phi_seq <- 0.8 + 0:10/100


## Simulation study
##

out_est <- numeric(length(phi_seq))

out_truth <- numeric(length(phi_seq))

for (ii in seq_along(phi_seq))
{
  phi <- phi_seq[ii]
  
  tavc <- get_ar1_tavc(ww, phi)
  
  out_est_tmp <- numeric(KK)
  
  out_truth_tmp <- numeric(KK)
  
  for (jj in 1:KK)
  {
    xx <- arima.sim(model = list(ar = phi), n = nn)
    
    fit_est <- diffInf(yy = xx, degree = 0, gaussian_noise = FALSE, independent_noise = FALSE)
    
    fit_truth <- diffInf(yy = xx, degree = 0, gaussian_noise = FALSE, independent_noise = FALSE, tau = tavc)
    
    out_est_tmp[jj] <- (nrow(fit_est$intervals) == 0)
    
    out_truth_tmp[jj] <- (nrow(fit_truth$intervals) == 0) 
    
    setTxtProgressBar(pb, jj)
  }
  
  out_est[ii] <- mean(out_est_tmp)
  
  out_truth[ii] <- mean(out_truth_tmp)
}

strong_serial_dep <- data.frame(rbind(out_est, out_truth))

rownames(strong_serial_dep) = c("Estimated TAVC", "True TAVC")

colnames(strong_serial_dep) = phi_seq

print(
  xtable(strong_serial_dep, align = "|l|c|c|c|c|c|c|c|c|c|c|c|"),
  file = "../tables/coverage-under-serial-dependence.tex",
  floating = FALSE, 
  include.rownames = FALSE
)


