
#--------------
#
# Setup
#
#--------------


library(doParallel)

library(foreach)

library(xtable)


## Simulation parameters
##

set.seed(100)

KK <- 500 # number of replications

alpha <- 0.1 # desired coverage


## Noise types
##

noise_types <- list(
  N1 = function(nn, sd) rnorm(nn, sd = sd), 
  N2 = function(nn, sd) sd * sqrt(3/5) * rt(nn, df = 5), 
  N3 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, sd = sd / sqrt(1-0.5**2)), 
  N4 = function(nn, sd) arima.sim(model = list(ar = 0.5), n = nn, rand.gen = function(n) sd * sqrt(3/5) * rt(n, df = 5))  
)


## Evaluate intervals 
##

metrics_names <- c("no. genuine", "prop. genuine", "length", "coverage")

methods_names <- c(
  "DIF1-MAD","DIF2-SD","DIF2-TAVC",
  "NSP","NSP-SN","NSP-AR",
  "B&P",
  "MOSUM (uniscale)","MOSUM (multiscale)", 
  "SMUCE","H-SMUCE","Dep-SMUCE"
  )


special_pad <- function(str_vec, pad1, pad2)
{
  out <- c()
  
  for (ii in seq_along(str_vec)) out <- c(out, rep("",pad1), str_vec[ii], rep("",pad2))
  
  return(out)
}

evaluate_ints <- function(cpt_ints, cpt_locs)
{
  
  out <- c()
  
  eval_list <- list(
    
    no.genuine = function(cpt_ints, cpt_locs) 
    {
      if (nrow(cpt_ints) == 0) return(0)
      apply(cpt_ints, 1, function(ii) any(cpt_locs %in% ii[1]:ii[2])) |> sum()
    },
    
    prop.genuine = function(cpt_ints, cpt_locs) 
    {
      if (nrow(cpt_ints) == 0) return(0)
      (apply(cpt_ints, 1, function(ii) any(cpt_locs %in% ii[1]:ii[2])) |> sum()) / nrow(cpt_ints)
    }, 
    
    int.length = function(cpt_ints, cpt_locs)
    {
      if (nrow(cpt_ints) == 0) return(0)
      (apply(cpt_ints, 1, diff) + 1) |> mean()
    }, 
    
    coverage = function(cpt_ints, cpt_locs)
    {
      if (nrow(cpt_ints) == 0) return(1)
      apply(cpt_ints, 1, function(ii) any(cpt_locs %in% ii[1]:ii[2])) |> all()
    }
  )
  
  for (eval in eval_list) out <- c(out, eval(cpt_ints,cpt_locs))
  
  return(out)
}


unpack <- function(lst)
{
  
  d1 <- nrow(lst[[1]])
  d2 <- ncol(lst[[2]])
  d3 <- length(lst)
  
  aa <- array(do.call(cbind, lst), dim = c(d1, d2, d3))
  
  mm <- apply(aa, c(1,2), mean)
  
  as.vector(mm)
}


## Do par 
##

cl <- makeCluster(6)

registerDoParallel(cl)



# ----------------------------
# 
# The blocks signal
# 
# ----------------------------


print("The blocks signal")


## Signal specific params
##

blocks_signal <- c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40))

blocks_cpt <- c(205, 267, 308, 472)

nn <- length(blocks_signal)

min_width <- floor(sqrt(nn))

blocks_deg <- 0

blocks_sd <- \(ii) ifelse(ii %in% 1:2, 10, 5)

blocks_out <- c()


## Simulation
##

for (ii in seq_along(noise_types))
{
  foreach(kk = 1:KK) %dopar%
    {
      source("methods-wrappers.R")

      progress_file <- file("progress.txt")

      writeLines(paste("running itteration ", kk, "; noise =", names(noise_types)[ii]),progress_file)

      xx <- blocks_signal + noise_types[[ii]](length(blocks_signal), blocks_sd(ii))

      out <- matrix(0,4,12)

      diff1_mad <- diffInf(xx, degree = blocks_deg, alpha = alpha, gaussian_noise = TRUE, independent_noise = TRUE)$intervals[,1:2]
      diff2_sd <- diffInf(xx, degree = blocks_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = TRUE, min_scale = min_width)$intervals[,1:2]
      diff2_tavc <- diffInf(xx, degree = blocks_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = FALSE, min_scale = min_width)$intervals[,1:2]

      out[,1] <- evaluate_ints(diff1_mad, blocks_cpt)
      out[,2] <- evaluate_ints(diff2_sd, blocks_cpt)
      out[,3] <- evaluate_ints(diff2_tavc, blocks_cpt)

      nsp_ <- nsp_poly(xx, deg = blocks_deg, alpha = alpha)$intervals[,1:2]
      nsp_sn <- nsp_poly_selfnorm(xx, deg = blocks_deg, alpha = alpha)$intervals[,1:2]
      nsp_ar <- nsp_poly_ar(xx, deg = blocks_deg, alpha = alpha)$intervals[,1:2]

      out[,4] <- evaluate_ints(nsp_, blocks_cpt)
      out[,5] <- evaluate_ints(nsp_sn, blocks_cpt)
      out[,6] <- evaluate_ints(nsp_ar, blocks_cpt)

      bp <- Bai_Perron_ints(xx, alpha = alpha, degree = blocks_deg, min_width = min_width)$intervals[,1:2]
      out[,7] <- evaluate_ints(bp, blocks_cpt)

      mosum_uni <- uniscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals[,1:2]
      mosum_multi <- multiscale_mosum_ints(xx, alpha = alpha, min_width = min_width)$intervals[,1:2]

      out[,8] <- evaluate_ints(mosum_uni, blocks_cpt)
      out[,9] <- evaluate_ints(mosum_multi, blocks_cpt)

      smuce <- smuce_ints(xx, alpha = alpha)$intervals[,1:2]
      dep_smuce <- dep_smuce_ints(xx, alpha = alpha)$intervals[,1:2]
      h_smuce <- h_smuce_ints(xx, alpha = alpha)$intervals[,1:2]

      out[,10] <- evaluate_ints(smuce, blocks_cpt)
      out[,11] <- evaluate_ints(dep_smuce, blocks_cpt)
      out[,12] <- evaluate_ints(h_smuce, blocks_cpt)

      out

    } -> out

  blocks_out<- cbind(blocks_out, unpack(out))
}


## Save outputs
##

blocks_out <- data.frame(blocks_out)

save(blocks_out, file = "RData/blocks-performance")

blocks_out <- cbind(
  special_pad(methods_names,1,2),
  rep(metrics_names, 12),
  blocks_out
)

colnames(blocks_out) <- c(rep("",2),names(noise_types))

print(
  xtable(blocks_out, align = "|l|c|c|c|c|c|c|"),
  file = "tables/blocks-performance.tex",
  floating = FALSE,
  include.rownames = FALSE
)



#---------------------
#
# The waves signal
#
#---------------------


print("The waves signal")


## Signal specific params
##

waves_signal <- c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3))

waves_cpt <- c(150,300,450)

nn <- length(waves_signal)

min_width <- floor(sqrt(nn))

waves_deg <- 1

waves_sd <- 5

waves_out <- c()


## Simulation
##

for (ii in seq_along(noise_types))
{
  foreach(kk = 1:KK) %dopar%
    {
      source("methods-wrappers.R")
      
      progress_file <- file("progress.txt")
      
      writeLines(paste("running itteration ", kk, "; noise =", names(noise_types)[ii]),progress_file)
      
      xx <- waves_signal + noise_types[[ii]](length(waves_signal), waves_sd) 
      
      out <- matrix(0,4,6)
      
      diff1_mad <- diffInf(xx, degree = waves_deg, alpha = alpha, gaussian_noise = TRUE, independent_noise = TRUE)$intervals[,1:2]
      diff2_sd <- diffInf(xx, degree = waves_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = TRUE, min_scale = min_width)$intervals[,1:2]
      diff2_tavc <- diffInf(xx, degree = waves_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = FALSE, min_scale = min_width)$intervals[,1:2]
      
      out[,1] <- evaluate_ints(diff1_mad, waves_cpt)
      out[,2] <- evaluate_ints(diff2_sd, waves_cpt)
      out[,3] <- evaluate_ints(diff2_tavc, waves_cpt)
      
      nsp_ <- nsp_poly(xx, deg = waves_deg, alpha = alpha)$intervals[,1:2]
      nsp_sn <- nsp_poly_selfnorm(xx, deg = waves_deg, alpha = alpha)$intervals[,1:2]
      nsp_ar <- nsp_poly_ar(xx, deg = waves_deg, alpha = alpha)$intervals[,1:2]

      out[,4] <- evaluate_ints(nsp_, waves_cpt)
      out[,5] <- evaluate_ints(nsp_sn, waves_cpt)
      out[,6] <- evaluate_ints(nsp_ar, waves_cpt)

      out
      
    } -> out
  
  waves_out <- cbind(waves_out, unpack(out))
}


## Save outputs
##

waves_out <- data.frame(waves_out)

save(waves_out, file = "RData/waves-performance")

waves_out <- cbind(
  special_pad(methods_names[1:6],1,2), 
  rep(metrics_names, 6), 
  waves_out
  ) 

colnames(waves_out) <- c(rep("",2),names(noise_types))

print(
  xtable(waves_out, align = "|l|c|c|c|c|c|c|"),
  file = "tables/waves-performance.tex",
  floating = FALSE, 
  include.rownames = FALSE
)


#------------------------
#
# The hills signal
#
#------------------------


print("The hills signal")


## Signal specific params
##

hills_signal <- 10 * rep((1:100/100) * (1 - (1:100/100)), 4)

hills_cpt <- c(100,200,300)

nn <- length(hills_signal)

min_width <- floor(sqrt(nn))

hills_deg <- 2

hills_sd <- 1

hills_out <- c()


## Simulation
##

for (ii in seq_along(noise_types))
{
  foreach(kk = 1:KK) %dopar%
    {
      source("methods-wrappers.R")
      
      progress_file <- file("progress.txt")
      
      writeLines(paste("running itteration ", kk, "; noise =", names(noise_types)[ii]),progress_file)
      
      xx <- hills_signal + noise_types[[ii]](length(hills_signal), hills_sd) 
      
      out <- matrix(0,4,6)
      
      diff1_mad <- diffInf(xx, degree = hills_deg, alpha = alpha, gaussian_noise = TRUE, independent_noise = TRUE)$intervals[,1:2]
      diff2_sd <- diffInf(xx, degree = hills_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = TRUE, min_scale = min_width)$intervals[,1:2]
      diff2_tavc <- diffInf(xx, degree = hills_deg, alpha = alpha, gaussian_noise = FALSE, independent_noise = FALSE, min_scale = min_width)$intervals[,1:2]
      
      out[,1] <- evaluate_ints(diff1_mad, hills_cpt)
      out[,2] <- evaluate_ints(diff2_sd, hills_cpt)
      out[,3] <- evaluate_ints(diff2_tavc, hills_cpt)
      
      nsp_ <- nsp_poly(xx, deg = hills_deg, alpha = alpha)$intervals[,1:2]
      nsp_sn <- nsp_poly_selfnorm(xx, deg = hills_deg, alpha = alpha)$intervals[,1:2]
      nsp_ar <- nsp_poly_ar(xx, deg = hills_deg, alpha = alpha)$intervals[,1:2]

      out[,4] <- evaluate_ints(nsp_, hills_cpt)
      out[,5] <- evaluate_ints(nsp_sn, hills_cpt)
      out[,6] <- evaluate_ints(nsp_ar, hills_cpt)

      out
      
    } -> out
  
  hills_out <- cbind(hills_out, unpack(out))
}


## Save outputs
##

hills_out <- data.frame(hills_out)

save(hills_out, file = "RData/hills-performance")

hills_out <- cbind(
  special_pad(methods_names[1:6],1,2), 
  rep(metrics_names, 6), 
  hills_out
) 

colnames(hills_out) <- c(rep("",2),names(noise_types))

print(
  xtable(hills_out, align = "|l|c|c|c|c|c|c|"),
  file = "tables/hills-performance.tex",
  floating = FALSE, 
  include.rownames = FALSE
)


stopCluster(cl)

rm(cl)

