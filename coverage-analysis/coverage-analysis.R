
#-------------------------------
#
# Coverage simulations
#
#-------------------------------




## Simulation parameters
##

set.seed(100)

N_draws <- 1000

signal_length <- 500

coverage_levels <- c(0.10, 0.05, 0.01)

pb <- txtProgressBar(min = 0, max = N_draws, style = 3)


## simulated thresholds for NSP
##

nsp_sims_deg0 <- cov_dep_multi_norm_poly(signal_length, 0, 1000)

nsp_sims_deg1 <- cov_dep_multi_norm_poly(signal_length, 1, 1000)

nsp_sims_deg2 <- cov_dep_multi_norm_poly(signal_length, 2, 1000)


## Degree 0 polynomials
##

print("Degree 0 polynomials")

zz <- matrix(0, 9, length(coverage_levels))

for (ii in 1:N_draws)
{
  xx <- rnorm(signal_length)
  
  for (jj in seq_along(coverage_levels))
  {
    # zz[1,jj] <- zz[1,jj] + as.numeric(nrow(nsp_poly(xx, alpha = coverage_levels[jj])$intervals) > 0)
    # zz[2,jj] <- zz[2,jj] + as.numeric(nrow(nsp_poly(xx, thresh.val = quantile(nsp_sims_deg0, 1 - coverage_levels[jj]))$intervals) > 0)
    zz[3,jj] <- zz[3,jj] + as.numeric(nrow(mosum_ints(xx, alpha = coverage_levels[jj], bandwidths = 4*floor(log(signal_length)))$intervals) > 0)
    zz[4,jj] <- zz[4,jj] + as.numeric(nrow(mosum_ints(xx, alpha = coverage_levels[jj])$intervals) > 0)
    zz[5,jj] <- zz[5,jj] + as.numeric(nrow(smuce_ints(xx, coverage_levels[jj])$intervals) > 0)
    zz[6,jj] <- zz[6,jj] + as.numeric(nrow(cpt_diff_inf(xx, 0, coverage_levels[jj], aa = 2, LL = 0)$intervals) > 0)
    zz[7,jj] <- zz[7,jj] + as.numeric(nrow(cpt_diff_inf(xx, 0, coverage_levels[jj], aa = 2, LL = 5)$intervals) > 0)
    zz[8,jj] <- zz[8,jj] + as.numeric(nrow(cpt_diff_inf(xx, 0, coverage_levels[jj], aa = sqrt(2), LL = 0)$intervals) > 0)
    zz[9,jj] <- zz[9,jj] + as.numeric(nrow(cpt_diff_inf(xx, 0, coverage_levels[jj], aa = sqrt(2), LL = 5)$intervals) > 0)
  }
  
  setTxtProgressBar(pb, ii)
}

deg_0_coverage <- zz / N_draws

colnames(deg_0_coverage) <- coverage_levels

save(deg_0_coverage, file = "deg_0_coverage.RData")



## Degree 1 polynomials
##

print("Degree 1 polynomials")

zz <- matrix(0, 6, length(coverage_levels))

for (ii in 1:N_draws)
{
  xx <- rnorm(signal_length)
  
  for (jj in seq_along(coverage_levels))
  {
    # zz[1,jj] <- zz[1,jj] + as.numeric(nrow(nsp_poly(xx, deg = 1, alpha = coverage_levels[jj])$intervals) > 0)
    # zz[2,jj] <- zz[2,jj] + as.numeric(nrow(nsp_poly(xx, deg = 1 ,thresh.val = quantile(nsp_sims_deg1, 1 - coverage_levels[jj]))$intervals) > 0)
    zz[3,jj] <- zz[3,jj] + as.numeric(nrow(cpt_diff_inf(xx, 1, coverage_levels[jj], aa = 2, LL = 0)$intervals) > 0)
    zz[4,jj] <- zz[4,jj] + as.numeric(nrow(cpt_diff_inf(xx, 1, coverage_levels[jj], aa = 2, LL = 5)$intervals) > 0)
    zz[5,jj] <- zz[5,jj] + as.numeric(nrow(cpt_diff_inf(xx, 1, coverage_levels[jj], aa = sqrt(2), LL = 0)$intervals) > 0)
    zz[6,jj] <- zz[6,jj] + as.numeric(nrow(cpt_diff_inf(xx, 1, coverage_levels[jj], aa = sqrt(2), LL = 5)$intervals) > 0)
  }
  
  setTxtProgressBar(pb, ii)
}

deg_1_coverage <- zz / N_draws

colnames(deg_1_coverage) <- coverage_levels

save(deg_1_coverage, file = "deg_1_coverage.RData")



## Degree 2 polynomials
##

print("Degree 2 polynomials")

zz <- matrix(0, 6, length(coverage_levels))

for (ii in 1:N_draws)
{
  xx <- rnorm(signal_length)
  
  for (jj in seq_along(coverage_levels))
  {
    # zz[1,jj] <- zz[1,jj] + as.numeric(nrow(nsp_poly(xx, deg = 2, alpha = coverage_levels[jj])$intervals) > 0)
    # zz[2,jj] <- zz[2,jj] + as.numeric(nrow(nsp_poly(xx, deg = 2 ,thresh.val = quantile(nsp_sims_deg2, 1 - coverage_levels[jj]))$intervals) > 0)
    zz[3,jj] <- zz[3,jj] + as.numeric(nrow(cpt_diff_inf(xx, 2, coverage_levels[jj], aa = 2, LL = 0)$intervals) > 0)
    zz[4,jj] <- zz[4,jj] + as.numeric(nrow(cpt_diff_inf(xx, 2, coverage_levels[jj], aa = 2, LL = 5)$intervals) > 0)
    zz[5,jj] <- zz[5,jj] + as.numeric(nrow(cpt_diff_inf(xx, 2, coverage_levels[jj], aa = sqrt(2), LL = 0)$intervals) > 0)
    zz[6,jj] <- zz[6,jj] + as.numeric(nrow(cpt_diff_inf(xx, 2, coverage_levels[jj], aa = sqrt(2), LL = 5)$intervals) > 0)
  }
  
  setTxtProgressBar(pb, ii)
}

deg_2_coverage <- zz / N_draws

colnames(deg_2_coverage) <- coverage_levels

save(deg_2_coverage, file = "deg_1_coverage.RData")
