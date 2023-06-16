
library(ChangePointInference)

library(not)

library(IDetect)

library(strucchange)

library(cpop)

spbdm <- read.csv("spnbmd.csv")

spbdm <- aggregate(spnbmd ~ age + sex, data = spbdm, FUN = mean)

set.seed(42)


#--------------------------------------
#
# Analysis of male bone density
#
#--------------------------------------


## Plot raw data
##

mal_spbdm <- spbdm[spbdm$sex == "mal", 3]

png("../plots/male-spbdm.png", width = 600, height = 600)

plot(mal_spbdm, 
     type = "l", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n"
     )

axis(1, at = seq(1, length(mal_spbdm), length.out = 17), label = 9:25)

dev.off()


## Change point analysis
##

png("../plots/male-spbdm-cpt.png", width = 600, height = 600)

diff_obj <- diffInf(mal_spbdm,1, gaussian_noise = FALSE, min_scale = 2)

plot(diff_obj,
     type = "l", 
     col = "grey", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n"
     )

axis(1, at = seq(1, length(mal_spbdm), length.out = 17), label = 9:25)

abline(v = features(not(mal_spbdm, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(mal_spbdm, contrast = "pcwsLinContMean"))$cpt, lty = 3, lwd = 3, col = "blue")

abline(v = ID(mal_spbdm, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")

abline(v = breakpoints(mal_spbdm ~ poly(seq_along(mal_spbdm),1))$breakpoints, lty = 2, lwd = 3, col = "green4")

abline(v = changepoints(cpop(mal_spbdm))[,1], lty = 2, lwd = 3, col = "purple")

dev.off()


## Change point analysis on a sub-interval
## 

png("../plots/male-spbdm-cpt-subint.png", width = 600, height = 600)

ss <- diff_obj$intervals[1,1]

ee <- diff_obj$intervals[1,2]

mal_spbdm_subint <- mal_spbdm[ss:ee]

plot(mal_spbdm_subint,
     type = "l", 
     col = "grey", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n")

axis(1, at = seq(1, length(mal_spbdm_subint), length.out = 10), label = 12:21)

abline(v = features(not(mal_spbdm_subint, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(mal_spbdm_subint, contrast = "pcwsLinContMean"))$cpt, lty = 3, lwd = 3, col = "blue")

abline(v = ID(mal_spbdm_subint, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")

abline(v = breakpoints(mal_spbdm_subint ~ poly(seq_along(mal_spbdm_subint),1))$breakpoints, lty = 2, lwd = 3, col = "green4")

abline(v = changepoints(cpop(mal_spbdm_subint))[,1], lty = 2, lwd = 3, col = "purple")

dev.off()


#--------------------------------------
#
# Analysis of female bone density
#
#--------------------------------------


## Plot raw data
##

fem_spbdm <- spbdm[spbdm$sex == "fem", 3]

png("../plots/female-spbdm.png", width = 600, height = 600)

plot(fem_spbdm,
     type = "l", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n"
     )

axis(1, at = seq(1, length(fem_spbdm), length.out = 17), label = 9:25)

dev.off()


## Change point analysis
##

png("../plots/female-spbdm-cpt.png", width = 600, height = 600)

diff_obj <- diffInf(fem_spbdm,1)

plot(diff_obj, 
     type = "l", 
     col = "grey", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n"
)

axis(1, at = seq(1, length(fem_spbdm), length.out = 17), label = 9:25)

abline(v = features(not(fem_spbdm, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(fem_spbdm, contrast = "pcwsLinContMean"))$cpt, lty = 2, lwd = 3, col = "green")

abline(v = ID(fem_spbdm, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")

abline(v = breakpoints(fem_spbdm ~ poly(seq_along(fem_spbdm),1))$breakpoints, lty = 2, lwd = 3, col = "green4")

abline(v = changepoints(cpop(fem_spbdm))[,1], lty = 2, lwd = 3, col = "purple")

dev.off()


## Change point analysis on a sub-interval
##

png("../plots/fem-spbdm-cpt-subint.png", width = 600, height = 600)

ss <- diff_obj$intervals[1,1]

ee <- diff_obj$intervals[1,2]

fem_spbdm_subint <- fem_spbdm[ss:ee]

plot(fem_spbdm_subint,
     type = "l", 
     col = "grey", 
     lwd = 3, 
     xlab = "age", 
     ylab = "bone denisty", 
     xaxt = "n")

axis(1, at = seq(1, length(fem_spbdm_subint), length.out = 10), label = 9:18)

abline(v = features(not(fem_spbdm_subint, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(fem_spbdm_subint, contrast = "pcwsLinContMean"))$cpt, lty = 2, lwd = 3, col = "green")

abline(v = ID(fem_spbdm_subint, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")

abline(v = breakpoints(fem_spbdm_subint ~ poly(seq_along(fem_spbdm_subint),1))$breakpoints, lty = 2, lwd = 3, col = "green4")

abline(v = changepoints(cpop(fem_spbdm_subint))[,1], lty = 2, lwd = 3, col = "purple")

dev.off()
