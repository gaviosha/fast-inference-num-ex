
library(ChangePointInference)

library(not)

library(IDetect)

library(strucchange)

spbdm <- read.csv("data/spnbmd.csv")

spbdm <- aggregate(spnbmd ~ age + sex, data = spbdm, FUN = mean)

set.seed(42)


## Analysis of male bone denisty
##

mal_spbdm <- spbdm[spbdm$sex == "mal", 3]

diff_obj <- diffInf(mal_spbdm,1)

plot(diff_obj, type = "l", col = "grey", lwd = 3, xlab = "", ylab = "bone denisty")


abline(v = features(not(mal_spbdm, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(mal_spbdm, contrast = "pcwsLinContMean"))$cpt, lty = 3, lwd = 3, col = "blue")

abline(v = ID(mal_spbdm, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")

abline(v = breakpoints(mal_spbdm ~ poly(seq_along(mal_spbdm),1))$breakpoints, lty = 2, lwd = 3, col = "green4")

abline(v = changepoints(cpop(mal_spbdm)), lty = 2, lwd = 3, col = "purple")


## Analysis of female bone denisty
##

fem_spbdm <- spbdm[spbdm$sex == "fem", 3]

diff_obj <- diffInf(fem_spbdm,1)

plot(diff_obj, type = "l", col = "grey", lwd = 3)


abline(v = features(not(fem_spbdm, contrast = "pcwsLinMean"))$cpt, lty = 2, lwd = 3, col = "blue")

abline(v = features(not(fem_spbdm, contrast = "pcwsLinContMean"))$cpt, lty = 2, lwd = 3, col = "green")

abline(v = ID(fem_spbdm, contrast = "slope")$cpt, lty = 2, lwd = 3, col = "orange")


