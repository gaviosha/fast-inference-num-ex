
#---------------------------------------
#
# Setup
#
#---------------------------------------

set.seed(42)

library(ChangePointInference)

nn <- 500       

cpt <- nn/2

t1 <- 1:cpt; t2 <- (cpt+1):nn


#---------------------------------------
#
# Change in d0 and d1
#
#---------------------------------------

a0 <- 1/2; b0 <- -1/2; a1 <- -2; b1 <- 2

ff <- c(
  a0 + a1 * (t1/nn - cpt/nn),
  b0 + b1 * (t2/nn - cpt/nn)
)

yy <- ff + rnorm(length(ff), sd = 0.5)


png("../plots/d0-d1-signal.png", width = 600, height = 600)

plot(
  yy, 
  type = "l", 
  col = "grey", 
  lwd = 3, 
  xlab = "", 
  ylab = ""
)

lines(
  ff, 
  lty = 2, 
  lwd = 3
)

dev.off()

png("../plots/d0-d1-ints.png", width = 600, height = 600)

plot(diffInf(yy, 1), type = "l", col = "grey", lwd = 3, cpt_loc_est = TRUE, xlab = "", ylab = "")

dev.off()


#---------------------------------------
#
# Change in d1 and d2
#
#---------------------------------------

a1 <- 6; b1 <- -6; a2 <- -12; b2 <- 12

ff <- c(
  a1 * (t1/nn - cpt/nn) + a2 * (t1/nn - cpt/nn)^2,
  b1 * (t2/nn - cpt/nn) + b2 * (t2/nn - cpt/nn)^2
)

yy <- ff + rnorm(length(ff), sd = 0.5)

png("../plots/d1-d2-signal.png", width = 600, height = 600)

plot(
  yy, 
  type = "l", 
  col = "grey", 
  lwd = 3, 
  xlab = "", 
  ylab = ""
)

lines(
  ff, 
  lty = 2, 
  lwd = 3
)

dev.off()

png("../plots/d1-d2-ints.png", width = 600, height = 600)

plot(diffInf(yy, 2), type = "l", col = "grey", lwd = 3, cpt_loc_est = TRUE, xlab = "", ylab = "")

dev.off()


#---------------------------------------
#
# Change in d0 and d2
#
#---------------------------------------

a0 <- 1/2; b0 <- -1/2; a2 <- -2; b2 <- 2

ff <- c(
  a0 + a2 * (t1/nn - cpt/nn)^2,
  b0 + b2 * (t2/nn - cpt/nn)^2
)

yy <- ff + rnorm(length(ff), sd = 0.5)

png("../plots/d0-d2-signal.png", width = 600, height = 600)

plot(
  yy, 
  type = "l", 
  col = "grey", 
  lwd = 3, 
  xlab = "", 
  ylab = ""
)

lines(
  ff, 
  lty = 2, 
  lwd = 3
)

dev.off()

png("../plots/d0-d2-ints.png", width = 600, height = 600)

plot(diffInf(yy, 2), type = "l", col = "grey", lwd = 3, cpt_loc_est = TRUE, xlab = "", ylab = "")

dev.off()
