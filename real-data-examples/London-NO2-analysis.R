
library(ChangePointInference)

library(DeCAFS)

library(AR1seg)

source("wcm-gsa.R")

source("estimate-signal-from-breaks.R")


#---------------------------------------
#
# Analysis following Cho & Fryzlewicz 
#
#---------------------------------------

no <- read.csv("AirQualityData.csv", header = TRUE) 

hols <- read.csv("ukbankholidays.csv", header = TRUE)

allhols <- as.Date(hols$UK.BANK.HOLIDAYS, "%d-%b-%Y")

day <- rep(1:7, ceiling(length(no$NO2)/7))[1:length(no$NO2)]

no$day <- day

no <- no[!is.na(no$NO2), ]

dates <- as.Date(no$Date)

xx <- sqrt(no$NO2)

dat <- data.frame(
  months = format.Date(dates, "%m"), 
  days = factor(no$day), 
  is.hol = dates %in% allhols, 
  y = xx
  )

fit <- lm(y ~ ., data = dat[format.Date(dates, "%Y") < "2011" & format.Date(dates, "%Y") > "2003", ])

x <- xx - predict(fit, newdata = dat)


#---------------------------------------
#
# Plot raw data
#
#---------------------------------------

png("../plots/London-NO2.png", width = 600, height = 300)

plot(
  x, 
  type = "l", 
  xlab = "", 
  ylab = "", 
  xaxt = "n"
)

axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])

rect(
  xleft = which(dates == "2003-02-01"),
  ybottom = -100,
  xright =  which(dates == "2003-02-26"),
  ytop = 100,
  col = "red",
  border = NA
)

abline(v = which(dates == "2019-04-08"), lty = 2, col = "red")

abline(v = which(dates == "2020-03-23"), lty = 2, col = "red")

dev.off()


#---------------------------------------
#
# DIF2-TAVC 
#
#---------------------------------------

diff_inf_obj <- diffInf(x, degree = 0, independent_noise = FALSE)


#---------------------------------------
#
# AR1seg
#
#---------------------------------------

ar1_seg_obj <- AR1seg_func(x)

png("../plots/London-NO2-AR1-seg.png", width = 600, height = 300)

plot(
  diff_inf_obj, 
  type = "l", 
  col = "grey",
  xlab = "", 
  ylab = "", 
  xaxt = "n")

axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])

lines(estimate_signal_from_breaks(x, ar1_seg_obj$PPSelectedBreaks), col = "blue", lwd = 1.5)

abline(v = ar1_seg_obj$PPSelectedBreaks[-length(ar1_seg_obj$PPSelectedBreaks)], lty = 2, col = "blue")

dev.off()

## DeCAFS

p <- estimateParameters(x, model = "AR")

decafs_obj <- DeCAFS(x, modelParam = p, beta = 2 * log(length(x)))

png("../plots/London-NO2-decafs.png", width = 600, height = 300)

plot(
  diff_inf_obj, 
  type = "l", 
  col = "grey",
  xlab = "", 
  ylab = "", 
  xaxt = "n")

axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])

lines(estimate_signal_from_breaks(x, decafs_obj$changepoints), col = "blue", lwd = 1.5)

abline(v = decafs_obj$changepoints, lty = 2, col = "blue")

dev.off()


#---------------------------------------
#
# WGS
#
#---------------------------------------

wgs_obj <- wcm.gsa(x, max.iter = 10, min.len = 10)

png("../plots/London-NO2-wgs.png", width = 600, height = 300)

plot(
  diff_inf_obj, 
  cpt_loc_est = NULL, 
  type = "l", 
  col = "grey",
  xlab = "", 
  ylab = "", 
  xaxt = "n")

axis(1, at = seq(1, length(x), length.out = 20), label = format.Date(dates, "%Y-%m")[seq(1, length(x), length.out = 20)])

lines(estimate_signal_from_breaks(x, wgs_obj$cp), col = "blue", lwd = 1.5)

abline(v = wgs_obj$cp, lty = 2, col = "blue")

dev.off()



#---------------------------------------
#
# Verify claims in Section 5.2
#
#---------------------------------------

## WCM.gSa identifies only one cpt in third interval
## 

y <- x[6083:6806]

wcm.gsa(y)$cp

## WCM.gSa and DeCAFS find no cpts between first and second interval
##

y <- x[1007:2478]

p <- estimateParameters(y, model = "AR")

DeCAFS(y, modelParam = p, beta = 2 * log(length(x)))$changepoints

wcm.gsa(y)$cp

