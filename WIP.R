
set.seed(42)

library(ChangePointInference)



## Blocks
## 

png("plots/blocks-signal.png", width = 600, height = 600)

blocks_signal <- c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40))

yy <- blocks_signal + rnorm(length(blocks_signal), sd = 5)

dev.off()

yy |> plot(type = "l", col = "grey")

blocks_signal |> lines(lty = 2, lwd = 3)



diff_out <- diffInf(yy, 0)

diff_out |> plot(type = "l", col = "grey")

blocks_signal |> lines(lty = 2, lwd = 3)


## Waves
##


## Hills
##

hills_signal <- 10 * rep((1:100/100) * (1 - (1:100/100)), 4)

yy <- hills_signal + rnorm(length(hills_signal), sd = .25)


yy |> plot(type = "l", col = "grey")

hills_signal |> lines(lty = 2, lwd = 3)



diff_out <- diffInf(yy, 2)

diff_out |> plot(type = "l", col = "grey")

hills_signal |> lines(lty = 2, lwd = 3)


