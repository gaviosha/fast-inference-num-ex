#### WIP 


yy <- blocks_signal + noise_types$N4(length(blocks_signal),5)

yy |> plot(type = "l", col = "grey")

blocks_signal |> lines(lty = 2, lwd = 2, col = "red")




yy <- waves_signal + rnorm(length(waves_signal), sd = 10)

yy |> plot(type = "l", col = "grey")

waves_signal |> lines(lty = 2, lwd = 2, col = "red")



yy <- hills_signal + rnorm(length(hills_signal), sd = 2)

yy |> plot(type = "l", col = "grey")

hills_signal |> lines(lty = 2, lwd = 2, col = "red")

# 
# library(StepRNonparam)
# 
# out <- StepRNonparam::stepFitNonparam(yy, sd = 1, alpha = 0.1)
# 
# cbind(out$leftEnd[-1], out$rightEnd[-length(out$rightEnd)])
# 
# 
# dep_smuce_ints(yy)
# 
# cbind(foo$leftEnd[-1], foo$rightEnd[-length(foo$rightEnd)])
