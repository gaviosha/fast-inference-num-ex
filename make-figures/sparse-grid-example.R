
## params
##

nn <- 20

aa <- 2

eps <- 0.02

n_scales <- floor(logb(nn/2,aa))


## base plot
##

png("../plots/sparse-grid.png", width = 600, height = 300)

plot(NULL, xlim = c(0,nn+1), ylim = c(0,1), axes = FALSE, ann = FALSE)

axis(1)


## plot scales
##

for (kk in 1:n_scales)
{
  ww <- floor(aa**kk)
  
  yy <- seq(kk/n_scales,(kk-1)/n_scales, length.out = nn - ww)
  
  axis(2, at = c((kk-1)/n_scales + eps, kk/n_scales - eps), ann = FALSE, labels = FALSE)
  
  for (ll in 1:(nn - ww))
  {
    segments(x0 = ll, x1 = ll+ww, y0 = yy[ll], y1 = yy[ll], lwd = 3)
  }
}

dev.off()