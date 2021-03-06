## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE,
  fig.width=7, 
  fig.height=5.5
)

## ----setup, warning=FALSE, message=FALSE---------------------------------
library(convexhull)
library(scales)

## ------------------------------------------------------------------------
set.seed(2)
x <- rnorm(100)
y <- rnorm(100)

## ------------------------------------------------------------------------
op <- par(mfrow=c(2,3), mar=c(3,2,3,1))

plot(x,y, main="convex hull")
convex_hull(x,y)

plot(x,y, main="colored hull")
convex_hull(x,y, col="#FF000050", density=10, angle=45, border="red")

plot(x,y, main="20% of most distant\npoints removed")
convex_hull(x,y, alpha=.2, col="#FF000050", density=10, angle=45, border="red")

plot(x,y, main="initial hull points\nremoved (peeled off)")
convex_hull(x,y, peel=TRUE, col="#00FF0050", border="red")

plot(x,y, main="smoothed hull\n(spline interpolated)")
convex_hull(x,y, smooth=2, col="#FF000050", border="red")

plot(x,y, main="another smooting\nmethod")
convex_hull(x,y, peel=TRUE, smooth=1, shape=1, col="#00FF0050", border="red")

## ------------------------------------------------------------------------
library(MASS)

# example 1
n <- 100
set.seed(0)
x <- rnorm(n, 3) 
y <- rnorm(n, 3)
plot(x,y)
den <- kde2d(x, y, n=100)     # estimate non-parameteric density surface via kernel smoothing

cols = alpha(heat.colors(12), .1)
crit_dens_image(den, col=cols, prob=.1)
crit_contour(den, n=100, col="red", drawlabels=FALSE, prob=.1)

## ------------------------------------------------------------------------
# example 2
n <- 200
set.seed(0)
x <- rnorm(n, c(1,3)) 
y <- rnorm(n, c(1,3))
d <- data.frame(x,y, g=1:2)
plot(x,y, pch=16, cex=.7, col=1:2)
g1 <- subset(d, g==1)
g2 <- subset(d, g==2)
den1 <- kde2d(g1$x, g1$y, n=100)     # estimate non-parameteric density surface via kernel smoothing
den2 <- kde2d(g2$x, g2$y, n=100)     # estimate non-parameteric density surface via kernel smoothing

reds <- alpha("red", seq(.05, 1, len=20))
blues <- alpha("blue", seq(.05, 1, len=20))
# cols <- colorRampPalette(c("white", "darkred"))(10)
# cols <- set_alpha_color_value(cols, alpha=.5)
crit_dens_image(den1, col=blues, prob=.3)
crit_contour(den1, n=100, lty=3, drawlabels=FALSE, prob=.3)
crit_dens_image(den2, col=reds, prob=.3)
crit_contour(den2, n=100, lty=3, drawlabels=FALSE, prob=.3)

