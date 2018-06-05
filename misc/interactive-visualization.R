

library(manipulate) # RStudio is required

# Beta

x <- seq(0, 1, by=0.01)
manipulate(
  { ci <- qprop(c(0.05, .5, 0.95), size, mean)
  plot(x, dprop(x, size, mean),
       col="blue", lwd=2, type="l", las=1, bty="n",
       ylab="density", xlab="", ylim=c(0, 12),
       main="Beta distribution")
  box()
  mtext(paste(c("95% CI:", round(ci, 2)),
              collapse="  "), cex=0.8, side=3)
  if (printci) abline(v=ci, lty=c(3,2,3))
  },
  size=slider(0.001, 100, step=0.001, initial=100),
  mean=slider(0, 1, step=0.001, initial=0.5),
  printci=checkbox(TRUE, "Show 95% CI"))


# GEV

x <- seq(-5, 5, by=0.01)
manipulate(
  { ci <- qgev(c(0.05, .5, 0.95), mu, sigma, xi)
  plot(x, dgev(x, mu, sigma, xi),
       col="blue", lwd=2, type="l", las=1, bty="n",
       ylab="density", xlab="", ylim=c(0, 2),
       main="GEV distribution")
  box()
  mtext(paste(c("95% CI:", round(ci, 2)),
              collapse="  "), cex=0.8, side=3)
  if (printci) abline(v=ci, lty=c(3,2,3))
  },
  mu=slider(-2, 2, step=0.001, initial=0),
  sigma=slider(0.001, 5, step=0.001, initial=1),
  xi=slider(-2, 2, step=0.001, initial=0),
  printci=checkbox(TRUE, "Show 95% CI"))


x <- seq(-5, 5, by=0.01)
manipulate(
  { ci <- qgev(c(0.05, .5, 0.95), mu, sigma, xi)
  plot(x, pgev(x, mu, sigma, xi),
       col="blue", lwd=2, type="l", las=1, bty="n",
       ylab="density", xlab="", ylim=c(0, 1),
       main="GEV distribution")
  box()
  mtext(paste(c("95% CI:", round(ci, 2)),
              collapse="  "), cex=0.8, side=3)
  if (printci) abline(v=ci, lty=c(3,2,3))
  },
  mu=slider(-2, 2, step=0.001, initial=0),
  sigma=slider(0.001, 5, step=0.001, initial=1),
  xi=slider(-2, 2, step=0.001, initial=0),
  printci=checkbox(TRUE, "Show 95% CI"))




# GPD

x <- seq(-5, 5, by=0.01)
manipulate(
  { ci <- qgpd(c(0.05, .5, 0.95), mu, sigma, xi)
  plot(x, dgpd(x, mu, sigma, xi),
       col="blue", lwd=2, type="l", las=1, bty="n",
       ylab="density", xlab="", ylim=c(0, 2),
       main="GPD distribution")
  box()
  mtext(paste(c("95% CI:", round(ci, 2)),
              collapse="  "), cex=0.8, side=3)
  if (printci) abline(v=ci, lty=c(3,2,3))
  },
  mu=slider(-2, 2, step=0.001, initial=0),
  sigma=slider(0.001, 5, step=0.001, initial=1),
  xi=slider(-2, 2, step=0.001, initial=0),
  printci=checkbox(TRUE, "Show 95% CI"))



x <- seq(-5, 5, by=0.01)
manipulate(
  { ci <- qgpd(c(0.05, .5, 0.95), mu, sigma, xi)
  plot(x, pgpd(x, mu, sigma, xi),
       col="blue", lwd=2, type="l", las=1, bty="n",
       ylab="density", xlab="", ylim=c(0, 1),
       main="GPD distribution")
  box()
  mtext(paste(c("95% CI:", round(ci, 2)),
              collapse="  "), cex=0.8, side=3)
  if (printci) abline(v=ci, lty=c(3,2,3))
  },
  mu=slider(-2, 2, step=0.001, initial=0),
  sigma=slider(0.001, 5, step=0.001, initial=1),
  xi=slider(-2, 2, step=0.001, initial=0),
  printci=checkbox(TRUE, "Show 95% CI"))



# Inverse Gamma

manipulate({
  q <- qinvgamma(c(0.05, .5, 0.95, 0.99), alpha, beta)
  x <- seq(0, q[4], length.out=5000)
  u <- rinvgamma(50000, alpha, beta)
  d <- dinvgamma(x, alpha, beta)
  mxd <- max(d)
  hist(u, 1000, freq=FALSE, col='darkgrey', border='darkgrey',
       xlim=c(0, q[4]), ylim=c(0, mxd+0.1*mxd),
       main="Inverse Gamma distribution", las=1, bty="n")
  lines(x, d, col='blue', lwd=2, type="l")
  box()
  mtext(paste(c("95% CI:", round(q[1:3], 2)),
              collapse="  "), cex=0.8, side=3)
  },
  alpha=slider(0.01, 10, step=0.001, initial=1),
  beta=slider(0.01, 5, step=0.001, initial=1))
