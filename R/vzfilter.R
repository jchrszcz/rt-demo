## ======================================================================
## Copyright 2007--2010, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## The use of this software is permitted for academic purposes only,
## provided its use is acknowledged.  Commercial use is not allowed
## without the permission of the author.  This software is made available
## AS IS, and no warranty -- about the software, its performance, or its
## conformity to any specification -- is given or implied.
##
## Supported by the National Science Foundation under award numbers
## BCS-0738059, DMS-0604963, DMS-0605052, SES-0214574 and SES-0437251
## ======================================================================
## ======================================================================
## File     : detrending_RT.R
## Contains : R code to detrend RT sequences.
## Version  : 0.1
## Updated  : pfc@stat.osu.edu, Sep 2008.
##
## Reference:
##
## P. F. Craigmile, M. Peruggia, T. Van Zandt,
## "Several Methods for Detrending Response Time Series"
## ======================================================================

transform.to.normal <- function (x) {
## ======================================================================
## transform the time series 'x' to normal scores.
## ======================================================================  
  n <- length(x)
  qnorm((rank(x)-0.5)/n)
}

transform.back <- function (x.norm.fitted, x,
  x.norm=transform.to.normal(x)) {
## ======================================================================
## Transform back from normal scores.
##
## Arguments:
## 'x.norm.fitted' : the estimated trend on the normal score scale.
## 'x'             : the original time series on the original scale
##                   before detrending was carried out.
## 'x.norm'        : 'x' transformed to normal scores
##
## Returns the estimated trend on the original data scale.
## ======================================================================  
  reverse.quantile <- function (p, p1, p2, x1, x2) {
    (p - p1) / (p2 - p1) * (x2 - x1) + x1
  }
  sorted.x <- sort(x)
  sorted.x.norm <- sort(x.norm)
  x.fitted <- NULL
  for (i in 1:length(x)) {
    j <- sum(sorted.x.norm <= x.norm.fitted[i])
    x.fitted[i] <- reverse.quantile(pnorm(x.norm.fitted[i]),
      pnorm(sorted.x.norm[j]), pnorm(sorted.x.norm[j+1]),
      sorted.x[j], sorted.x[j+1])
  }
  x.fitted
}

est.poly.trend <- function (x, deg=2) {
## ======================================================================
## estimate a polynomial trend of degree 'deg' (the default is 2,
## a quadratic trend) for the time series 'x'.
## ======================================================================
require(signal)
  ts <- seq(length(x))/length(x)
  fitted(lm(x ~ stats::poly(ts, deg)))
}

est.loess <- function (x, span=0.75) {
## ======================================================================
## estimate the trend from a time series 'x' using a loess fit with
## a given 'span' (default span is to use 0.75 of the data in a window).
## ======================================================================
  ts <- seq(length(x))/length(x)
  fitted(loess(x ~ ts, span=span))
}

est.MA.filter <- function (x, wlen=32) {
## ======================================================================
## estimate the trend from a time series 'x' by using a moving average
## smoothing window of length 'wlen' (default value for wlen is 32).
## The boundary problem is circumvented by assuming that the
## time series 'x' is periodic, which may not be true in practice.
## ======================================================================
  as.numeric(stats::filter(x, rep(1,wlen)/wlen, circular=T))
}

est.smoothing.splines <- function (x, AR1.errors=TRUE) {
## ======================================================================
## estimate the trend from a time series 'x' by using a smoothing spline.
## If 'AR1.errors' is TRUE use an AR(1) error structure, otherwise
## assume the errors are uncorrelated.
##
## Requires the 'assist' R library to be installed.
##
## NOTE: Because of an error in the assist library, this function
##       creates and deletes two global variables called 'est.ss.x'
##       and 'est.ss.ts'.
## ======================================================================
  est.ss.x  <<- x
  est.ss.ts <<- seq(length(x))/length(x)
  require(assist)
  if (AR1.errors) {
    fit <- ssr(est.ss.x ~ est.ss.ts, rk=cubic(est.ss.ts),
    corr=corAR1(), spar="m")$fit
  }
  else {
    fit <- ssr(est.ss.x ~ est.ss.ts, rk=cubic(est.ss.ts), spar="m")$fit
  }
  rm(est.ss.x, est.ss.ts, envir = globalenv())
  fit
}

est.norm.smoothing.splines <- function (x, AR1.errors=TRUE) {
## ======================================================================
## estimate the trend from a time series 'x' after transforming to
## normal scores by using a smoothing spline.  The estimated trend is
## transformed back to the original scale.
## If 'AR1.errors' is TRUE use an AR(1) error structure, otherwise  
## assume the errors are uncorrelated.
##
## Requires the 'assist' R library to be installed.
##
## NOTE: Because of an error in the assist library, this function
##       creates and deletes two global variables called
##       'est.norm.ss.x' and 'est.norm.ss.ts'.
## ======================================================================
est.norm.ss.x.norm <<- transform.to.normal(x)
est.norm.ss.ts     <<- seq(length(x))/length(x)
require(assist)
if (AR1.errors) {
  fitted.norm <- ssr(est.norm.ss.x.norm ~ est.norm.ss.ts,
    rk=cubic(est.norm.ss.ts),
    corr=corAR1(), spar="m")$fit
  }
else {
  fitted.norm <- ssr(est.norm.ss.x.norm ~ est.norm.ss.ts,
    rk=cubic(est.norm.ss.ts), spar="m")$fit
  }
  fit <- transform.back(fitted.norm, x, est.norm.ss.x.norm)
  rm(est.norm.ss.x.norm, est.norm.ss.ts, envir = globalenv())
  fit  
}