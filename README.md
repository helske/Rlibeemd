[![Build Status](https://travis-ci.org/helske/Rlibeemd.png?branch=master)](https://travis-ci.org/helske/Rlibeemd)
[![codecov.io](http://codecov.io/github/helske/Rlibeemd/coverage.svg?branch=master)](http://codecov.io/github/helske/Rlibeemd?branch=master)
[![downloads](http://cranlogs.r-pkg.org/badges/Rlibeemd)](http://cranlogs.r-pkg.org/badges/Rlibeemd)
[![cran version](http://www.r-pkg.org/badges/version/Rlibeemd)](http://cran.r-project.org/package=Rlibeemd)

# Rlibeemd #

An R interface for [libeemd C library](https://bitbucket.org/luukko/libeemd) for ensemble empirical mode decomposition (EEMD) and its complete variant (CEEMDAN). These methods decompose possibly nonlinear and/or nonstationary time series data into a finite amount of components (called IMFs, insintric mode functions) separated by instantaneous frequencies. This decomposition provides a powerful method to look into the different processes behind a given time series, and provides a way to separate short time-scale events from a general trend.

### Example ###
Here a CEEMDAN decomposition is performed for the UK gas consumption series (length n = 108). 
By default, `ceemdan` extracts [log_2(n)] components, so here we get five IMFs and the residual.

```{r, fig.height = 4, fig.width = 8}
library("Rlibeemd")
data(UKgas, package = "datasets")
imfs <- ceemdan(UKgas, ensemble_size = 1000)
plot(imfs, main = "Five IMFs and residual extracted by CEEMDAN algorithm")
```
![imfs](https://github.com/helske/Rlibeemd/blob/master/imfs.png)

The residual components shows clear trend whereas the first IMF see so contain clear multiplicative trend. The remaining IMFs are bit more complex, and one could argue that they are partly seasonal, trend or just some irregularity i.e. noise. 

Lets compare the decomposition with basic structural time series model fit from `StructTS` (for smoothing of more complex state space models, one could use [KFAS](https://github.com/helske/KFAS))

```{r, fig.height = 4, fig.width = 8}
bsm <- tsSmooth(StructTS(UKgas))
plot(bsm[, c(1, 3)], main = "Local linear trend and seasonal components by StructTS")
```
![bsm](https://github.com/helske/Rlibeemd/blob/master/bsm.png)

``StructTS`` decomposes the data for three components, where one of the components is (possibly time varying) slope, which has no direct effect to overall signal (it is the slope of the level component).

```{r, fig.height=4, fig.width=8}
ts.plot(cbind(UKgas, imfs[, ncol(imfs)], rowSums(imfs[, 5:6]), bsm[,"level"]), col = 1:4,
  main = "Quarterly UK gas consumption", ylab = "Million therms")
legend("topleft", c("Observations", "Residual", "Last IMF + residual", "Trend from BSM"),
  col = 1:4, lty = 1)
```
![ceemdan_and_bsm](https://github.com/helske/Rlibeemd/blob/master/ceemdan_and_bsm.png)

The IMF_5 + residual is quite close to the trend obtained by structural time series model of `StructTS`.




### Installing libeemd ###

Package is available at [CRAN](http://cran.r-project.org/web/packages/Rlibeemd/index.html).

If you want to compile Rlibeemd from source, you will also need GNU Scientific Library.

You can also install the latest development version from the github using the devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/Rlibeemd")
```
