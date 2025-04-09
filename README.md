[![R-CMD-check](https://github.com/helske/Rlibeemd/workflows/R-CMD-check/badge.svg)](https://github.com/helske/Rlibeemd/actions)
[![codecov.io](https://codecov.io/github/helske/Rlibeemd/coverage.svg?branch=main)](https://app.codecov.io/github/helske/Rlibeemd?branch=main)
[![downloads](http://cranlogs.r-pkg.org:443/badges/Rlibeemd](https//cranlogs.r-pkg.org:443/badges/Rlibeemd)
[![cran version](https://www.r-pkg.org/badges/version/Rlibeemd)](https://CRAN.R-project.org/package=Rlibeemd)

# Rlibeemd #

An R interface for [libeemd C library](https://bitbucket.org/luukko/libeemd) for 
ensemble empirical mode decomposition (EEMD) and its complete variant (CEEMDAN). 
These methods decompose possibly nonlinear and/or nonstationary time series data 
into a finite amount of components (called IMFs, insintric mode functions) 
separated by instantaneous frequencies. This decomposition provides a 
powerful method to look into the different processes behind a given time 
series, and provides a way to separate short time-scale events from a general trend.

If you use Rlibeemd/libeemd for scientific work please cite 
[*Luukko, P.J.J., Helske, J., Räsänen, E., Comput. Stat. **31**, 545 (2016)*](https://dx.doi.org/10.1007/s00180-015-0603-9) ([also on arXiv](https://arxiv.org/abs/1707.00487)). 
This article also describes in detail what `libeemd` actually computes. 
You should definitely read it if you are unsure about what EMD, EEMD and CEEMDAN are.
 
# OpenMP parallel computing support

Current CRAN policies do not allow the use of `SHLIB_OPENMP_CFLAGS` combined
with linking with C++. Therefore the CRAN version does not use OpenMP at all 
anymore (OpenMP flags have been removed from `Makevars`), but the the version 
on GitHub version does. So if you want to use parallel version of the `Rlibeemd`, please install the package via

```
install.packages('Rlibeemd', repos = 'https://helske.r-universe.dev')
```


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

The residual components shows smooth trend whereas the first IMF contains clear multiplicative trend. The remaining IMFs are bit more complex, and one could argue that they are partly seasonal, trend or just some irregularity i.e. noise. 

Let us compare the decomposition with basic structural time series model fit from `StructTS` (for smoothing of more complex state space models, one could use [KFAS](https://github.com/helske/KFAS))

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

