# tvp-qr

## Description
> These codes come without technical support of any kind. The code is free to use, provided that the paper is cited properly.

Codes based on M. Pfarrhofer (2022): "[Modeling tail risks of inflation using unobserved component quantile regressions](https://doi.org/10.1016/j.jedc.2022.104493)" _Journal of Economic Dynamics and Control_ *143* 104493, using time-varying parameter quantile regressions (TVP-QR) with dynamic shrinkage priors and a time-varying scale parameter.

The file `data_raw.rda` contains the quarterly inflation series for the United States (`US`, 1947–2021), the United Kingdom (`UK`, 1960–2021) and the euro area (`EA`, 1990–2021).

## Source files
- `!example` is the main file which loads the data, estimates unobserved component quantile regressions on a grid of quantiles and plots the estimated quantiles over time (requires user input)
- `!qrdhs` contains the MCMC sampler `tvpqr()` for a single quantile and the wrapper `tvpqr.grid()` which estimates a grid of quantiles in parallel (is sourced automatically from `!example`)
- `aux`, `ffbs.cpp` and `jpr_qr.cpp` contain helper functions, including forward filtering backward sampling for the time-varying parameters and the JPR (Jacquier, Polson & Rossi) sampler for the time-varying scale parameter (are sourced automatically from `!qrdhs`)
- `dsp_aux` implements the dynamic horseshoe prior, adapted from the [dsp](https://github.com/drkowal/dsp) R-package by D. R. Kowal (see Kowal, Matteson & Ruppert, 2019, "[Dynamic shrinkage processes](https://doi.org/10.1111/rssb.12325)", _JRSS-B_ *81*(4) 781–804)

All files are sourced by relative path, so the working directory must be the repository root.

### Estimation options in `!example`
- `sl.cn` selects the country, choose from `US`, `UK` or `EA`
- `grid.p` sets the grid of quantiles
- `prior` refers to the prior on the time-varying parameters, choose from:
  - `dhs` dynamic horseshoe
  - `shs` horseshoe prior on the state innovations (time-specific local scales)
  - `iG` inverse gamma prior on constant state innovation variances
- MCMC settings (`nburn`, `nsave`, `thinfac`) are set in the file, and the number of cores via `cpu` in `tvpqr.grid()`

Both specifications are estimated, with a time-invariant (`UCQR-TIS`, `sv=FALSE`) and a time-varying scale parameter (`UCQR-TVS`, `sv=TRUE`). The posterior means of the quantiles are plotted over time for both models.

### Requirements
R packages `Rcpp`, `RcppArmadillo` (with a working C++ compiler), `Matrix`, `MASS`, `spam`, `pgdraw`, `stochvol`, `GIGrvg`, `invgamma`, `LaplacesDemon`, `doParallel` and `foreach` for estimation; `lubridate`, `reshape2`, `dplyr`, `tidyr`, `ggplot2`, `cowplot` and `lemon` for the output.

## License
GPL-2, since the repository includes code adapted from the GPL-2 licensed [dsp](https://github.com/drkowal/dsp) package.
