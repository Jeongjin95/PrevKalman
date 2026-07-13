# PrevKalman

`PrevKalman` is an R package for Horvitz-Thompson prevalence estimation under repeated testing, block jackknife uncertainty quantification, and Kalman filtering/smoothing of daily prevalence trajectories.

## Installation

```r
remotes::install_github("Jeongjin95/PrevKalman")
```

## Main functions

- `ht_prevalence()` computes the daily Horvitz-Thompson prevalence estimate and the default 20-block jackknife confidence interval.
- `fit_kalman_model()` fits the joint Kalman model and returns estimates for either the filter or smoother with three CI options: `model_based`, `ht_se`, and `ht_se_adjusted`.

## Example

```r
library(PrevKalman)

#' @param D i*t logical matrix indicating whether unit i is tested at time t
#' @param Y i*t logical matrix indicating whether unit i tests positive
#'          at time t; FALSE includes negative or untested observations
#' @param R i*t logical matrix indicating whether unit i is removed from the
#'          risk set at time t
#' @param Z.next i*t numeric matrix indicating next observed test time of unit i
#'               at time t (next test time may be t; Inf if no later test)
#' @param C.t i*t numeric matrix indicating last clearance time of unit i
#'            at time t
#' @param Sympt.t i*t numeric matrix indicating last symptom/contact-tracing
#'                trigger time of unit i at time t
#' @param gamma Numeric vector giving support values of C.t, typically computed
#'        as sort(unique(as.vector(C.t)))
#' @param SPEC Test specificity (0--1 scale)
#' @param SENS Test sensitivity (0--1 scale)

ht_fit <- ht_prevalence(
  D = D,
  Y = Y,
  R = R,
  Z.next = Z.next,
  C.t = C.t,
  Sympt.t = Sympt.t,
  gamma = gamma,
  SPEC = 1,
  SENS = 0.832,
  max.days = ncol(D),
  n_blocks = 20,
  n_cores = 1
)

kf_fit <- fit_kalman_model(
  x = ht_fit,
  state = "filter",
  ci_method = "all"
)

ks_fit <- fit_kalman_model(
  x = ht_fit,
  state = "smoother",
  ci_method = "all"
)
```

## GitHub

The package metadata currently assumes the GitHub repository will be:

- `https://github.com/Jeongjin95/PrevKalman`
