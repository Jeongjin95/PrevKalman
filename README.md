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

# D: testing indicator matrix
# Y: test-result matrix
# R: removed/exempt-state matrix
# Z.next: next testing-time matrix
# C.t: most recent clearance-time matrix
# Sympt.t: most recent symptom/contact-time matrix
# gamma: support values of C.t

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

If you choose a different repository name, update the `URL` and `BugReports` fields in `DESCRIPTION`.
