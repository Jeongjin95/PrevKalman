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

# D: binary testing-indicator matrix.
#    D[i, t] = 1 if individual i was tested on day t, and 0 otherwise.
#
# Y: binary test-result matrix.
#    Y[i, t] = 1 if individual i tested positive on day t, 0 if negative,
#    and NA when no test was administered.
#
# R: removed/exempt-state matrix.
#    R[i, t] = 1 if individual i was no longer in the risk set on day t
#    because of removal, exemption, or leaving surveillance; 0 otherwise.
#
# Z.next: next observed testing-time matrix.
#    Z.next[i, t] is the next day at or after day t on which individual i
#    is observed to be tested, computed from D. It is Inf if individual i
#    has no observed test on or after day t.
#
# C.t: most recent clearance-time matrix.
#    C.t[i, t] records the most recent day up to day t on which individual i
#    was cleared or became eligible again after a previous positive test,
#    isolation period, or other temporary ineligibility.
#
# Sympt.t: most recent symptom/contact-time matrix.
#    Sympt.t[i, t] records the most recent day up to day t on which individual i
#    had a symptom-driven or contact-tracing-driven testing trigger.
#
# gamma: support values of C.t.
#    gamma is usually the sorted unique set of values observed in C.t,
#    for example gamma = sort(unique(as.vector(C.t))).

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
