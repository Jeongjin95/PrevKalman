# PrevKalman

`PrevKalman` is an R package for Horvitz-Thompson prevalence estimation under repeated testing, block jackknife uncertainty quantification, and Kalman filtering/smoothing of daily prevalence trajectories.

## Installation

```r
remotes::install_github("Jeongjin95/PrevKalman")
```

## Main functions

- `prepare_surveillance_data()` converts long-form repeated testing data into the matrices needed by the estimator.
- `ht_prevalence()` computes the daily Horvitz-Thompson prevalence estimate and the default 20-block jackknife confidence interval.
- `fit_kalman_model()` fits the joint Kalman model and returns estimates for either the filter or smoother with three CI options: `model_based`, `ht_se`, and `ht_se_adjusted`.

## Example

```r
library(PrevKalman)

prep <- prepare_surveillance_data(
  osu_surveillance,
  low_test_threshold = 100,
  symptom_unobserved_days = 1:12
)

ht_fit <- ht_prevalence(
  prepared_data = prep,
  spec = 1,
  sens = 0.832,
  n_blocks = 20,
  n_cores = 1
)

kf_fit <- fit_kalman_model(
  ht_fit,
  state = "filter",
  ci_method = "all"
)

head(ht_fit$results)
head(kf_fit$results)
```

## GitHub

The package metadata currently assumes the GitHub repository will be:

- `https://github.com/Jeongjin95/PrevKalman`

If you choose a different repository name, update the `URL` and `BugReports` fields in `DESCRIPTION`.
