# PrevKalman

`PrevKalman` is an R package for Horvitz-Thompson prevalence estimation under repeated testing, block jackknife uncertainty quantification, and Kalman filtering/smoothing of daily prevalence trajectories.

## Installation

```r
remotes::install_github("Jeongjin95/PrevKalman")
```

## Main functions

- `generate_surveillance_data()` generates simulated repeated-testing surveillance data with scheduled, symptomatic, and contact-tracing testing components.
- `ht_prevalence()` computes the daily Horvitz-Thompson prevalence estimate and the default 20-block jackknife confidence interval.
- `fit_kalman_model()` fits the joint Kalman model and returns estimates for either the filter or smoother with three CI options: `model_based`, `ht_se`, and `ht_se_adjusted`.

## Example

```r
library(PrevKalman)

# This example shows the full workflow:
# 1. simulate repeated-testing surveillance data,
# 2. construct the inputs required by ht_prevalence(),
# 3. estimate daily prevalence using the HT estimator with block jackknife CIs,
# 4. fit the joint Kalman filter/smoother with three CI options.

# 1. Generate toy surveillance data.
# Scheduled testing uses once-per-period testing via sched_fun_opp.
# The symptomatic/contact-tracing testing component is included through test_fun().
sim <- generate_surveillance_data(
  n_clusters = 5000,
  n_per_cluster = 2,
  max_days = 21,
  initial_prev = 0.02,
  test_control = list(
    test.prob = 1 / 6,
    max.gap = 10,
    period.length = 7,
    min.gap = 5,
    symptomatic.prob = 1 / 4,
    symptom.false.prob = 1 / 100
  ),
  expose_control = list(
    hazard.scale = 1 / 10,
    resistance.scale = 1 / 2,
    interval.length = 21,
    max.hazard = 0.1,
    min.hazard = 0.02
  ),
  test_fun = test_fun,
  sched_fun = sched_fun_opp,
  expose_fun = expose_fun_hazard_quadratic,
  recovery_allowed = FALSE,
  recovery_days = 30,
  clearance_days = 5,
  sens = 0.832,
  spec = 0.992,
  seed = 1
)

# 2. Prepare inputs for ht_prevalence().
D <- sim$D          # testing indicator matrix
Y <- sim$Y          # test-result matrix
R <- sim$R          # removed/exempt-state matrix
C.t <- sim$C.t      # most recent clearance-time matrix
max.days <- ncol(D)

# Z.next[i, t] is the next observed test day for unit i at or after day t.
Z.next <- matrix(NA_real_, nrow = nrow(D), ncol = ncol(D))
Z.next[, max.days] <- ifelse(D[, max.days], max.days, Inf)

if (max.days > 1) {
  for (t in seq.int(max.days - 1, 1)) {
    Z.next[, t] <- ifelse(D[, t], t, Z.next[, t + 1])
  }
}

# Sympt.t[i, t] is the most recent symptom/contact-tracing trigger day
# for unit i up to day t.
Sympt.t <- matrix(0L, nrow = nrow(D), ncol = ncol(D))

for (i in seq_len(nrow(D))) {
  last_time <- 0L

  for (t in seq_len(max.days)) {
    if (sim$symptomatic.mat[i, t] || sim$tracing.mat[i, t]) {
      last_time <- t
    }

    Sympt.t[i, t] <- last_time
  }
}

# Support values for the clearance-time history.
gamma <- sort(unique(as.vector(C.t)))

# 3. Horvitz-Thompson prevalence estimate + block jackknife CI.
# Use n_blocks = 5 for a faster example. Use n_blocks = 20 for the paper.
ht_fit <- ht_prevalence(
  D = D,
  Y = Y,
  R = R,
  Z.next = Z.next,
  C.t = C.t,
  Sympt.t = Sympt.t,
  gamma = gamma,
  SPEC = 0.992,
  SENS = 0.832,
  max.days = max.days,
  n_blocks = 5,
  n_cores = 1
)

head(ht_fit$results)

# 4. Joint Kalman filter/smoother with three CI options:
# model_based, ht_se, and ht_se_adjusted.
kalman_fit <- fit_kalman_model(
  x = ht_fit,
  state = "both",
  ci_method = "all"
)

head(kalman_fit$results)
```

## GitHub

The package metadata currently assumes the GitHub repository will be:

- `https://github.com/Jeongjin95/PrevKalman`
