test_that("kalman model fitting returns three joint CI options", {
  y <- c(0.01, 0.015, NA, 0.03, 0.025, 0.02)
  r_t <- c(0.0001, 0.0001, NA, 0.0002, 0.00015, 0.00012)

  fit <- fit_kalman_model(
    y = y,
    r_t = r_t,
    state = "filter",
    ci_method = "all"
  )

  expect_equal(unique(fit$results$model), "Joint")
  expect_equal(
    sort(unique(fit$results$ci_method)),
    sort(c("model_based", "ht_se", "ht_se_adjusted"))
  )
  expect_equal(unique(fit$results$state), "filter")
  expect_equal(nrow(fit$results), 3 * length(y))
  expect_true(all(fit$results$estimate >= 0 & fit$results$estimate <= 1))
  expect_true(all(fit$results$lower <= fit$results$upper, na.rm = TRUE))
})

test_that("fit_mixture supports convex-combination CI-6", {
  sim <- generate_surveillance_data(
    n_clusters = 4,
    n_per_cluster = 2,
    max_days = 5,
    initial_prev = 0.05,
    test_control = list(
      test.prob = 1 / 6,
      max.gap = 5,
      period.length = 7,
      min.gap = 3,
      symptomatic.prob = 0.25,
      symptom.false.prob = 0.01
    ),
    expose_control = list(
      hazard.scale = 0.1,
      resistance.scale = 0.5,
      interval.length = 5,
      max.hazard = 0.1,
      min.hazard = 0.02
    ),
    recovery_allowed = FALSE,
    recovery_days = 30,
    clearance_days = 3,
    sens = 0.832,
    spec = 0.992,
    seed = 1
  )

  Z_next <- matrix(NA_real_, nrow = nrow(sim$D), ncol = ncol(sim$D))
  Z_next[, ncol(sim$D)] <- ifelse(sim$D[, ncol(sim$D)] == 1, ncol(sim$D), Inf)
  if (ncol(sim$D) > 1) {
    for (t in seq.int(ncol(sim$D) - 1, 1)) {
      Z_next[, t] <- ifelse(sim$D[, t] == 1, t, Z_next[, t + 1])
    }
  }

  Sympt_t <- matrix(0L, nrow = nrow(sim$D), ncol = ncol(sim$D))
  for (i in seq_len(nrow(sim$D))) {
    last_time <- 0L
    for (t in seq_len(ncol(sim$D))) {
      if (sim$symptomatic.mat[i, t] | sim$tracing.mat[i, t]) {
        last_time <- t
      }
      Sympt_t[i, t] <- last_time
    }
  }

  ht_fit <- ht_prevalence(
    D = sim$D,
    Y = sim$Y,
    R = sim$R,
    Z.next = Z_next,
    C.t = sim$C.t,
    Sympt.t = Sympt_t,
    gamma = unique(as.vector(sim$C.t)),
    SPEC = 0.992,
    SENS = 0.832,
    max.days = ncol(sim$D),
    n_blocks = 4,
    n_cores = 1
  )

  fit <- fit_mixture(
    x = ht_fit,
    state = "smoother"
  )

  expect_equal(unique(fit$results$ci_method), "convex_combination")
  expect_equal(unique(fit$results$state), "smoother")
  expect_true(all(is.finite(fit$results$weight), na.rm = TRUE))
  expect_true(all(fit$results$weight >= 0 & fit$results$weight <= 1, na.rm = TRUE))
})

test_that("kalman model fitting supports smoother output", {
  y <- c(0.01, 0.015, NA, 0.03, 0.025, 0.02)
  r_t <- c(0.0001, 0.0001, NA, 0.0002, 0.00015, 0.00012)

  fit <- fit_kalman_model(
    y = y,
    r_t = r_t,
    state = "smoother",
    ci_method = "model_based"
  )

  expect_equal(unique(fit$results$model), "Joint")
  expect_equal(unique(fit$results$state), "smoother")
  expect_equal(nrow(fit$results), length(y))
})

test_that("ht_prevalence returns estimate and jackknife CI", {
  sim <- generate_surveillance_data(
    n_clusters = 4,
    n_per_cluster = 2,
    max_days = 5,
    initial_prev = 0.05,
    test_control = list(
      test.prob = 1 / 6,
      max.gap = 5,
      period.length = 7,
      min.gap = 3,
      symptomatic.prob = 0.25,
      symptom.false.prob = 0.01
    ),
    expose_control = list(
      hazard.scale = 0.1,
      resistance.scale = 0.5,
      interval.length = 5,
      max.hazard = 0.1,
      min.hazard = 0.02
    ),
    recovery_allowed = FALSE,
    recovery_days = 30,
    clearance_days = 3,
    sens = 0.832,
    spec = 0.992,
    seed = 1
  )

  Z_next <- matrix(NA_real_, nrow = nrow(sim$D), ncol = ncol(sim$D))
  Z_next[, ncol(sim$D)] <- ifelse(sim$D[, ncol(sim$D)] == 1, ncol(sim$D), Inf)
  if (ncol(sim$D) > 1) {
    for (t in seq.int(ncol(sim$D) - 1, 1)) {
      Z_next[, t] <- ifelse(sim$D[, t] == 1, t, Z_next[, t + 1])
    }
  }

  Sympt_t <- matrix(0L, nrow = nrow(sim$D), ncol = ncol(sim$D))
  for (i in seq_len(nrow(sim$D))) {
    last_time <- 0L
    for (t in seq_len(ncol(sim$D))) {
      if (sim$symptomatic.mat[i, t] | sim$tracing.mat[i, t]) {
        last_time <- t
      }
      Sympt_t[i, t] <- last_time
    }
  }

  est <- ht_prevalence(
    D = sim$D,
    Y = sim$Y,
    R = sim$R,
    Z.next = Z_next,
    C.t = sim$C.t,
    Sympt.t = Sympt_t,
    gamma = unique(as.vector(sim$C.t)),
    SPEC = 0.992,
    SENS = 0.832,
    max.days = ncol(sim$D),
    n_blocks = 4,
    n_cores = 1
  )

  expect_equal(nrow(est$results), ncol(sim$D))
  expect_true(all(is.finite(est$estimate)))
  expect_true(all(is.finite(est$lower)))
  expect_true(all(is.finite(est$upper)))
  expect_true(all(est$lower <= est$upper))
  expect_true(all(c("D", "Y", "R", "Z_next", "C_t", "Sympt_t", "gamma", "days") %in% names(est$inputs)))
})

test_that("ht_prevalence returns NA on days with no testing", {
  prep <- list(
    D = matrix(c(1L, 0L, 1L,
                 1L, 0L, 1L), nrow = 2, byrow = TRUE),
    Y = matrix(0L, nrow = 2, ncol = 3),
    R = matrix(FALSE, nrow = 2, ncol = 3),
    Z_next = matrix(c(1, 3, 3,
                      1, 3, 3), nrow = 2, byrow = TRUE),
    C_t = matrix(0L, nrow = 2, ncol = 3),
    Sympt_t = matrix(0L, nrow = 2, ncol = 3),
    gamma = 0L,
    days = 1:3
  )

  est <- ht_prevalence(
    prepared_data = prep,
    spec = 0.992,
    sens = 0.832,
    n_blocks = 2,
    n_cores = 1
  )

  expect_equal(est$no_test_days, 2)
  expect_true(is.na(est$results$estimate[2]))
  expect_true(is.na(est$results$se[2]))
  expect_true(is.na(est$results$lower[2]))
  expect_true(is.na(est$results$upper[2]))
  expect_true(all(is.na(est$jackknife$theta_leave_block_out[, 2])))
})
