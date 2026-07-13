safe_mc_cores <- function() {
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 2) {
    return(1L)
  }
  max(1L, cores - 1L)
}

#' Simple Random Testing Schedule
#'
#' Randomly tests eligible individuals with a fixed probability.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list containing `test.prob`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
sched_fun_srt <- function(t, state, test_control) {
  test_prob <- test_control$test.prob
  (stats::rbinom(nrow(state$R), 1, test_prob) == 1) & (state$R[, t] == 0)
}

#' Once-Per-Period Testing Schedule
#'
#' Tests individuals approximately once in each period among those who are
#' eligible and not currently removed.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list containing `period.length`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
sched_fun_opp <- function(t, state, test_control) {
  period_length <- test_control$period.length

  last_sched_test <- state$Z.prev.sched[, t]
  last_clear <- state$C.t[, t]
  last_period_boundary <- ((t - 1) %/% period_length) * period_length

  test_prob <- (1 / period_length) /
    (1 - (t - last_period_boundary - 1) / period_length)
  eligible <- (last_sched_test <= last_period_boundary |
    last_clear > last_sched_test) &
    state$R[, t] == FALSE

  stats::rbinom(nrow(state$D), 1, test_prob) == 1 & eligible
}

#' Maximum-Gap Testing Schedule
#'
#' Tests individuals with a probability that increases quadratically with the
#' elapsed time since their last scheduled test or clearance.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list containing `max.gap`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
sched_fun_maxgap <- function(t, state, test_control) {
  max_gap <- test_control$max.gap

  last_sched_test <- state$Z.prev.sched[, t]
  last_clearance <- state$C.t[, t]
  never_sched_tested <- (last_sched_test == 0) & (last_clearance == 0)

  test_prob <- ((t - pmax(last_sched_test, last_clearance)) / max_gap)^2
  test_prob[never_sched_tested] <-
    (1 / max_gap) / (1 - (t - 1) / max_gap)
  test_prob[state$R[, t] == 1] <- 0
  test_prob <- pmax(0, pmin(test_prob, 1))

  stats::rbinom(nrow(state$D), 1, test_prob) == 1
}

#' Min-Max Testing Schedule
#'
#' Uses a minimum testing gap and a quadratic maximum-gap schedule thereafter.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list containing `min.gap` and `max.gap`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
sched_fun_minmax <- function(t, state, test_control) {
  min_gap <- test_control$min.gap
  max_gap <- test_control$max.gap

  last_sched_test <- state$Z.prev.sched[, t]
  last_clearance <- state$C.t[, t]
  never_sched_tested <- (last_sched_test == 0) & (last_clearance == 0)

  test_prob <- ((t - pmax(last_sched_test, last_clearance)) / max_gap)^2
  test_prob[t - last_sched_test <= min_gap] <- 0
  test_prob[never_sched_tested] <-
    (1 / max_gap) / (1 - (t - 1) / max_gap)
  test_prob[state$R[, t] == 1] <- 0
  test_prob <- pmax(0, pmin(test_prob, 1))

  stats::rbinom(nrow(state$D), 1, test_prob) == 1
}

#' Symptom-Based Testing Component
#'
#' Generates additional tests induced by symptoms among infectious and
#' non-infectious individuals.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list containing `symptomatic.prob` and
#'   `symptom.false.prob`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
sympt_fun <- function(t, state, test_control) {
  test <- rep(FALSE, nrow(state$D))

  if (test_control$symptomatic.prob > 0) {
    if (t > 1) {
      newly_infectious <- state$compartment[, t] == "I" &
        state$compartment[, t - 1] == "W"
      symptomatic_true <- newly_infectious &
        stats::rbinom(
          nrow(state$D),
          1,
          prob = test_control$symptomatic.prob
        ) == 1
    } else {
      symptomatic_true <- FALSE
    }

    not_infectious <- state$compartment[, t] == "W"
    symptomatic_false <- not_infectious &
      stats::rbinom(
        nrow(state$D),
        1,
        prob = test_control$symptom.false.prob
      ) == 1

    test <- symptomatic_true | symptomatic_false
  }

  test
}

#' Contact-Tracing Testing Component
#'
#' Tests cluster members after a positive test in the preceding day.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control Unused placeholder for compatibility with `test_fun`.
#'
#' @return A logical vector indicating who is tested on day `t`.
#' @export
trace_fun <- function(t, state, test_control) {
  test <- rep(FALSE, nrow(state$D))
  clusters <- unique(state$cluster)

  if (length(clusters) < nrow(state$D) && t > 1) {
    for (i in clusters) {
      cluster <- state$cluster == i
      if (any(state$Y[cluster, t - 1] == TRUE)) {
        test[cluster & state$R[, t] == FALSE &
          state$tracing.mat[, t - 1] == FALSE] <- TRUE
      }
    }
  }

  test
}

#' Combined Testing Mechanism
#'
#' Combines scheduled, contact-tracing, and symptom-based testing.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param test_control A list of testing parameters.
#' @param sched_fun Scheduled-testing function. Defaults to `sched_fun_opp`.
#'
#' @return A logical vector indicating who is tested on day `t`, with testing
#'   components attached as attributes.
#' @export
test_fun <- function(t, state, test_control, sched_fun = sched_fun_opp) {
  test_sched <- sched_fun(t, state, test_control)
  test_trace <- trace_fun(t, state, test_control)
  test_sympt <- sympt_fun(t, state, test_control)

  test <- test_sched | test_trace | test_sympt

  attr(test, "sched") <- test_sched
  attr(test, "trace") <- test_trace
  attr(test, "sympt") <- test_sympt
  test
}

#' Quadratic Hazard Exposure Mechanism
#'
#' Generates new infections using a time-varying hazard with cluster effects.
#'
#' @param t Integer day index.
#' @param state A simulation state list produced inside the data generator.
#' @param expose_control A list containing `hazard.scale`,
#'   `resistance.scale`, `interval.length`, `max.hazard`, and `min.hazard`.
#'
#' @return A logical vector indicating who is newly exposed on day `t`.
#' @export
expose_fun_hazard_quadratic <- function(t, state, expose_control) {
  hazard_scale <- expose_control$hazard.scale
  resistance_scale <- expose_control$resistance.scale
  interval_length <- expose_control$interval.length
  max_hazard <- expose_control$max.hazard
  min_hazard <- expose_control$min.hazard

  resistant <- rowSums(state$X > -Inf) > 0
  n_people <- nrow(state$D)

  cluster_infectious_count <- rep(0, n_people)
  clusters <- unique(state$cluster)

  if (length(clusters) < n_people) {
    for (i in clusters) {
      cluster <- state$cluster == i
      cluster_infectious_count[cluster] <-
        sum(state$compartment[cluster, t] == "I")
    }
  }

  exposure_prob <- (1 -
    (1 - (t * (interval_length - t) / (interval_length / 2)^2 *
      (max_hazard - min_hazard) + min_hazard) * hazard_scale) *
      (1 - pmin(cluster_infectious_count / 10, 3 / 10))) *
    resistance_scale^resistant

  exposure <- as.logical(stats::rbinom(n_people, 1, exposure_prob))
  exposure & state$compartment[, t] == "W"
}

#' Generate Repeated-Testing Surveillance Data
#'
#' Simulates repeated testing, isolation, and infection dynamics under the
#' schedule used in the manuscript.
#'
#' @param n_clusters Number of clusters.
#' @param n_per_cluster Number of individuals per cluster.
#' @param max_days Number of study days.
#' @param initial_prev Initial prevalence.
#' @param test_fun Testing mechanism. Defaults to `test_fun`.
#' @param sched_fun Scheduled-testing component. Defaults to `sched_fun_opp`.
#' @param test_control A list of testing parameters.
#' @param expose_fun Exposure mechanism. Defaults to
#'   `expose_fun_hazard_quadratic`.
#' @param expose_control A list of exposure parameters.
#' @param recovery_allowed Logical indicating whether recovery is allowed.
#' @param recovery_days Recovery delay in days.
#' @param clearance_days Isolation clearance length in days.
#' @param sens Test sensitivity.
#' @param spec Test specificity.
#' @param time_varying_sens Logical; if `TRUE`, sensitivity varies with time
#'   since infection.
#' @param seed Random seed.
#'
#' @return A list of simulation matrices and metadata.
#' @export
#'
#' @examples
#' sim <- generate_surveillance_data(
#'   n_clusters = 5,
#'   n_per_cluster = 2,
#'   max_days = 7,
#'   initial_prev = 0.05,
#'   test_control = list(
#'     test.prob = 1 / 6,
#'     max.gap = 5,
#'     period.length = 7,
#'     min.gap = 3,
#'     symptomatic.prob = 0.25,
#'     symptom.false.prob = 0.01
#'   ),
#'   expose_control = list(
#'     hazard.scale = 0.1,
#'     resistance.scale = 0.5,
#'     interval.length = 7,
#'     max.hazard = 0.1,
#'     min.hazard = 0.02
#'   ),
#'   recovery_allowed = FALSE,
#'   recovery_days = 30,
#'   clearance_days = 5,
#'   sens = 0.832,
#'   spec = 0.992,
#'   seed = 1
#' )
#' names(sim)
generate_surveillance_data <- function(
    n_clusters,
    n_per_cluster,
    max_days,
    initial_prev,
    test_fun = NULL,
    sched_fun = sched_fun_opp,
    test_control = list(),
    expose_fun = expose_fun_hazard_quadratic,
    expose_control = list(),
    recovery_allowed,
    recovery_days,
    clearance_days,
    sens,
    spec,
    time_varying_sens = FALSE,
    seed) {
  if (is.null(test_fun)) {
    test_fun <- match.fun("test_fun")
  } else {
    test_fun <- match.fun(test_fun)
  }
  sched_fun <- match.fun(sched_fun)
  expose_fun <- match.fun(expose_fun)
  set.seed(seed)

  n_people <- n_clusters * n_per_cluster
  cluster <- rep(seq_len(n_clusters), each = n_per_cluster)
  cluster_virtual_test <- rep(0, n_people)

  D <- Y <- Z.prev <- Z.prev.sched <- R <- C.t <- X <- kt <- lt <-
    compartment <- matrix(NA, nrow = n_people, ncol = max_days)

  R[, 1] <- FALSE
  C.t[, 1] <- 0
  kt[, 1] <- 0
  lt[, 1] <- 1

  initial_infection <- stats::rbinom(n_people, 1, initial_prev) == 1
  compartment[, 1] <- ifelse(initial_infection, "I", "W")
  X[, 1] <- ifelse(
    initial_infection,
    sample(-(0:(recovery_days - 1)), size = n_people, replace = TRUE),
    -Inf
  )
  Z.prev[, 1] <- -Inf
  Z.prev.sched[, 1] <- 0

  schedule.mat <- symptomatic.mat <- tracing.mat <-
    matrix(FALSE, nrow = n_people, ncol = max_days)

  for (t in seq_len(max_days)) {
    current_state <- list(
      D = D[, 1:t, drop = FALSE],
      Y = Y[, 1:t, drop = FALSE],
      Z.prev = Z.prev[, 1:t, drop = FALSE],
      Z.prev.sched = Z.prev.sched[, 1:t, drop = FALSE],
      R = R[, 1:t, drop = FALSE],
      C.t = C.t[, 1:t, drop = FALSE],
      X = X[, 1:t, drop = FALSE],
      kt = kt[, 1:t, drop = FALSE],
      lt = lt[, 1:t, drop = FALSE],
      compartment = compartment[, 1:t, drop = FALSE],
      cluster = cluster,
      cluster.virtual.test = cluster_virtual_test,
      schedule.mat = schedule.mat,
      symptomatic.mat = symptomatic.mat,
      tracing.mat = tracing.mat
    )

    tested <- test_fun(
      t = t,
      state = current_state,
      test_control = test_control,
      sched_fun = sched_fun
    )

    D[, t] <- tested
    schedule.mat[, t] <- attr(tested, "sched")
    tracing.mat[, t] <- attr(tested, "trace")
    symptomatic.mat[, t] <- attr(tested, "sympt")

    current_state$D <- D[, 1:t, drop = FALSE]
    current_state$schedule.mat <- schedule.mat[, 1:t, drop = FALSE]
    current_state$symptomatic.mat <- symptomatic.mat[, 1:t, drop = FALSE]
    current_state$tracing.mat <- tracing.mat[, 1:t, drop = FALSE]

    if (time_varying_sens) {
      sens_all <- sens * pmax((t - X[, t]) * (10 - t + X[, t]) / (10 / 2)^2, 1 / 10)
      sens_all <- pmax(0, pmin(sens_all, sens))
    } else {
      sens_all <- sens
    }

    Y[, t] <- D[, t] & (
      (compartment[, t] == "I" &
        stats::rbinom(n_people, 1, prob = sens_all) == 1) |
        (compartment[, t] == "W" &
          stats::rbinom(n_people, 1, prob = 1 - spec) == 1)
    )

    current_state$Y <- Y[, 1:t, drop = FALSE]

    if (t < max_days) {
      Z.prev[, t + 1] <- ifelse(tested, t, Z.prev[, t])
      Z.prev.sched[, t + 1] <- ifelse(schedule.mat[, t], t, Z.prev.sched[, t])
      kt[, t + 1] <- ifelse(tested, kt[, t] + 1, kt[, t])
      lt[, t + 1] <- lt[, t]
      X[, t + 1] <- X[, t]
      compartment[, t + 1] <- compartment[, t]
      R[, t + 1] <- R[, t]
      C.t[, t + 1] <- C.t[, t]
    }

    if (t < max_days) {
      tested_positive <- Y[, t]
      R[tested_positive, t + 1] <- TRUE
      compartment[tested_positive, t + 1] <- "R"

      if (t > clearance_days) {
        cleared_today <- rowMeans(R[, (t - clearance_days + 1):t, drop = FALSE]) == 1
        R[cleared_today, t + 1] <- FALSE
        C.t[cleared_today, t + 1] <- t
        lt[cleared_today, t + 1] <- lt[cleared_today, t] + 1
        compartment[cleared_today, t + 1] <- "W"
      }
    }

    exposed <- expose_fun(t = t, state = current_state, expose_control = expose_control)
    exposed <- exposed & !R[, t]

    if (t < max_days) {
      compartment[exposed, t + 1] <- "I"
      X[exposed, t + 1] <- t
    }

    if (recovery_allowed && t < max_days) {
      recovered_today <- X[, t] < t - recovery_days &
        rowMeans(
          compartment[, max(t - recovery_days + 1, 1):t, drop = FALSE] == "I"
        ) == 1
      compartment[recovered_today, t + 1] <- "W"
    }
  }

  list(
    D = D,
    Y = Y,
    Z.prev = Z.prev,
    Z.prev.sched = Z.prev.sched,
    R = R,
    C.t = C.t,
    X = X,
    kt = kt,
    lt = lt,
    compartment = compartment,
    cluster = cluster,
    schedule.mat = schedule.mat,
    symptomatic.mat = symptomatic.mat,
    tracing.mat = tracing.mat
  )
}
