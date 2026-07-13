subset_data_list_by_keep <- function(data_list, keep) {
  lapply(data_list, function(x) {
    if (is.matrix(x)) {
      x[keep, , drop = FALSE]
    } else {
      x
    }
  })
}

normalize_prepared_data_names <- function(prepared_data) {
  if (is.null(prepared_data)) {
    return(NULL)
  }

  if (is.null(prepared_data$Z_next) && !is.null(prepared_data$Z.next)) {
    prepared_data$Z_next <- prepared_data$Z.next
  }
  if (is.null(prepared_data$C_t) && !is.null(prepared_data$C.t)) {
    prepared_data$C_t <- prepared_data$C.t
  }
  if (is.null(prepared_data$Sympt_t) && !is.null(prepared_data$Sympt.t)) {
    prepared_data$Sympt_t <- prepared_data$Sympt.t
  }
  if (is.null(prepared_data$days) && !is.null(prepared_data$max.days)) {
    prepared_data$days <- seq_len(prepared_data$max.days)
  }

  prepared_data
}

resolve_ht_inputs <- function(
    prepared_data = NULL,
    D = NULL,
    Y = NULL,
    R = NULL,
    Z_next = NULL,
    C_t = NULL,
    Sympt_t = NULL,
    gamma = NULL,
    days = NULL) {
  prepared_data <- normalize_prepared_data_names(prepared_data)

  if (!is.null(prepared_data)) {
    D <- prepared_data$D
    Y <- prepared_data$Y
    R <- prepared_data$R
    Z_next <- prepared_data$Z_next
    C_t <- prepared_data$C_t
    Sympt_t <- prepared_data$Sympt_t
    gamma <- prepared_data$gamma
    if (is.null(days) && !is.null(prepared_data$days)) {
      days <- prepared_data$days
    }
  }

  inputs <- list(
    D = D,
    Y = Y,
    R = R,
    Z_next = Z_next,
    C_t = C_t,
    Sympt_t = Sympt_t,
    gamma = gamma
  )

  missing_names <- names(inputs)[vapply(inputs, is.null, logical(1))]
  if (length(missing_names) > 0) {
    stop(
      "Missing required HT inputs: ",
      paste(missing_names, collapse = ", ")
    )
  }

  if (is.null(days)) {
    days <- seq_len(ncol(D))
  }

  c(inputs, list(days = days))
}

extract_kalman_inputs <- function(x = NULL, y = NULL, r_t = NULL, days = NULL) {
  ht_se <- NULL
  if (!is.null(x)) {
    if (is.list(x) && !is.null(x$results)) {
      if (is.null(y)) {
        y <- x$results$estimate
      }
      if (is.null(r_t)) {
        r_t <- x$results$se^2
      }
      if ("se" %in% names(x$results)) {
        ht_se <- x$results$se
      }
      if (is.null(days) && "day" %in% names(x$results)) {
        days <- x$results$day
      }
    } else if (is.data.frame(x)) {
      if (is.null(y) && "estimate" %in% names(x)) {
        y <- x$estimate
      }
      if (is.null(r_t) && "se" %in% names(x)) {
        r_t <- x$se^2
      }
      if ("se" %in% names(x)) {
        ht_se <- x$se
      }
      if (is.null(days) && "day" %in% names(x)) {
        days <- x$day
      }
    }
  }

  if (is.null(y) || is.null(r_t)) {
    stop("Provide either `x` with estimate/SE columns or both `y` and `r_t`.")
  }

  if (is.null(days)) {
    days <- seq_along(y)
  }

  if (is.null(ht_se)) {
    ht_se <- sqrt(r_t)
  }

  list(y = y, r_t = r_t, ht_se = ht_se, days = days)
}

ci5_coverage_probability <- function(angle, z_value) {
  left_term <- z_value / cos(angle) - tan(angle)
  right_term <- -z_value / cos(angle) - tan(angle)
  stats::pnorm(left_term) - stats::pnorm(right_term)
}

ci5_critical_value_from_angle <- function(angle, conf_level = 0.95) {
  alpha <- 1 - conf_level
  z_standard <- stats::qnorm(1 - alpha / 2)

  stats::uniroot(
    function(z_value) ci5_coverage_probability(angle, z_value) - conf_level,
    interval = c(0, z_standard)
  )$root
}

compute_ci_adjusted_critical <- function(ht_se, model_se, conf_level = 0.95) {
  z_standard <- stats::qnorm(1 - (1 - conf_level) / 2)

  vapply(seq_along(ht_se), function(i) {
    if (!is.finite(ht_se[i]) || ht_se[i] <= 0 ||
        !is.finite(model_se[i]) || model_se[i] <= 0 ||
        model_se[i] >= ht_se[i]) {
      return(z_standard)
    }

    angle <- acos(pmin(1, pmax(model_se[i] / ht_se[i], 0)))
    ci5_critical_value_from_angle(angle, conf_level = conf_level)
  }, numeric(1))
}

compute_ci6_angle <- function(se_ht, se_biased) {
  acos(pmin(1, pmax(se_biased / se_ht, 0)))
}

compute_ci6_se <- function(se_ht, se_biased, weight, rho) {
  sqrt(
    (1 - weight)^2 * se_ht^2 +
      weight^2 * se_biased^2 +
      2 * weight * (1 - weight) * rho * se_ht * se_biased
  )
}

ci6_coverage_probability <- function(angle, z_value, weight, rho) {
  denominator <- sqrt(
    (1 - weight)^2 +
      weight^2 * cos(angle)^2 +
      2 * rho * weight * (1 - weight) * cos(angle)
  )

  left_term <- (z_value - weight * sin(angle)) / denominator
  right_term <- (-z_value - weight * sin(angle)) / denominator
  stats::pnorm(left_term) - stats::pnorm(right_term)
}

ci6_critical_value <- function(angle, weight, rho, conf_level = 0.95) {
  stats::uniroot(
    function(z_value) {
      ci6_coverage_probability(angle, z_value, weight, rho) - conf_level
    },
    interval = c(0, 10)
  )$root
}

choose_ci6_weight <- function(se_ht, se_biased, rho, conf_level = 0.95) {
  z_standard <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (!is.finite(se_ht) || se_ht <= 0 ||
      !is.finite(se_biased) || se_biased <= 0 ||
      !is.finite(rho)) {
    return(0)
  }

  angle <- compute_ci6_angle(se_ht, se_biased)
  weight_grid <- seq(0, 1, by = 0.01)

  z_values <- vapply(weight_grid, function(weight) {
    tryCatch(
      ci6_critical_value(angle, weight, rho, conf_level = conf_level),
      error = function(...) Inf
    )
  }, numeric(1))

  if (all(!is.finite(z_values))) {
    return(0)
  }

  weight_grid[which.min(z_values)]
}

compute_ci6_weighted_results <- function(
    ht_estimate,
    model_estimate,
    ht_se,
    model_se,
    rho,
    conf_level = 0.95) {
  z_standard <- stats::qnorm(1 - (1 - conf_level) / 2)

  out <- lapply(seq_along(ht_estimate), function(i) {
    if (!is.finite(ht_estimate[i]) || !is.finite(model_estimate[i]) ||
        !is.finite(ht_se[i]) || ht_se[i] <= 0 ||
        !is.finite(model_se[i]) || model_se[i] <= 0) {
      return(list(
        estimate = NA_real_,
        se = NA_real_,
        critical_value = z_standard,
        lower = NA_real_,
        upper = NA_real_,
        weight = NA_real_,
        rho = NA_real_,
        combo_se = NA_real_
      ))
    }

    rho_i <- rho[i]
    if (!is.finite(rho_i)) {
      rho_i <- 0
    }
    rho_i <- pmin(1, pmax(-1, rho_i))

    weight_i <- choose_ci6_weight(
      se_ht = ht_se[i],
      se_biased = model_se[i],
      rho = rho_i,
      conf_level = conf_level
    )

    angle_i <- compute_ci6_angle(ht_se[i], model_se[i])
    critical_i <- tryCatch(
      ci6_critical_value(angle_i, weight_i, rho_i, conf_level = conf_level),
      error = function(...) z_standard
    )

    estimate_i <- (1 - weight_i) * ht_estimate[i] + weight_i * model_estimate[i]
    combo_se_i <- compute_ci6_se(ht_se[i], model_se[i], weight_i, rho_i)
    lower_i <- estimate_i - critical_i * ht_se[i]
    upper_i <- estimate_i + critical_i * ht_se[i]

    list(
      estimate = estimate_i,
      se = ht_se[i],
      critical_value = critical_i,
      lower = lower_i,
      upper = upper_i,
      weight = weight_i,
      rho = rho_i,
      combo_se = combo_se_i
    )
  })

  do.call(rbind, lapply(out, function(x) {
    data.frame(
      estimate = x$estimate,
      se = x$se,
      critical_value = x$critical_value,
      lower = x$lower,
      upper = x$upper,
      weight = x$weight,
      rho = x$rho,
      combo_se = x$combo_se
    )
  }))
}

build_joint_ci_results <- function(fit, state, ci_method, days, ht_se, conf_level) {
  if (state == "filter") {
    estimate <- fit$filtered
    model_se <- fit$filtered_se
  } else {
    estimate <- fit$smoothed
    model_se <- fit$smoothed_se
  }

  z_standard <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (ci_method == "model_based") {
    se <- model_se
    critical_value <- rep(z_standard, length(se))
  } else if (ci_method == "ht_se") {
    se <- ht_se
    critical_value <- rep(z_standard, length(se))
  } else if (ci_method == "ht_se_adjusted") {
    se <- ht_se
    critical_value <- compute_ci_adjusted_critical(
      ht_se = ht_se,
      model_se = model_se,
      conf_level = conf_level
    )
  } else {
    stop("Unknown ci_method: ", ci_method)
  }

  lower <- estimate - critical_value * se
  upper <- estimate + critical_value * se
  lower[!is.finite(lower)] <- NA_real_
  upper[!is.finite(upper)] <- NA_real_

  data.frame(
    day = days,
    model = "Joint",
    state = state,
    ci_method = ci_method,
    estimate = estimate,
    se = se,
    critical_value = critical_value,
    lower = lower,
    upper = upper
  )
}

#' Horvitz-Thompson Daily Prevalence Estimate
#'
#' Computes the daily prevalence estimator under repeated testing using the
#' manuscript's weighting construction.
#'
#' @param D Binary testing indicator matrix.
#' @param Y Binary test-result matrix.
#' @param R Logical removed-state matrix.
#' @param Z_next Matrix of next testing times.
#' @param C_t Matrix of most recent clearance times.
#' @param Sympt_t Matrix of most recent symptom times.
#' @param gamma Clearance support values, typically `unique(as.vector(C_t))`.
#' @param spec Test specificity.
#' @param sens Test sensitivity.
#' @param max_days Number of study days. Defaults to `ncol(D)`.
#'
#' @return A numeric vector of daily prevalence estimates.
#' @keywords internal
#' @noRd
#'
#' @examples
#' sim <- generate_surveillance_data(
#'   n_clusters = 4,
#'   n_per_cluster = 2,
#'   max_days = 6,
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
#'     interval.length = 6,
#'     max.hazard = 0.1,
#'     min.hazard = 0.02
#'   ),
#'   recovery_allowed = FALSE,
#'   recovery_days = 30,
#'   clearance_days = 3,
#'   sens = 0.832,
#'   spec = 0.992,
#'   seed = 1
#' )
#' Z_next <- matrix(NA_real_, nrow = nrow(sim$D), ncol = ncol(sim$D))
#' Z_next[, ncol(sim$D)] <- ifelse(sim$D[, ncol(sim$D)] == 1, ncol(sim$D), Inf)
#' for (t in seq.int(ncol(sim$D) - 1, 1)) {
#'   Z_next[, t] <- ifelse(sim$D[, t] == 1, t, Z_next[, t + 1])
#' }
#' Sympt_t <- matrix(0, nrow = nrow(sim$D), ncol = ncol(sim$D))
#' for (i in seq_len(nrow(sim$D))) {
#'   last_time <- 0
#'   for (t in seq_len(ncol(sim$D))) {
#'     if (sim$symptomatic.mat[i, t] | sim$tracing.mat[i, t]) {
#'       last_time <- t
#'     }
#'     Sympt_t[i, t] <- last_time
#'   }
#' }
#' ht_prevalence(
#'   D = sim$D,
#'   Y = sim$Y,
#'   R = sim$R,
#'   Z_next = Z_next,
#'   C_t = sim$C.t,
#'   Sympt_t = Sympt_t,
#'   gamma = unique(as.vector(sim$C.t)),
#'   spec = 0.992,
#'   sens = 0.832
#' )
ht_prevalence_core <- function(
    D,
    Y,
    R,
    Z_next,
    C_t,
    Sympt_t,
    gamma,
    spec,
    sens,
    max_days = ncol(D)) {
  gamma_vals <- sort(unique(gamma))
  n_gamma <- length(gamma_vals)
  weights_unknown <- array(0, dim = c(n_gamma, max_days, max_days + 1))

  prob_cache <- new.env(hash = TRUE, parent = emptyenv())
  filter_cache <- new.env(hash = TRUE, parent = emptyenv())

  for (g_idx in seq_along(gamma_vals)) {
    c_val <- gamma_vals[g_idx]

    for (t in seq_len(max_days)) {
      if (t <= c_val) {
        next
      }

      for (s in 0:t) {
        idx <- s + 1
        if (s == t) {
          weights_unknown[g_idx, t, idx] <- 1
          next
        }

        numer <- 0
        denom <- 0

        for (u in c_val:(t - 1)) {
          key_tested <- paste0("tested_and_negative_u", u)
          DY_u <- get0(key_tested, filter_cache, ifnotfound = NULL)
          if (is.null(DY_u)) {
            DY_u <- D[, u] == 1 & Y[, u] == 0
            assign(key_tested, DY_u, envir = filter_cache)
          }

          if (s > 0 && s <= c_val) {
            key_filter <- paste0("symptom_le_c", c_val, "_clearance_c", c_val, "_t", t)
            filter <- get0(key_filter, filter_cache)
            if (is.null(filter)) {
              filter <- (Sympt_t[, t] <= c_val) & (C_t[, t] == c_val)
              assign(key_filter, filter, envir = filter_cache)
            }
          } else {
            key_filter <- paste0("symptom_s", s, "_clearance_c", c_val, "_t", t)
            filter <- get0(key_filter, filter_cache)
            if (is.null(filter)) {
              filter <- (Sympt_t[, t] == s) & (C_t[, t] == c_val)
              assign(key_filter, filter, envir = filter_cache)
            }
          }

          key_num <- paste0("num_prob_", key_filter, "_u", u)
          key_den <- paste0("den_prob_", key_filter, "_u", u)

          if (!exists(key_num, envir = prob_cache, inherits = FALSE)) {
            if (u == c_val) {
              if (!any(filter)) {
                num_prob <- 0
                den_prob <- 0
              } else {
                num_prob <- mean((Z_next[, u + 1] == t)[filter], na.rm = TRUE)
                den_prob <- mean((Z_next[, u + 1] >= t)[filter], na.rm = TRUE)
              }
            } else {
              idx_filter <- DY_u & filter
              if (!any(idx_filter)) {
                num_prob <- 0
                den_prob <- 0
              } else {
                num_prob <- spec * mean((Z_next[, u + 1] == t)[idx_filter], na.rm = TRUE)
                den_prob <- spec * mean((Z_next[, u + 1] >= t)[idx_filter], na.rm = TRUE)
              }
            }
            assign(key_num, num_prob, envir = prob_cache)
            assign(key_den, den_prob, envir = prob_cache)
          } else {
            num_prob <- get(key_num, envir = prob_cache, inherits = FALSE)
            den_prob <- get(key_den, envir = prob_cache, inherits = FALSE)
          }

          prev_weight <- if (u == c_val) 1 else weights_unknown[g_idx, u, idx]
          if (!is.finite(prev_weight)) {
            prev_weight <- 0
          }

          multiplier <- if (prev_weight == 0) 0 else 1 / prev_weight
          numer <- numer + num_prob * multiplier
          denom <- denom + den_prob * multiplier
        }

        weights_unknown[g_idx, t, idx] <- if (numer == 0) 0 else denom / numer
      }
    }
  }

  W_est_unknown <- numeric(max_days)

  for (t in seq_len(max_days)) {
    sum_result <- 0
    cs_pairs <- unique(cbind(C_t[, t], Sympt_t[, t]))
    colnames(cs_pairs) <- c("c", "s")

    for (k in seq_len(nrow(cs_pairs))) {
      c_val <- cs_pairs[k, "c"]
      s_val <- cs_pairs[k, "s"]
      g_idx <- match(c_val, gamma_vals)

      if (is.na(g_idx)) {
        next
      }

      weight <- weights_unknown[g_idx, t, s_val + 1]
      if (is.null(weight) || is.na(weight)) {
        weight <- 0
      }

      if (weight == 0) {
        sum_result <- sum_result + sum(
          Sympt_t[, t] == s_val & C_t[, t] == c_val & R[, t] == 0
        )
      } else {
        sum_result <- sum_result +
          weight * sum(
            D[, t] == 1 & Y[, t] == 0 & Sympt_t[, t] == s_val & C_t[, t] == c_val
          ) / (sens + spec - 1) -
          (1 - sens) / (sens + spec - 1) * weight * sum(
            D[, t] == 1 & Sympt_t[, t] == s_val & C_t[, t] == c_val
          )
      }
    }

    W_est_unknown[t] <- sum_result
  }

  (nrow(D) - colSums(R) - W_est_unknown) / (nrow(D) - colSums(R))
}

#' Horvitz-Thompson Prevalence with Block Jackknife Confidence Intervals
#'
#' Computes the daily Horvitz-Thompson prevalence estimate and, by default,
#' 20-block jackknife standard errors and confidence intervals.
#'
#' @param prepared_data Optional preprocessed data list containing `D`, `Y`,
#'   `R`, `Z_next`, `C_t`, `Sympt_t`, and `gamma`. Lists using the legacy
#'   names `Z.next`, `C.t`, and `Sympt.t` are also accepted.
#' @param D Binary testing indicator matrix.
#' @param Y Binary test-result matrix.
#' @param R Logical removed-state matrix.
#' @param Z.next Legacy alias for `Z_next`.
#' @param C.t Legacy alias for `C_t`.
#' @param Sympt.t Legacy alias for `Sympt_t`.
#' @param Z_next Matrix of next testing times.
#' @param C_t Matrix of most recent clearance times.
#' @param Sympt_t Matrix of most recent symptom times.
#' @param gamma Clearance support values, typically `unique(as.vector(C_t))`.
#' @param SPEC Legacy alias for `spec`.
#' @param SENS Legacy alias for `sens`.
#' @param max.days Optional legacy alias for the number of study days.
#' @param spec Test specificity.
#' @param sens Test sensitivity.
#' @param n_blocks Number of delete-a-group jackknife blocks. Defaults to `20`.
#' @param conf_level Confidence level for the interval.
#' @param seed Random seed used for jackknife block assignment.
#' @param n_cores Number of CPU cores used by [block_jackknife()].
#' @param days Optional day labels. Defaults to `prepared_data$days` or
#'   `seq_len(ncol(D))`.
#'
#' @return A list with `results`, `estimate`, `se`, `lower`, `upper`, the
#'   full `jackknife` output, and the resolved matrix inputs used internally.
#' @export
#'
#' @examples
#' sim <- generate_surveillance_data(
#'   n_clusters = 4,
#'   n_per_cluster = 2,
#'   max_days = 6,
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
#'     interval.length = 6,
#'     max.hazard = 0.1,
#'     min.hazard = 0.02
#'   ),
#'   recovery_allowed = FALSE,
#'   recovery_days = 30,
#'   clearance_days = 3,
#'   sens = 0.832,
#'   spec = 0.992,
#'   seed = 1
#' )
#' Z_next <- matrix(NA_real_, nrow = nrow(sim$D), ncol = ncol(sim$D))
#' Z_next[, ncol(sim$D)] <- ifelse(sim$D[, ncol(sim$D)] == 1, ncol(sim$D), Inf)
#' for (t in seq.int(ncol(sim$D) - 1, 1)) {
#'   Z_next[, t] <- ifelse(sim$D[, t] == 1, t, Z_next[, t + 1])
#' }
#' Sympt_t <- matrix(0, nrow = nrow(sim$D), ncol = ncol(sim$D))
#' for (i in seq_len(nrow(sim$D))) {
#'   last_time <- 0
#'   for (t in seq_len(ncol(sim$D))) {
#'     if (sim$symptomatic.mat[i, t] | sim$tracing.mat[i, t]) {
#'       last_time <- t
#'     }
#'     Sympt_t[i, t] <- last_time
#'   }
#' }
#' ht_fit <- ht_prevalence(
#'   D = sim$D,
#'   Y = sim$Y,
#'   R = sim$R,
#'   Z_next = Z_next,
#'   C_t = sim$C.t,
#'   Sympt_t = Sympt_t,
#'   gamma = unique(as.vector(sim$C.t)),
#'   spec = 0.992,
#'   sens = 0.832,
#'   n_blocks = 4,
#'   n_cores = 1
#' )
#' head(ht_fit$results)
ht_prevalence <- function(
    D = NULL,
    Y = NULL,
    R = NULL,
    Z.next = NULL,
    C.t = NULL,
    Sympt.t = NULL,
    gamma = NULL,
    SPEC = NULL,
    SENS = NULL,
    max.days = NULL,
    Z_next = NULL,
    C_t = NULL,
    Sympt_t = NULL,
    spec = NULL,
    sens = NULL,
    prepared_data = NULL,
    n_blocks = 20,
    conf_level = 0.95,
    seed = 12345,
    n_cores = safe_mc_cores(),
    days = NULL) {
  if (is.null(Z_next)) {
    Z_next <- Z.next
  }
  if (is.null(C_t)) {
    C_t <- C.t
  }
  if (is.null(Sympt_t)) {
    Sympt_t <- Sympt.t
  }
  if (is.null(spec)) {
    spec <- SPEC
  }
  if (is.null(sens)) {
    sens <- SENS
  }

  if (is.null(prepared_data) &&
      is.list(D) &&
      is.null(Y) &&
      is.null(R) &&
      is.null(Z_next) &&
      is.null(C_t) &&
      is.null(Sympt_t) &&
      is.null(gamma)) {
    prepared_data <- D
    D <- NULL
  }

  if (is.null(spec) || is.null(sens)) {
    stop("Provide both `spec`/`sens` or the legacy aliases `SPEC`/`SENS`.")
  }

  if (!is.null(max.days)) {
    if (!is.null(D) && ncol(D) != max.days) {
      stop("`max.days` must match `ncol(D)` when both are supplied.")
    }
    if (is.null(days)) {
      days <- seq_len(max.days)
    }
  }

  inputs <- resolve_ht_inputs(
    prepared_data = prepared_data,
    D = D,
    Y = Y,
    R = R,
    Z_next = Z_next,
    C_t = C_t,
    Sympt_t = Sympt_t,
    gamma = gamma,
    days = days
  )

  data_list <- inputs[c("D", "Y", "R", "Z_next", "C_t", "Sympt_t", "gamma")]
  max_days <- ncol(inputs$D)
  test_counts <- unname(colSums(inputs$D))
  no_test_days <- unname(which(test_counts == 0))

  jk_out <- block_jackknife(
    estimator_fun = ht_prevalence_core,
    data_list = data_list,
    n_blocks = n_blocks,
    conf_level = conf_level,
    seed = seed,
    n_cores = n_cores,
    spec = spec,
    sens = sens,
    max_days = max_days
  )

  if (length(no_test_days) > 0) {
    jk_out$estimate[no_test_days] <- NA_real_
    jk_out$se[no_test_days] <- NA_real_
    jk_out$lower[no_test_days] <- NA_real_
    jk_out$upper[no_test_days] <- NA_real_

    if (!is.null(jk_out$theta_leave_block_out)) {
      jk_out$theta_leave_block_out[, no_test_days] <- NA_real_
    }
  }

  results <- data.frame(
    day = inputs$days,
    estimate = jk_out$estimate,
    se = jk_out$se,
    lower = jk_out$lower,
    upper = jk_out$upper
  )

  list(
    results = results,
    estimate = jk_out$estimate,
    se = jk_out$se,
    lower = jk_out$lower,
    upper = jk_out$upper,
    jackknife = jk_out,
    inputs = inputs,
    days = inputs$days,
    low_test_days = if (!is.null(prepared_data$low_test_days)) unname(prepared_data$low_test_days) else integer(0),
    no_test_days = no_test_days,
    test_counts = test_counts
  )
}

#' Delete-a-Group Block Jackknife
#'
#' Computes block jackknife estimates, standard errors, and confidence
#' intervals for an arbitrary matrix-based estimator.
#'
#' @param estimator_fun Estimator function.
#' @param data_list Named list of inputs to `estimator_fun`.
#' @param n_blocks Number of jackknife blocks.
#' @param block_ids Optional vector of block memberships.
#' @param conf_level Confidence level.
#' @param seed Random seed used when `block_ids` is not supplied.
#' @param n_cores Number of cores. On Windows or when `n_cores <= 1`, the
#'   function falls back to sequential evaluation.
#' @param ... Additional arguments passed to `estimator_fun`.
#'
#' @return A list with the point estimate, standard errors, confidence
#'   intervals, leave-one-block-out estimates, and block metadata.
#' @export
block_jackknife <- function(
    estimator_fun,
    data_list,
    n_blocks = 10,
    block_ids = NULL,
    conf_level = 0.95,
    seed = 1,
    n_cores = safe_mc_cores(),
    ...) {
  D <- data_list$D
  n <- nrow(D)

  if (is.null(block_ids)) {
    n_blocks <- min(n_blocks, n)
    set.seed(seed)
    perm <- sample.int(n)
    block_ids <- integer(n)
    block_labels <- rep(seq_len(n_blocks), length.out = n)
    block_ids[perm] <- block_labels
  } else {
    if (length(block_ids) != n) {
      stop("Length of block_ids must match the number of rows in D.")
    }
    block_ids <- as.integer(as.factor(block_ids))
    n_blocks <- length(unique(block_ids))
  }

  if (n_blocks < 2) {
    stop("At least two jackknife blocks are required.")
  }

  theta_hat <- do.call(estimator_fun, c(data_list, list(...)))

  worker <- function(b) {
    keep <- block_ids != b
    data_b <- subset_data_list_by_keep(data_list, keep)
    do.call(estimator_fun, c(data_b, list(...)))
  }

  if (.Platform$OS.type == "windows" || n_cores <= 1) {
    theta_list <- lapply(seq_len(n_blocks), worker)
  } else {
    theta_list <- parallel::mclapply(
      X = seq_len(n_blocks),
      FUN = worker,
      mc.cores = n_cores
    )
  }

  theta_j <- do.call(rbind, theta_list)
  theta_bar <- colMeans(theta_j, na.rm = TRUE)
  diffs <- sweep(theta_j, 2, theta_bar)
  var_jk <- (n_blocks - 1) / n_blocks * colSums(diffs^2, na.rm = TRUE)
  se_jk <- sqrt(var_jk)

  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  list(
    estimate = theta_hat,
    se = se_jk,
    lower = theta_hat - z * se_jk,
    upper = theta_hat + z * se_jk,
    theta_leave_block_out = theta_j,
    block_ids = block_ids,
    n_blocks = n_blocks,
    seed = seed
  )
}

#' Joint Kalman Filter for Daily Prevalence
#'
#' Fits the local linear trend state-space model with fixed process variances.
#'
#' @param y Numeric vector of noisy prevalence estimates.
#' @param r_t Numeric vector of observation variances.
#' @param q_level Level process variance.
#' @param q_slope Slope process variance.
#' @param init_level Optional initial level.
#' @param init_slope Optional initial slope.
#' @param init_cov Optional initial state covariance matrix.
#'
#' @return A list containing predicted and filtered states and covariances.
#' @export
kalman_filter_joint <- function(
    y,
    r_t,
    q_level,
    q_slope,
    init_level = NULL,
    init_slope = 0,
    init_cov = NULL) {
  n <- length(y)
  obs_idx <- which(is.finite(y) & is.finite(r_t) & r_t >= 0)

  if (length(obs_idx) == 0) {
    stop("No observed y values are available for Kalman filtering.")
  }

  first_obs <- obs_idx[1]
  if (is.null(init_level)) {
    init_level <- y[first_obs]
  }
  if (is.null(init_cov)) {
    init_cov <- diag(c(r_t[first_obs] + q_level + q_slope, max(q_slope, 1e-6)), nrow = 2)
  }

  F_mat <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
  H_mat <- matrix(c(1, 0), nrow = 1)
  Q_mat <- matrix(c(q_level, 0, 0, q_slope), nrow = 2, byrow = TRUE)
  I_mat <- diag(2)

  alpha_pred <- matrix(NA_real_, nrow = 2, ncol = n)
  P_pred <- array(NA_real_, dim = c(2, 2, n))
  alpha_filt <- matrix(NA_real_, nrow = 2, ncol = n)
  P_filt <- array(NA_real_, dim = c(2, 2, n))

  alpha_pred[, 1] <- c(init_level, init_slope)
  P_pred[, , 1] <- init_cov

  for (t in seq_len(n)) {
    if (t > 1) {
      alpha_pred[, t] <- F_mat %*% alpha_filt[, t - 1]
      P_pred[, , t] <- F_mat %*% P_filt[, , t - 1] %*% t(F_mat) + Q_mat
    }

    observed_t <- is.finite(y[t]) & is.finite(r_t[t]) & r_t[t] >= 0

    if (observed_t) {
      S_t <- as.numeric(H_mat %*% P_pred[, , t] %*% t(H_mat) + r_t[t])
      K_vec <- P_pred[, , t] %*% t(H_mat) / S_t
      innov <- y[t] - as.numeric(H_mat %*% alpha_pred[, t])

      alpha_filt[, t] <- alpha_pred[, t] + as.numeric(K_vec) * innov
      P_filt[, , t] <- (I_mat - K_vec %*% H_mat) %*% P_pred[, , t]
    } else {
      alpha_filt[, t] <- alpha_pred[, t]
      P_filt[, , t] <- P_pred[, , t]
    }
  }

  list(
    alpha_pred = alpha_pred,
    alpha_filt = alpha_filt,
    level_pred = alpha_pred[1, ],
    slope_pred = alpha_pred[2, ],
    level_filt = alpha_filt[1, ],
    level_filt_trunc = pmin(1, pmax(0, alpha_filt[1, ])),
    slope_filt = alpha_filt[2, ],
    P_pred = P_pred,
    P_filt = P_filt
  )
}

#' Joint Kalman Smoother for Daily Prevalence
#'
#' Applies fixed-interval smoothing to the local linear trend model.
#'
#' @inheritParams kalman_filter_joint
#'
#' @return A list containing the filtered fit and smoothed state trajectories.
#' @export
kalman_smoother_joint <- function(
    y,
    r_t,
    q_level,
    q_slope,
    init_level = NULL,
    init_slope = 0,
    init_cov = NULL) {
  filt <- kalman_filter_joint(
    y = y,
    r_t = r_t,
    q_level = q_level,
    q_slope = q_slope,
    init_level = init_level,
    init_slope = init_slope,
    init_cov = init_cov
  )

  n <- length(y)
  F_mat <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)

  alpha_filt <- filt$alpha_filt
  alpha_pred <- filt$alpha_pred

  alpha_smooth <- matrix(NA_real_, nrow = 2, ncol = n)
  P_smooth <- array(NA_real_, dim = c(2, 2, n))

  alpha_smooth[, n] <- alpha_filt[, n]
  P_smooth[, , n] <- filt$P_filt[, , n]

  if (n >= 2) {
    for (t in seq.int(n - 1, 1)) {
      J_t <- filt$P_filt[, , t] %*% t(F_mat) %*% solve(filt$P_pred[, , t + 1])
      alpha_smooth[, t] <- alpha_filt[, t] +
        J_t %*% (alpha_smooth[, t + 1] - alpha_pred[, t + 1])
      P_smooth[, , t] <- filt$P_filt[, , t] +
        J_t %*% (P_smooth[, , t + 1] - filt$P_pred[, , t + 1]) %*% t(J_t)
    }
  }

  list(
    filtered_fit = filt,
    alpha_smooth = alpha_smooth,
    level_smooth = alpha_smooth[1, ],
    level_smooth_trunc = pmin(1, pmax(0, alpha_smooth[1, ])),
    slope_smooth = alpha_smooth[2, ],
    P_smooth = P_smooth
  )
}

get_model_q <- function(model, free_log_q) {
  if (model == "Level-Only") {
    c(Q_level = exp(free_log_q[1]), Q_slope = 0)
  } else if (model == "Slope-Only") {
    c(Q_level = 0, Q_slope = exp(free_log_q[1]))
  } else if (model == "Joint") {
    c(Q_level = exp(free_log_q[1]), Q_slope = exp(free_log_q[2]))
  } else {
    stop("Unknown model: ", model)
  }
}

negloglik_joint_q <- function(
    free_log_q,
    y,
    r_t,
    model,
    init_level = NULL,
    init_slope = 0,
    init_cov = NULL) {
  if (length(y) < 2) {
    return(Inf)
  }

  q_vals <- get_model_q(model, free_log_q)

  fit <- kalman_filter_joint(
    y = y,
    r_t = r_t,
    q_level = q_vals["Q_level"],
    q_slope = q_vals["Q_slope"],
    init_level = init_level,
    init_slope = init_slope,
    init_cov = init_cov
  )

  idx <- which(
    seq_along(y) >= 2 &
      is.finite(y) &
      is.finite(r_t) &
      r_t >= 0 &
      is.finite(fit$level_pred)
  )

  if (length(idx) == 0) {
    return(Inf)
  }

  S_t <- vapply(idx, function(t) as.numeric(fit$P_pred[1, 1, t] + r_t[t]), numeric(1))
  nu_t <- y[idx] - fit$level_pred[idx]

  if (any(!is.finite(S_t)) || any(S_t <= 0)) {
    return(Inf)
  }

  0.5 * sum(log(2 * pi) + log(S_t) + (nu_t^2) / S_t)
}

#' Estimate Kalman Process Variances
#'
#' Estimates process variances by maximizing the Gaussian innovation
#' likelihood.
#'
#' @param y Numeric vector of noisy prevalence estimates.
#' @param r_t Numeric vector of observation variances.
#' @param model One of `"Joint"`, `"Level-Only"`, or `"Slope-Only"`.
#' @param init_level Optional initial level.
#' @param init_slope Optional initial slope.
#' @param init_cov Optional initial covariance matrix.
#' @param q_level_start Starting value for the level variance.
#' @param q_slope_start Starting value for the slope variance.
#' @param q_lower Lower bound used during optimization.
#' @param q_upper Upper bound used during optimization.
#'
#' @return A list with estimated process variances and optimization details.
#' @export
estimate_kalman_q <- function(
    y,
    r_t,
    model = c("Joint", "Level-Only", "Slope-Only"),
    init_level = NULL,
    init_slope = 0,
    init_cov = NULL,
    q_level_start = 1e-6,
    q_slope_start = 1e-5,
    q_lower = 1e-10,
    q_upper = 1e-2) {
  model <- match.arg(model)
  obs_idx <- which(is.finite(y) & is.finite(r_t) & r_t >= 0)

  if (length(obs_idx) < 2) {
    stop("At least two observed days are needed to estimate process variances.")
  }

  if (is.null(init_cov)) {
    init_cov <- diag(
      c(r_t[obs_idx[1]] + q_level_start + q_slope_start, max(q_slope_start, 1e-6)),
      nrow = 2
    )
  }
  if (is.null(init_level)) {
    init_level <- y[obs_idx[1]]
  }

  if (model %in% c("Level-Only", "Slope-Only")) {
    opt <- stats::optimize(
      f = function(par, ...) negloglik_joint_q(par, ..., model = model),
      interval = c(log(q_lower), log(q_upper)),
      y = y,
      r_t = r_t,
      init_level = init_level,
      init_slope = init_slope,
      init_cov = init_cov
    )
    q_hat <- get_model_q(model, opt$minimum)
    return(list(
      Q_level_hat = unname(q_hat["Q_level"]),
      Q_slope_hat = unname(q_hat["Q_slope"]),
      negloglik = opt$objective,
      convergence = TRUE
    ))
  }

  opt <- stats::optim(
    par = log(c(q_level_start, q_slope_start)),
    fn = negloglik_joint_q,
    y = y,
    r_t = r_t,
    model = model,
    init_level = init_level,
    init_slope = init_slope,
    init_cov = init_cov,
    method = "L-BFGS-B",
    lower = rep(log(q_lower), 2),
    upper = rep(log(q_upper), 2)
  )

  q_hat <- get_model_q(model, opt$par)

  list(
    Q_level_hat = unname(q_hat["Q_level"]),
    Q_slope_hat = unname(q_hat["Q_slope"]),
    negloglik = opt$value,
    convergence = opt$convergence == 0
  )
}

#' Fit a Kalman Prevalence Model
#'
#' Fits one of the manuscript's Kalman model variants and returns filtered and
#' smoothed prevalence estimates with pointwise normal intervals.
#'
#' @param y Numeric vector of noisy prevalence estimates.
#' @param r_t Numeric vector of observation variances.
#' @param model One of `"Joint"`, `"Level-Only"`, or `"Slope-Only"`.
#' @param init_level Optional initial level. Defaults to the first observed
#'   value.
#' @param init_slope Optional initial slope.
#' @param q_level_start Starting value for the level variance.
#' @param q_slope_start Starting value for the slope variance.
#' @param q_lower Lower bound used during optimization.
#' @param q_upper Upper bound used during optimization.
#'
#' @return A list containing estimated process variances and filtered/smoothed
#'   prevalence trajectories.
#' @keywords internal
#' @noRd
#'
#' @examples
#' y <- c(0.01, 0.02, NA, 0.03, 0.025)
#' r_t <- c(0.0001, 0.0001, NA, 0.0002, 0.00015)
#' fit_kalman_model(y, r_t, model = "Joint")
fit_kalman_model_single <- function(
    y,
    r_t,
    model = c("Joint", "Level-Only", "Slope-Only"),
    init_level = NULL,
    init_slope = 0,
    q_level_start = 1e-6,
    q_slope_start = 1e-5,
    q_lower = 1e-10,
    q_upper = 1e-2) {
  model <- match.arg(model)

  if (is.null(init_level)) {
    init_idx <- which(is.finite(y) & is.finite(r_t) & r_t >= 0)[1]
    if (is.na(init_idx)) {
      stop("No valid observations are available to initialize the model.")
    }
    init_level <- y[init_idx]
  }

  q_fit <- estimate_kalman_q(
    y = y,
    r_t = r_t,
    model = model,
    init_level = init_level,
    init_slope = init_slope,
    q_level_start = q_level_start,
    q_slope_start = q_slope_start,
    q_lower = q_lower,
    q_upper = q_upper
  )

  ks_fit <- kalman_smoother_joint(
    y = y,
    r_t = r_t,
    q_level = q_fit$Q_level_hat,
    q_slope = q_fit$Q_slope_hat,
    init_level = init_level,
    init_slope = init_slope
  )
  kf_fit <- ks_fit$filtered_fit

  filtered_se <- sqrt(pmax(0, vapply(
    seq_along(kf_fit$level_filt),
    function(t) kf_fit$P_filt[1, 1, t],
    numeric(1)
  )))
  filtered_lower <- pmax(0, kf_fit$level_filt_trunc - 1.96 * filtered_se)
  filtered_upper <- pmin(1, kf_fit$level_filt_trunc + 1.96 * filtered_se)

  smoothed_se <- sqrt(pmax(0, vapply(
    seq_along(ks_fit$level_smooth),
    function(t) ks_fit$P_smooth[1, 1, t],
    numeric(1)
  )))
  smoothed_lower <- pmax(0, ks_fit$level_smooth_trunc - 1.96 * smoothed_se)
  smoothed_upper <- pmin(1, ks_fit$level_smooth_trunc + 1.96 * smoothed_se)

  list(
    label = model,
    Q_level = q_fit$Q_level_hat,
    Q_slope = q_fit$Q_slope_hat,
    negloglik = q_fit$negloglik,
    filtered = kf_fit$level_filt_trunc,
    filtered_se = filtered_se,
    filtered_lower = filtered_lower,
    filtered_upper = filtered_upper,
    smoothed = ks_fit$level_smooth_trunc,
    smoothed_se = smoothed_se,
    smoothed_lower = smoothed_lower,
    smoothed_upper = smoothed_upper
  )
}

build_joint_ci6_results <- function(
    fit,
    state,
    days,
    ht_estimate,
    ht_se,
    ht_jackknife,
    conf_level,
    low_test_days = integer(0),
    init_slope = 0,
    q_level_start = 1e-6,
    q_slope_start = 1e-5,
    q_lower = 1e-10,
    q_upper = 1e-2) {
  if (state == "filter") {
    model_estimate <- fit$filtered
    model_se <- fit$filtered_se
  } else {
    model_estimate <- fit$smoothed
    model_se <- fit$smoothed_se
  }

  if (is.null(ht_jackknife) || nrow(ht_jackknife) < 2) {
    stop(
      "CI-6 requires `x` from ht_prevalence() so leave-one-block-out jackknife ",
      "estimates are available."
    )
  }

  r_t <- ht_se^2
  if (length(low_test_days) > 0) {
    r_t[low_test_days] <- NA_real_
  }

  biased_jk <- t(vapply(seq_len(nrow(ht_jackknife)), function(i) {
    y_i <- ht_jackknife[i, ]
    if (length(low_test_days) > 0) {
      y_i[low_test_days] <- NA_real_
    }

    valid_obs <- which(is.finite(y_i) & is.finite(r_t) & r_t >= 0)
    if (length(valid_obs) < 2) {
      return(rep(NA_real_, length(days)))
    }

    fit_i <- tryCatch(
      fit_kalman_model_single(
        y = y_i,
        r_t = r_t,
        model = "Joint",
        init_slope = init_slope,
        q_level_start = q_level_start,
        q_slope_start = q_slope_start,
        q_lower = q_lower,
        q_upper = q_upper
      ),
      error = function(...) NULL
    )

    if (is.null(fit_i)) {
      return(rep(NA_real_, length(days)))
    }

    if (state == "filter") fit_i$filtered else fit_i$smoothed
  }, numeric(length(days))))

  rho <- vapply(seq_along(days), function(j) {
    keep <- is.finite(ht_jackknife[, j]) & is.finite(biased_jk[, j])
    if (sum(keep) < 2) {
      return(NA_real_)
    }

    x_j <- ht_jackknife[keep, j]
    y_j <- biased_jk[keep, j]
    if (stats::sd(x_j) == 0 || stats::sd(y_j) == 0) {
      return(NA_real_)
    }

    stats::cor(x_j, y_j)
  }, numeric(1))

  ci6 <- compute_ci6_weighted_results(
    ht_estimate = ht_estimate,
    model_estimate = model_estimate,
    ht_se = ht_se,
    model_se = model_se,
    rho = rho,
    conf_level = conf_level
  )

  data.frame(
    day = days,
    model = "Joint",
    state = state,
    ci_method = "convex_combination",
    estimate = ci6$estimate,
    se = ci6$se,
    critical_value = ci6$critical_value,
    weight = ci6$weight,
    rho = ci6$rho,
    combo_se = ci6$combo_se,
    lower = ci6$lower,
    upper = ci6$upper
  )
}

#' Fit the Joint Kalman Prevalence Model
#'
#' Fits the joint local linear trend Kalman model and returns daily estimates
#' with one or more confidence-interval options for either the filter or the
#' smoother.
#'
#' @param x Optional output from [ht_prevalence()] or a data frame containing
#'   `estimate` and `se` columns.
#' @param y Optional numeric vector of prevalence estimates.
#' @param r_t Optional numeric vector of observation variances.
#' @param state One of `"filter"`, `"smoother"`, or `"both"`.
#' @param ci_method Confidence-interval option. Choose from
#'   `"model_based"` (default), `"ht_se"`, `"ht_se_adjusted"`, or `"all"`.
#' @param conf_level Confidence level for the interval.
#' @param low_test_days Optional integer vector of days to treat as missing
#'   before fitting the state-space model.
#' @param init_level Optional initial level. Defaults to the first observed
#'   estimate after masking missing days.
#' @param init_slope Optional initial slope.
#' @param q_level_start Starting value for the level variance.
#' @param q_slope_start Starting value for the slope variance.
#' @param q_lower Lower bound used during optimization.
#' @param q_upper Upper bound used during optimization.
#' @param days Optional day labels.
#'
#' @return A list with `results`, one row per day and CI option, and
#'   `joint_fit`, the full fitted joint model object.
#' @export
#'
#' @examples
#' y <- c(0.01, 0.02, NA, 0.03, 0.025)
#' r_t <- c(0.0001, 0.0001, NA, 0.0002, 0.00015)
#' fit_kalman_model(
#'   y = y,
#'   r_t = r_t,
#'   state = "filter",
#'   ci_method = "all"
#' )
fit_kalman_model <- function(
    y = NULL,
    r_t = NULL,
    x = NULL,
    state = c("filter", "smoother", "both"),
    ci_method = c("model_based", "ht_se", "ht_se_adjusted", "all"),
    conf_level = 0.95,
    low_test_days = NULL,
    init_level = NULL,
    init_slope = 0,
    q_level_start = 1e-6,
    q_slope_start = 1e-5,
    q_lower = 1e-10,
    q_upper = 1e-2,
    days = NULL) {
  if (is.null(x) && (is.list(y) || is.data.frame(y)) && is.null(r_t)) {
    x <- y
    y <- NULL
  }

  state <- match.arg(state)
  ci_method <- match.arg(ci_method)

  inputs <- extract_kalman_inputs(x = x, y = y, r_t = r_t, days = days)

  if (is.null(low_test_days) && is.list(x) && !is.null(x$low_test_days)) {
    low_test_days <- x$low_test_days
  }

  y_fit <- inputs$y
  r_t_fit <- inputs$r_t

  if (!is.null(low_test_days) && length(low_test_days) > 0) {
    y_fit[low_test_days] <- NA_real_
    r_t_fit[low_test_days] <- NA_real_
  }

  if (is.null(init_level)) {
    init_idx <- which(is.finite(y_fit) & is.finite(r_t_fit) & r_t_fit >= 0)[1]
    if (is.na(init_idx)) {
      stop("No valid observations are available to initialize the model.")
    }
    init_level <- y_fit[init_idx]
  }

  joint_fit <- fit_kalman_model_single(
    y = y_fit,
    r_t = r_t_fit,
    model = "Joint",
    init_level = init_level,
    init_slope = init_slope,
    q_level_start = q_level_start,
    q_slope_start = q_slope_start,
    q_lower = q_lower,
    q_upper = q_upper
  )

  ci_methods <- if (ci_method == "all") {
    c("model_based", "ht_se", "ht_se_adjusted")
  } else {
    ci_method
  }

  if (state == "both") {
    results <- do.call(
      rbind,
      lapply(ci_methods, function(method) {
        rbind(
          build_joint_ci_results(joint_fit, "filter", method, inputs$days, inputs$ht_se, conf_level),
          build_joint_ci_results(joint_fit, "smoother", method, inputs$days, inputs$ht_se, conf_level)
        )
      })
    )
  } else {
    results <- do.call(
      rbind,
      lapply(ci_methods, function(method) {
        build_joint_ci_results(joint_fit, state, method, inputs$days, inputs$ht_se, conf_level)
      })
    )
  }

  rownames(results) <- NULL

  list(
    results = results,
    joint_fit = joint_fit,
    state = state,
    ci_method = ci_methods,
    low_test_days = if (is.null(low_test_days)) integer(0) else low_test_days
  )
}

#' Fit the CI-6 Convex-Combination Mixture
#'
#' Combines the Horvitz-Thompson estimator with the joint Kalman filter or
#' smoother using the CI-6 convex-combination construction.
#'
#' @param x Output from [ht_prevalence()]. This function requires the
#'   leave-one-block-out jackknife estimates stored in `x$jackknife`.
#' @param state One of `"filter"`, `"smoother"`, or `"both"`.
#' @param conf_level Confidence level for the interval.
#' @param low_test_days Optional integer vector of days to treat as missing
#'   before fitting the state-space model. Defaults to `x$low_test_days`.
#' @param init_level Optional initial level. Defaults to the first observed
#'   estimate after masking missing days.
#' @param init_slope Optional initial slope.
#' @param q_level_start Starting value for the level variance.
#' @param q_slope_start Starting value for the slope variance.
#' @param q_lower Lower bound used during optimization.
#' @param q_upper Upper bound used during optimization.
#' @param days Optional day labels.
#'
#' @return A list with `results`, one row per day and requested state, and
#'   `joint_fit`, the fitted joint Kalman model. The `results` data frame
#'   includes the CI-6 mixture estimate, weight, correlation estimate, and
#'   combination standard error.
#' @export
#'
#' @examples
#' \dontrun{
#' # D: testing indicator matrix
#' # Y: test-result matrix
#' # R: removed/exempt-state matrix
#' # Z.next: next testing-time matrix
#' # C.t: most recent clearance-time matrix
#' # Sympt.t: most recent symptom/contact-time matrix
#' # gamma: support values of C.t
#'
#' ht_fit <- ht_prevalence(
#'   D = D,
#'   Y = Y,
#'   R = R,
#'   Z.next = Z.next,
#'   C.t = C.t,
#'   Sympt.t = Sympt.t,
#'   gamma = gamma,
#'   SPEC = 1,
#'   SENS = 0.832,
#'   max.days = ncol(D)
#' )
#'
#' fit_mixture(
#'   x = ht_fit,
#'   state = "smoother"
#' )
#' }
fit_mixture <- function(
    x,
    state = c("filter", "smoother", "both"),
    conf_level = 0.95,
    low_test_days = NULL,
    init_level = NULL,
    init_slope = 0,
    q_level_start = 1e-6,
    q_slope_start = 1e-5,
    q_lower = 1e-10,
    q_upper = 1e-2,
    days = NULL) {
  if (missing(x) || !is.list(x) || is.null(x$jackknife$theta_leave_block_out)) {
    stop(
      "`fit_mixture()` requires the output from ht_prevalence() so ",
      "leave-one-block-out jackknife estimates are available."
    )
  }

  state <- match.arg(state)
  inputs <- extract_kalman_inputs(x = x, days = days)

  if (is.null(low_test_days) && !is.null(x$low_test_days)) {
    low_test_days <- x$low_test_days
  }

  y_fit <- inputs$y
  r_t_fit <- inputs$r_t
  ht_estimate <- inputs$y
  ht_se <- inputs$ht_se

  if (!is.null(low_test_days) && length(low_test_days) > 0) {
    y_fit[low_test_days] <- NA_real_
    r_t_fit[low_test_days] <- NA_real_
    ht_estimate[low_test_days] <- NA_real_
    ht_se[low_test_days] <- NA_real_
  }

  if (is.null(init_level)) {
    init_idx <- which(is.finite(y_fit) & is.finite(r_t_fit) & r_t_fit >= 0)[1]
    if (is.na(init_idx)) {
      stop("No valid observations are available to initialize the model.")
    }
    init_level <- y_fit[init_idx]
  }

  joint_fit <- fit_kalman_model_single(
    y = y_fit,
    r_t = r_t_fit,
    model = "Joint",
    init_level = init_level,
    init_slope = init_slope,
    q_level_start = q_level_start,
    q_slope_start = q_slope_start,
    q_lower = q_lower,
    q_upper = q_upper
  )

  build_results_for_state <- function(state_name) {
    build_joint_ci6_results(
      fit = joint_fit,
      state = state_name,
      days = inputs$days,
      ht_estimate = ht_estimate,
      ht_se = ht_se,
      ht_jackknife = x$jackknife$theta_leave_block_out,
      conf_level = conf_level,
      low_test_days = if (is.null(low_test_days)) integer(0) else low_test_days,
      init_slope = init_slope,
      q_level_start = q_level_start,
      q_slope_start = q_slope_start,
      q_lower = q_lower,
      q_upper = q_upper
    )
  }

  results <- if (state == "both") {
    rbind(
      build_results_for_state("filter"),
      build_results_for_state("smoother")
    )
  } else {
    build_results_for_state(state)
  }

  rownames(results) <- NULL

  list(
    results = results,
    joint_fit = joint_fit,
    state = state,
    low_test_days = if (is.null(low_test_days)) integer(0) else low_test_days
  )
}

