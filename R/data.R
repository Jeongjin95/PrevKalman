#' Preprocessed OSU Fall 2020 SARS-CoV-2 Surveillance Data
#'
#' Preprocessed repeated-testing matrices used in the manuscript's real-data
#' example. The package ships the prepared object directly, so users can run
#' the HT and Kalman estimators without a separate preprocessing step.
#'
#' @format A named list with 12 elements:
#' \describe{
#'   \item{D}{Binary testing indicator matrix.}
#'   \item{Y}{Binary test-result matrix.}
#'   \item{A_sympt}{Binary symptom/contact indicator matrix.}
#'   \item{R}{Logical removed-state matrix.}
#'   \item{C_t}{Most recent clearance-time matrix.}
#'   \item{Z_next}{Next testing-time matrix.}
#'   \item{Sympt_t}{Most recent symptom/contact-time matrix.}
#'   \item{gamma}{Clearance support values.}
#'   \item{days}{Study-day labels.}
#'   \item{low_test_days}{Days flagged as having low testing volume.}
#'   \item{test_counts}{Number of tests observed on each day.}
#'   \item{total_removed_days}{Clearance plus exemption period used in the
#'   prepared example.}
#' }
#'
#' @source Ohio State University Fall 2020 surveillance dataset as provided in
#'   the submission workspace, preprocessed with `clearance_days = 10`,
#'   `exempt_days = 80`, `symptom_unobserved_days = 1:12`, and
#'   `low_test_threshold = 100`.
"osu_surveillance_preprocessed"
