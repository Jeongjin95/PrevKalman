#' OSU Fall 2020 SARS-CoV-2 Surveillance Data
#'
#' Long-form repeated-testing data used in the manuscript's real-data example.
#'
#' @format A data frame with 137779 rows and 5 variables:
#' \describe{
#'   \item{name_n}{Integer participant identifier.}
#'   \item{test_day}{Integer study day.}
#'   \item{provider}{Integer provider/source code.}
#'   \item{result}{Binary test result indicator.}
#'   \item{symp_cont}{Binary symptom/contact-trigger indicator.}
#' }
#'
#' @source Ohio State University Fall 2020 surveillance dataset as provided in
#'   the submission workspace.
"osu_surveillance"
