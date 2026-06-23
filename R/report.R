# 135 characters ######################################################################################################################
#' @title Report the top-K models by posterior probability
#' @description Extracts and formats the \code{top} models with the highest
#' posterior probabilities from a completed MRBMA run. Causal estimates for each reported model are recomputed via
#' \code{beta_summary_cpp} from pre-computed summary-level statistics, giving results numerically identical to those stored in the
#' output object. A unified interface is provided: the function accepts both \code{MRBF} (exhaustive) and \code{mvMR_SSS}
#' (stochastic) output objects
#'
#' @param output     An object of class \code{MRBF} or \code{mvMR_SSS}, as
#'                   returned by \code{MRBMA_exhaustive} or \code{MRBMA_SSS}
#'                   respectively
#' @param top        Positive integer. Number of top models to report
#'                   (default \code{10}). If fewer models were visited, all
#'                   are returned with a warning
#' @param digits     Non-negative integer. Number of decimal places for
#'                   rounding of posterior probabilities and causal estimates
#'                   (default \code{3})
#' @param write_out  Logical. If \code{TRUE}, the result table is written to
#'                   a \code{.csv} file (default \code{FALSE})
#' @param file_name  Character scalar. File path and name (without extension)
#'                   for the CSV output when \code{write_out = TRUE}
#'                   (default \code{"best_models"})
#'
#' @details
#' A unified dispatcher detects the class of the input object automatically.
#' Results summarise the posterior distribution over models as defined in
#' \insertCite{Zuber2020;textual}{MRBMA}.
#'
#' @export
#'
#' @return A \code{data.frame} with \code{top} rows and three columns:
#' \describe{
#'   \item{\code{Exposure Combination}}{Comma-separated names of exposures
#'         included in the model}
#'   \item{\code{Marginal Posterior Probability of Inclusion (MPPI)}}{Posterior
#'         probability of the model, rounded to \code{digits} decimal places}
#'   \item{\code{BMA Causal Estimate}}{Comma-separated causal estimates for
#'         the exposures in the model, rounded to \code{digits} decimal places}
#' }
#' Rows are sorted by \code{Marginal Posterior Probability of Inclusion (MPPI)} in
#' descending order
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{report_MRBMA}}, \code{\link{MRBMA_exhaustive}},
#' \code{\link{MRBMA_SSS}}
#'
#' @examples
#' \dontrun{
#'
#' data(AMD_data)
#' betaX_ivw   <- as.matrix(AMD_data$betaX) / AMD_data$seAMD
#' betaAMD_ivw <- AMD_data$betaAMD / AMD_data$seAMD
#' rs          <- AMD_data$annotate[, 1]
#' rf          <- colnames(AMD_data$betaX)
#'
#' mvmr_input <- new("mvMRInput",
#'                   betaX    = betaX_ivw,
#'                   betaY    = as.matrix(betaAMD_ivw),
#'                   snps     = rs,
#'                   exposure = rf,
#'                   outcome  = "AMD")
#'
#' bma_out <- MRBMA_SSS(mvmr_input, kmin = 1, kmax = 12, max_iter = 100)
#'
#' # Print top 10 models to console
#' report_best_models(bma_out, top = 10, digits = 3)
#'
#' # Write top 10 models to CSV
#' report_best_models(bma_out, top = 10, write_out = TRUE,
#'                    file_name = "amd_best_models")
#'
#' }
report_best_models <- function(output,
                                top       = 10,
                                digits    = 3,
                                write_out = FALSE,
                                file_name = "best_models") {

  # ---- Validate class ----
  if (!is(output, "MRBF") && !is(output, "mvMR_SSS"))
    stop("'output' must be of class 'MRBF' or 'mvMR_SSS'.")

  pp     <- output@pp
  models <- output@tupel
  betaX  <- output@betaX
  betaY  <- output@betaY
  rf     <- output@Exposure
  sigma  <- output@sigma

  # ---- Guard: top cannot exceed number of models ----
  n_top <- min(top, length(pp))
  if (n_top < top)
    warning(sprintf(
      "Only %d models were visited (requested top = %d). Returning all %d.",
      n_top, top, n_top
    ))
  if (n_top == 0)
    stop("No models found in output object.")

  # ---- Sort by posterior probability ----
  sort_obj   <- sort.int(pp, index.return = TRUE, decreasing = TRUE)
  top_models <- models[sort_obj$ix][seq_len(n_top)]
  top_pp     <- pp[sort_obj$ix][seq_len(n_top)]

  # ---- Pre-compute sufficient statistics for beta recomputation ----
  sigma_vec  <- rep(sigma, ncol(betaX))
  XtX        <- t(betaX) %*% betaX
  XtY        <- as.vector(t(betaX) %*% betaY)
  d_Exposure <- ncol(betaX)

  # ---- Build output table ----
  rf_top    <- character(n_top)
  theta_top <- character(n_top)

  for (i in seq_len(n_top)) {
    tupel_i <- as.integer(unlist(strsplit(top_models[i], ",")))

    if (any(tupel_i < 1) || any(tupel_i > d_Exposure))
      stop(sprintf("Model %d contains out-of-range exposure indices.", i))

    rf_top[i] <- paste(rf[tupel_i], collapse = ", ")

    theta_vec  <- beta_summary_cpp(
      XtY        = XtY,
      XtX        = XtX,
      sigma_vec  = sigma_vec,
      gamma_idx  = as.integer(tupel_i),
      d_Exposure = d_Exposure
    )
    theta_top[i] <- paste(round(theta_vec[tupel_i], digits), collapse = ", ")
  }

  out <- data.frame(
    "Exposure Combination"                    = rf_top,
    "Marginal Posterior Probability of Inclusion (MPPI)" = round(top_pp, digits),
    "BMA Causal Estimate"                   = theta_top,
    stringsAsFactors                          = FALSE,
    check.names                               = FALSE
  )

  if (write_out) {
    fname <- if (grepl("\\.csv$", file_name, ignore.case = TRUE)) {
      file_name
    } else {
      paste0(file_name, ".csv")
    }
    write.csv(out, file = fname, row.names = FALSE)
    message("Best models table written to: ", fname)
  }

  return(out)
}


# 135 characters ######################################################################################################################
#' @title Report model-averaged causal effects (MRBMA table)
#' @description Extracts and formats the exposure-level summary of a
#' completed MRBMA run: Marginal Posterior Probabilities of Inclusion (MPPIs) and BMA-averaged causal estimates for the \code{top}
#' exposures. This is the primary output table of the MRBMA method as defined in \insertCite{Zuber2020;textual}{MRBMA}. A unified
#' interface is provided: the function accepts both \code{MRBF} and \code{mvMR_SSS} objects
#'
#' @param output     An object of class \code{MRBF} or \code{mvMR_SSS}, as
#'                   returned by \code{MRBMA_exhaustive} or \code{MRBMA_SSS}
#' @param top        Positive integer. Number of top exposures to report,
#'                   sorted by Marginal Posterior Probability of Inclusion (MPPI) (default \code{10})
#' @param digits     Non-negative integer. Number of decimal places for
#'                   rounding (default \code{3})
#' @param write_out  Logical. If \code{TRUE}, the result table is written to
#'                   a \code{.csv} file (default \code{FALSE})
#' @param file_name  Character scalar. File path and name (without extension)
#'                   for the CSV output when \code{write_out = TRUE}
#'                   (default \code{"MRBMA"})
#'
#' @details
#' The Marginal Posterior Probability of Inclusion (MPPI) for exposure \eqn{j}
#' is computed as \eqn{\text{MPPI}_j = \sum_{\gamma \ni j} p(\gamma \mid \text{data})},
#' where the sum runs over all visited models that include exposure \eqn{j}
#'
#' @export
#'
#' @return A \code{data.frame} with \code{top} rows and three columns:
#' \describe{
#'   \item{\code{Exposure}}{Name of the exposure}
#'   \item{\code{Marginal Posterior Probability of Inclusion (MPPI)}}{Marginal inclusion
#'         probability, rounded to \code{digits} decimal places}
#'   \item{\code{BMA Causal Estimate}}{BMA-averaged causal estimate, rounded
#'         to \code{digits} decimal places}
#' }
#' Rows are sorted by \code{Marginal Posterior Probability of Inclusion (MPPI)} in
#' descending order
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{report_best_models}}, \code{\link{MRBMA_pvalues}}
#'
#' @examples
#' \dontrun{
#'
#' data(AMD_data)
#' betaX_ivw   <- as.matrix(AMD_data$betaX) / AMD_data$seAMD
#' betaAMD_ivw <- AMD_data$betaAMD / AMD_data$seAMD
#' rs          <- AMD_data$annotate[, 1]
#' rf          <- colnames(AMD_data$betaX)
#'
#' mvmr_input <- new("mvMRInput",
#'                   betaX    = betaX_ivw,
#'                   betaY    = as.matrix(betaAMD_ivw),
#'                   snps     = rs,
#'                   exposure = rf,
#'                   outcome  = "AMD")
#'
#' bma_out <- MRBMA_SSS(mvmr_input, kmin = 1, kmax = 12, max_iter = 100)
#' report_MRBMA(bma_out, top = 10)
#'
#' }
report_MRBMA <- function(output,
                          top       = 10,
                          digits    = 3,
                          write_out = FALSE,
                          file_name = "MRBMA") {

  # ---- Validate class ----
  if (!is(output, "MRBF") && !is(output, "mvMR_SSS"))
    stop("'output' must be of class 'MRBF' or 'mvMR_SSS'.")

  pp_marginal <- output@pp_marginal
  bma         <- output@BMAve_Estimate
  rf          <- output@Exposure

  if (length(pp_marginal) == 0)
    stop("pp_marginal is empty in the output object.")

  # ---- Guard: top cannot exceed number of exposures ----
  n_top <- min(top, length(pp_marginal))
  if (n_top < top)
    warning(sprintf(
      "Only %d exposures available (requested top = %d).",
      n_top, top
    ))

  # ---- Sort by marginal inclusion ----
  sort_obj <- sort.int(pp_marginal, index.return = TRUE, decreasing = TRUE)
  idx      <- sort_obj$ix[seq_len(n_top)]

  out <- data.frame(
    "Exposure"                                = rf[idx],
    "Marginal Posterior Probability of Inclusion (MPPI)" = round(pp_marginal[idx], digits),
    "BMA Causal Estimate"                   = round(bma[idx], digits),
    stringsAsFactors                          = FALSE,
    check.names                               = FALSE
  )

  if (write_out) {
    fname <- if (grepl("\\.csv$", file_name, ignore.case = TRUE)) {
      file_name
    } else {
      paste0(file_name, ".csv")
    }
    write.csv(out, file = fname, row.names = FALSE)
    message("MRBMA table written to: ", fname)
  }

  return(out)
}
