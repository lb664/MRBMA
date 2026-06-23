########## .pp_from_log_evidence ##########
#
# Convert a vector of log10-scale evidence values (logBF + logprior) to a
# vector of normalised posterior probabilities, with numerical rescaling to
# prevent 10^x overflow when the maximum log10-evidence exceeds 308.
#
# Arguments:
#   log_evidence  numeric vector of log10(BF) + log10(prior) values
#
# Returns:
#   numeric vector of posterior probabilities summing to 1
#
.pp_from_log_evidence <- function(log_evidence) {
  max_ev <- max(log_evidence)
  if (max_ev > 308) {
    log_evidence <- log_evidence - (max_ev - 308 + 1)
    max_ev       <- 308
  }
  ev        <- 10^(log_evidence - max_ev)
  sum_calib <- sum(ev) * 10^max_ev
  pp        <- 10^log_evidence / sum_calib
  pp
}

########## .marginal_inclusion ##########
#
# Compute the Marginal Posterior Probability of Inclusion (MPPI) for each exposure by
# summing the posterior probabilities of all models that include it.
#
# Arguments:
#   keys   character vector of comma-separated column indices into betaX
#   pp     numeric vector of posterior probabilities (same length as keys)
#   d_Exposure   integer: total number of exposures
#
# Returns:
#   numeric vector of length d_Exposure containing the MPPIs
#
.marginal_inclusion <- function(keys, pp, d_Exposure) {
  index_mat <- matrix(0L, nrow = d_Exposure, ncol = length(pp))
  parsed    <- strsplit(keys, ",")
  for (i in seq_along(pp)) {
    idx               <- as.integer(parsed[[i]])
    index_mat[idx, i] <- 1L
  }
  pp_mat <- matrix(pp, nrow = d_Exposure, ncol = length(pp), byrow = TRUE)
  rowSums(pp_mat * index_mat)
}

# 135 characters ######################################################################################################################
#' @title Run MRBMA with exhaustive evaluation of all models
#' @description Evaluates all \eqn{2^d - 1} possible sub-models, where
#' \eqn{d} is the number of exposures, and returns the full posterior distribution over models together with model-averaged causal
#' estimates and Marginal Posterior Probabilities of Inclusion (MPPIs). The computation delegates to the Rcpp function
#' \code{exhaustive_bf_cpp}, which pre-computes summary-level statistics \eqn{X^\top X}, \eqn{X^\top Y} and \eqn{Y^\top Y} once and
#' evaluates the log10 Bayes Factor (BF) for every combination via Armadillo matrix routines. Recommended for \eqn{d \leq 12}
#' exposures; for larger \eqn{d} use \code{MRBMA_SSS}.
#'
#' @param object     An object of class \code{mvMRInput} containing the
#'                   SNP-exposure and SNP-outcome summary associations
#' @param sigma      Positive numeric scalar. Prior standard deviation for
#'                   causal effects, applied identically to all exposures
#'                   (default \code{0.5}). Corresponds to \eqn{\sigma} in
#'                   \insertCite{Zuber2020;textual}{MRBMA}
#' @param prior_prob Numeric scalar strictly between 0 and 1. Prior probability
#'                   of inclusion for each exposure independently
#'                   (default \code{0.5}). Smaller values favour sparser models
#'
#' @details
#' Throughout, \eqn{d} denotes the number of exposures (columns of
#' \code{betaX}), \eqn{n} the number of SNPs (rows of \code{betaX}), and
#' \eqn{k = |\gamma|} the model size (\eqn{1 \leq k \leq d}).
#'
#' The Bayes Factor (BF) for a model \eqn{\gamma} is computed from the closed-form
#' expression derived in \insertCite{Zuber2020;textual}{MRBMA} using summary-level
#' statistics only, which makes the exhaustive search feasible in pure
#' summary-level data. The prior on model size is a product Bernoulli with
#' parameter \code{prior_prob}. A log-sum-exp rescaling is applied to prevent
#' numerical overflow when normalising posterior probabilities.
#'
#' Data are always extracted from the input object slots.
#'
#' @export
#'
#' @return An object of class \code{MRBF} with the following slots:
#' \describe{
#'   \item{\code{Exposure}}{Character vector of exposure names}
#'   \item{\code{Outcome}}{Character scalar naming the outcome}
#'   \item{\code{BMAve_Estimate}}{Numeric vector of BMA-averaged causal
#'         estimates, one per exposure}
#'   \item{\code{BestModel_Estimate}}{Numeric vector of causal estimates from
#'         the highest-posterior model}
#'   \item{\code{BestModel}}{Comma-separated exposure indices (column positions
#'         in betaX) of exposures in the best model}
#'   \item{\code{tupel}}{Character vector of all \eqn{2^d - 1} model keys}
#'   \item{\code{pp}}{Numeric vector of posterior probabilities}
#'   \item{\code{pp_marginal}}{Numeric vector of marginal inclusion
#'         probabilities, one per exposure}
#'   \item{\code{betaX}}{SNP-exposure association matrix (stored for
#'         reporting)}
#'   \item{\code{betaY}}{SNP-outcome association matrix (stored for
#'         reporting)}
#'   \item{\code{sigma}}{Prior standard deviation used}
#'   \item{\code{prior_prob}}{Prior inclusion probability used}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{MRBMA_SSS}}, \code{\link{report_best_models}},
#' \code{\link{report_MRBMA}}, \code{\link{MRBMA_diagnostics}}
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
#' # Exhaustive search (feasible here because d_Exposure is small)
#' bma_out <- MRBMA_exhaustive(mvmr_input, sigma = 0.5, prior_prob = 0.5)
#'
#' # Report top 5 models and exposures
#' report_best_models(bma_out, top = 5)
#' report_MRBMA(bma_out, top = 5)
#'
#' }
MRBMA_exhaustive <- function(object, sigma = 0.5, prior_prob = 0.5) {

  # ---- Input validation ----
  if (!is(object, "mvMRInput"))
    stop("'object' must be of class 'mvMRInput'.")

  bX <- object@betaX
  bY <- object@betaY

  if (!is.matrix(bX) || !is.matrix(bY))
    stop("betaX and betaY must be matrices.")
  if (nrow(bX) != nrow(bY))
    stop(sprintf(
      "betaX (%d rows) and betaY (%d rows) must have the same number of SNPs.",
      nrow(bX), nrow(bY)
    ))
  if (ncol(bY) != 1)
    stop("betaY must be a single-column matrix (one outcome).")
  if (any(is.na(bX)) || any(is.na(bY)))
    stop("betaX and betaY must not contain NA values.")
  if (sigma <= 0)
    stop("sigma must be a positive number.")
  if (prior_prob <= 0 || prior_prob >= 1)
    stop("prior_prob must be strictly between 0 and 1.")

  n_SNPs     <- nrow(bX)
  d_Exposure <- ncol(bX)

  if (d_Exposure > 20)
    warning(sprintf(
      paste0("d_Exposure = %d. Exhaustive search over 2^%d - 1 = %.0f models may be ",
             "very slow. Consider MRBMA_SSS() instead."),
      d_Exposure, d_Exposure, 2^d_Exposure - 1
    ))

  sigma_vec <- rep(sigma, d_Exposure)

  # ---- Sufficient statistics (computed once) ----
  XtX <- t(bX) %*% bX   # d x d
  XtY <- t(bX) %*% bY   # d x 1
  YtY <- as.numeric(t(bY) %*% bY)

  if (YtY <= 0)
    stop("YtY (t(betaY) %*% betaY) must be positive. Check betaY values.")

  # ---- Call Rcpp exhaustive search ----
  res <- exhaustive_bf_cpp(
    XtY        = as.vector(XtY),
    XtX        = XtX,
    YtY        = YtY,
    sigma_vec  = sigma_vec,
    d_Exposure = d_Exposure,
    n_SNPs     = n_SNPs,
    prior_prob = prior_prob
  )

  keys      <- res$keys
  logBF     <- as.numeric(res$logBF)
  logprior  <- as.numeric(res$logprior)
  theta_mat <- res$theta_mat

  if (length(keys) == 0)
    stop("Rcpp exhaustive search returned no models. Check your input data.")

  # ---- Posterior probabilities ----
  pp <- .pp_from_log_evidence(logBF + logprior)

  if (any(!is.finite(pp)))
    stop("Non-finite posterior probabilities detected. Check for collinear exposures.")

  # ---- BMA-averaged causal estimates ----
  pp_mat       <- matrix(pp, nrow = d_Exposure, ncol = length(pp), byrow = TRUE)
  BMA_estimate <- rowSums(pp_mat * theta_mat)

  # ---- Marginal inclusion probabilities ----
  pp_marginal  <- .marginal_inclusion(keys, pp, d_Exposure)

  # ---- Best model ----
  best_idx       <- which.max(logBF + logprior)
  best_model     <- keys[best_idx]
  best_model_est <- theta_mat[, best_idx]

  # ---- Assemble output ----
  new("MRBF",
    Exposure           = object@exposure,
    Outcome            = object@outcome,
    BMAve_Estimate     = BMA_estimate,
    BestModel_Estimate = as.numeric(best_model_est),
    BestModel          = best_model,
    tupel              = keys,
    pp                 = pp,
    pp_marginal        = pp_marginal,
    betaX              = bX,
    betaY              = bY,
    sigma              = sigma,
    prior_prob         = prior_prob
  )
}
