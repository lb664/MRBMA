# 135 characters ######################################################################################################################
#' @title Run MRBMA with stochastic search (SSS)
#' @description Implements the shotgun stochastic search (SSS) algorithm of
#' \insertCite{Zuber2020;textual}{MRBMA} to explore the model space when an exhaustive enumeration of all \eqn{2^d - 1} models is
#' computationally infeasible (typically \eqn{d > 12}). The core search loop and memoisation hash table are implemented in C++ via
#' \code{sss_core_cpp}, replacing the R \code{hash} package for O(1)-average memoisation.
#'
#' When \code{k_min == k_max}, an exhaustive deterministic search of all models of size \eqn{1, \ldots, k_{\min}} is performed and the
#' SSS loop is skipped. When \code{kmin < kmax}, the algorithm alternates between neighbourhood expansion (add, delete, swap moves) and
#' a weighted random draw proportional to posterior evidence
#'
#' @param object     An object of class \code{mvMRInput} containing the
#'                   SNP-exposure and SNP-outcome summary associations
#' @param kmin       Positive integer \eqn{k_{\min}}. Minimum model size
#'                   explored (default \code{1}). All models of size
#'                   \eqn{1, \ldots, k_{\min}} are evaluated exhaustively in
#'                   the deterministic phase
#' @param kmax       Positive integer \eqn{k_{\max}}. Maximum model size
#'                   explored (default \code{20}). Cannot exceed the number of
#'                   exposures \eqn{d}. Note that computing all combinations of
#'                   more than 12 exposures exhaustively is infeasible; SSS is
#'                   therefore recommended when \code{kmax} is large.
#' @param max_iter   Positive integer. Number of SSS iterations (default
#'                   \code{1000}). Ignored when \code{k_min == k_max}. Start
#'                   with a small value (e.g., \code{100}) to assess runtime;
#'                   use at least \code{10000}, ideally \code{100000}, for
#'                   final results
#' @param sigma      Positive numeric scalar. Prior standard deviation for
#'                   causal effects (default \code{0.5})
#' @param prior_prob Numeric scalar strictly between 0 and 1. Prior inclusion
#'                   probability per exposure (default \code{0.5}). For
#'                   high-dimensional settings with many exposures, a
#'                   smaller value (e.g., \code{0.1}) is recommended to favour
#'                   sparse models
#' @param verbose    Logical. If \code{TRUE}, prints iteration progress to the
#'                   console every 1,000 iterations (default \code{FALSE})
#'
#' @details
#' Throughout, \eqn{d} denotes the number of exposures (columns of
#' \code{betaX}), \eqn{n} the number of SNPs (rows of \code{betaX}), and
#' \eqn{k = |\gamma|} the current model size (\eqn{k_{\min} \leq k \leq k_{\max}}).
#'
#' The stochastic search maintains a \code{std::unordered_map} hash table in
#' C++ to memoise Bayes Factors (BFs) and causal estimates for all visited models,
#' avoiding redundant computations when the same model is reached from multiple
#' paths. Neighbourhood moves are: (i) deletion of one exposure (model size
#' \eqn{- 1}), (ii) swap of one exposure (model size unchanged), and (iii)
#' addition of one exposure (model size \eqn{+ 1}). The current model is
#' also retained in the neighbourhood to allow zero-move transitions.
#'
#' For details on the SSS algorithm see \insertCite{Hans2007;textual}{MRBMA}. For
#' the Bayes Factor (BF) derivation and the MR-BMA model see
#' \insertCite{Zuber2020;textual}{MRBMA}. For the permutation procedure to
#' compute empirical p-values see \insertCite{Levin2021;textual}{MRBMA}
#'
#' @export
#'
#' @return An object of class \code{mvMR_SSS} with the following slots:
#' \describe{
#'   \item{\code{Exposure}}{Character vector of exposure names}
#'   \item{\code{Outcome}}{Character scalar naming the outcome}
#'   \item{\code{BMAve_Estimate}}{Numeric vector of BMA-averaged causal
#'         estimates, one per exposure}
#'   \item{\code{BestModel_Estimate}}{Numeric vector of causal estimates from
#'         the highest-posterior model}
#'   \item{\code{BestModel}}{Comma-separated exposure indices (column positions
#'         in betaX) of exposures in the best model}
#'   \item{\code{tupel}}{Character vector of all visited model keys}
#'   \item{\code{pp}}{Numeric vector of posterior probabilities for all
#'         visited models}
#'   \item{\code{pp_marginal}}{Numeric vector of marginal inclusion
#'         probabilities, one per exposure}
#'   \item{\code{betaX}}{SNP-exposure association matrix (stored for
#'         downstream use)}
#'   \item{\code{betaY}}{SNP-outcome association matrix (stored for
#'         downstream use)}
#'   \item{\code{kmin}}{Minimum model size used in the search}
#'   \item{\code{kmax}}{Maximum model size used in the search}
#'   \item{\code{max_iter}}{Number of SSS iterations performed}
#'   \item{\code{sigma}}{Prior standard deviation used}
#'   \item{\code{prior_prob}}{Prior inclusion probability used}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{MRBMA_exhaustive}}, \code{\link{report_best_models}},
#' \code{\link{report_MRBMA}}, \code{\link{MRBMA_diagnostics}},
#' \code{\link{MRBMA_permutations}}
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
#' # Stochastic search with 100 iterations (increase max_iter for final results)
#' bma_out <- MRBMA_SSS(mvmr_input, kmin = 1, kmax = 12,
#'                       max_iter = 100, sigma = 0.5, prior_prob = 0.1)
#'
#' # Exhaustive deterministic search (k_min == k_max)
#' bma_out_det <- MRBMA_SSS(mvmr_input, kmin = 5, kmax = 5,
#'                           sigma = 0.5, prior_prob = 0.1)
#'
#' }
MRBMA_SSS <- function(object,
                       kmin       = 1,
                       kmax       = 20,
                       max_iter   = 1000,
                       sigma      = 0.5,
                       prior_prob = 0.5,
                       verbose    = FALSE) {

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
  if (!is.numeric(kmin) || kmin < 1)
    stop("kmin must be a positive integer >= 1.")
  if (!is.numeric(kmax) || kmax < 1)
    stop("kmax must be a positive integer >= 1.")
  if (kmin > kmax)
    stop(sprintf("kmin (%d) must be <= kmax (%d).", kmin, kmax))
  if (kmax > ncol(bX))
    stop(sprintf(
      "kmax (%d) cannot exceed the number of exposures (%d).",
      kmax, ncol(bX)
    ))
  if (!is.numeric(max_iter) || max_iter < 1)
    stop("max_iter must be a positive integer.")
  if (sigma <= 0)
    stop("sigma must be a positive number.")
  if (prior_prob <= 0 || prior_prob >= 1)
    stop("prior_prob must be strictly between 0 and 1.")

  n_SNPs     <- nrow(bX)
  d_Exposure <- ncol(bX)
  kmin       <- as.integer(kmin)
  kmax       <- as.integer(min(kmax, d_Exposure))
  max_iter   <- as.integer(max_iter)

  # ---- Informative messages ----
  if (kmin == kmax) {
    message(sprintf(
      paste0("[MRBMA_SSS] k_min == k_max == %d: running deterministic exhaustive ",
             "search of all models of size 1..%d (SSS loop skipped)."),
      kmin, kmin
    ))
  } else {
    message(sprintf(
      "[MRBMA_SSS] SSS search: k_min=%d, k_max=%d, max_iter=%d, d_Exposure=%d.",
      kmin, kmax, max_iter, d_Exposure
    ))
    if (max_iter < 10000)
      message("  [Note] For stable final results, recommend max_iter >= 10,000.")
  }

  sigma_vec <- rep(sigma, d_Exposure)

  # ---- Sufficient statistics (computed once) ----
  XtX <- t(bX) %*% bX
  XtY <- t(bX) %*% bY
  YtY <- as.numeric(t(bY) %*% bY)

  if (YtY <= 0)
    stop("YtY (t(betaY) %*% betaY) must be positive. Check betaY values.")
  if (any(diag(XtX) == 0))
    stop("One or more exposures have zero variance (XtX diagonal == 0).")

  # ---- Call Rcpp SSS ----
  res <- sss_core_cpp(
    XtY        = as.vector(XtY),
    XtX        = XtX,
    YtY        = YtY,
    sigma_vec  = sigma_vec,
    d_Exposure = d_Exposure,
    n_SNPs     = n_SNPs,
    prior_prob = prior_prob,
    kmin       = kmin,
    kmax       = kmax,
    max_iter   = max_iter,
    verbose    = verbose
  )

  keys      <- res$keys
  logBF     <- as.numeric(res$logBF)
  logprior  <- as.numeric(res$logprior)
  theta_mat <- res$theta_mat

  if (length(keys) == 0)
    stop("SSS returned no visited models. Check your input data.")

  # ---- Posterior probabilities ----
  pp <- .pp_from_log_evidence(logBF + logprior)

  if (any(!is.finite(pp)))
    stop("Non-finite posterior probabilities. Check for collinear exposures.")

  # ---- BMA-averaged causal estimates ----
  pp_mat       <- matrix(pp, nrow = d_Exposure, ncol = length(pp), byrow = TRUE)
  BMA_estimate <- rowSums(pp_mat * theta_mat)

  # ---- Marginal inclusion probabilities ----
  pp_marginal  <- .marginal_inclusion(keys, pp, d_Exposure)

  # ---- Best model ----
  best_idx       <- which.max(logBF + logprior)
  best_model     <- keys[best_idx]
  best_model_est <- theta_mat[, best_idx]

  message(sprintf(
    "[MRBMA_SSS] Done. Total models visited: %d. Best model: {%s}.",
    length(keys), best_model
  ))

  # ---- Assemble output ----
  new("mvMR_SSS",
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
    kmin               = as.numeric(kmin),
    kmax               = as.numeric(kmax),
    max_iter           = as.numeric(max_iter),
    sigma              = sigma,
    prior_prob         = prior_prob
  )
}
