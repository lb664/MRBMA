# 135 characters ######################################################################################################################
#' @title Detect influential and outlying genetic instruments for MRBMA
#' @description For each model whose posterior probability exceeds a
#' user-defined threshold, computes two diagnostic statistics: (i) Cook's distance to identify SNPs that exert disproportionate
#' influence on the causal estimate, and (ii) a Q-statistic to detect heterogeneous or pleiotropic SNPs whose association with the
#' outcome deviates from what is predicted by the causal model. Both diagnostics are evaluated using the Rcpp function
#' \code{cooksD_cpp} for the hat matrix computations.
#'
#' SNP names and prior parameters are always extracted from the output object slots.
#'
#' @param output        An object of class \code{MRBF} or \code{mvMR_SSS}, as
#'                      returned by \code{MRBMA_exhaustive} or \code{MRBMA_SSS}
#' @param diag_ppthresh Numeric scalar between 0 and 1. Only models with
#'                      posterior probability \eqn{\geq} \code{diag_ppthresh}
#'                      are examined (default \code{0.02}). Increase to focus
#'                      on the most probable models; decrease to be more
#'                      inclusive
#' @param top           Positive integer. Maximum number of top models (by
#'                      posterior probability) to consider as candidates before
#'                      applying the threshold filter (default \code{100})
#' @param digits        Non-negative integer. Decimal places for rounding in
#'                      the internal model summary (default \code{3})
#'
#' @details
#' Cook's distance for SNP \eqn{i} in model \eqn{\gamma} is computed as defined
#' in \insertCite{Cook1977;textual}{MRBMA}:
#' \deqn{D_i = \frac{e_i^2}{s^2 \cdot k} \cdot \frac{h_{ii}}{(1 - h_{ii})^2}}
#' where \eqn{e_i} is the residual, \eqn{h_{ii}} is the hat-matrix diagonal,
#' \eqn{k = |\gamma|} is the model size, and \eqn{s^2} is the residual
#' variance. The threshold is the 50th percentile of the
#' \eqn{F(k, n - k)} distribution. The Q-statistic threshold applies a
#' Bonferroni correction: \eqn{\chi^2_{1, \, 1 - 0.05/n}}.
#'
#' A SNP flagged by Cook's distance is considered influential; one flagged by
#' the Q-statistic is considered an outlier. The union of both sets
#' (\code{rm_any}) is recommended for sensitivity analyses
#'
#' @export
#'
#' @return A named list (or \code{NULL} if no model exceeds the threshold):
#' \describe{
#'   \item{\code{rmCD}}{Character vector of SNP names flagged by Cook's
#'         distance (influential SNPs)}
#'   \item{\code{rmO}}{Character vector of SNP names flagged by the
#'         Q-statistic (outlier SNPs)}
#'   \item{\code{rm_any}}{Character vector: sorted union of \code{rmCD} and
#'         \code{rmO}}
#'   \item{\code{rmCD_idx}}{Integer vector of row indices of \code{rmCD} SNPs
#'         in \code{betaX}}
#'   \item{\code{rmO_idx}}{Integer vector of row indices of \code{rmO} SNPs
#'         in \code{betaX}}
#'   \item{\code{cD_mat}}{Numeric matrix of Cook's distances
#'         (\code{n_SNPs x nr_models})}
#'   \item{\code{Q_mat}}{Numeric matrix of Q-statistics
#'         (\code{n_SNPs x nr_models})}
#'   \item{\code{models_examined}}{Character vector of model keys that exceeded
#'         \code{diag_ppthresh}}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{MRBMA_SSS}}, \code{\link{MRBMA_exhaustive}},
#' \code{\link{report_best_models}}
#'
#' @importFrom stats qchisq qf
#'
#' @examples
#' \dontrun{
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
#' bma_out  <- MRBMA_SSS(mvmr_input, kmin = 1, kmax = 12, max_iter = 100)
#' diag_out <- MRBMA_diagnostics(bma_out, diag_ppthresh = 0.02)
#'
#' # SNPs to consider removing:
#' diag_out$rm_any
#' }
MRBMA_diagnostics <- function(output,
                               diag_ppthresh = 0.02,
                               top           = 100,
                               digits        = 3) {

  # ---- Validate class ----
  if (!is(output, "MRBF") && !is(output, "mvMR_SSS"))
    stop("'output' must be of class 'MRBF' or 'mvMR_SSS'.")

  # ---- Extract objects ----
  pp     <- output@pp
  models <- output@tupel
  betaX  <- output@betaX
  betaY  <- as.vector(output@betaY)
  rf     <- output@Exposure
  sigma  <- output@sigma

  n_SNPs <- nrow(betaX)

  # ---- SNP names: use rownames if available, else generic labels ----
  snp_names <- rownames(betaX)
  if (is.null(snp_names) || length(snp_names) == 0)
    snp_names <- paste0("SNP_", seq_len(n_SNPs))

  # ---- Parameter guards ----
  if (diag_ppthresh < 0 || diag_ppthresh > 1)
    stop("diag_ppthresh must be between 0 and 1.")
  if (top < 1)
    stop("top must be >= 1.")

  # ---- Rank models by posterior probability ----
  n_top    <- min(top, length(pp))
  sort_obj <- sort.int(pp, index.return = TRUE, decreasing = TRUE)
  top_keys <- models[sort_obj$ix][seq_len(n_top)]
  top_pp   <- pp[sort_obj$ix][seq_len(n_top)]

  # ---- Filter by threshold ----
  above   <- which(top_pp >= diag_ppthresh)
  nr_diag <- length(above)

  if (nr_diag == 0) {
    warning(sprintf(
      "No models found with posterior probability >= %.3f. Try lowering diag_ppthresh.",
      diag_ppthresh
    ))
    return(NULL)
  }

  models_diag <- top_keys[above]
  message(sprintf(
    "[MRBMA_diagnostics] Examining %d model(s) with pp >= %.3f.",
    nr_diag, diag_ppthresh
  ))

  # ---- Cook's distance and Q-statistic for each qualifying model ----
  cD_mat        <- matrix(NA_real_, nrow = n_SNPs, ncol = nr_diag)
  cD_thresh_vec <- numeric(nr_diag)
  Q_mat         <- matrix(NA_real_, nrow = n_SNPs, ncol = nr_diag)

  for (j in seq_len(nr_diag)) {
    tupel_j   <- as.integer(unlist(strsplit(models_diag[j], ",")))
    betaX_mod <- as.matrix(betaX[, tupel_j, drop = FALSE])
    sigma_vec <- rep(sigma, ncol(betaX_mod))

    # ---- Cook's distance via Rcpp ----
    cd_res           <- cooksD_cpp(y = betaY, x = betaX_mod, sigma_vec = sigma_vec)
    cD_mat[, j]      <- cd_res$cooksD
    cD_thresh_vec[j] <- cd_res$cooksD_thresh

    # ---- Q-statistic (hat matrix) ----
    sigma_diag_inv <- diag(sigma_vec^{-2}, nrow = ncol(betaX_mod))
    H_fm           <- betaX_mod %*%
                      solve(t(betaX_mod) %*% betaX_mod + sigma_diag_inv) %*%
                      t(betaX_mod)
    pred_y         <- as.vector(H_fm %*% betaY)
    Q_mat[, j]     <- (betaY - pred_y)^2
  }

  # ---- Flag influential SNPs (Cook's distance) ----
  maxCD            <- apply(cD_mat, 1, max)
  cD_global_thresh <- max(cD_thresh_vec)
  rmCD_idx         <- which(maxCD > cD_global_thresh)
  rmCD             <- snp_names[rmCD_idx]

  # ---- Flag outlier SNPs (Q-statistic, Bonferroni corrected) ----
  p_adj_thresh <- 0.05 / n_SNPs
  q_thresh     <- qchisq(p_adj_thresh, df = 1, lower.tail = FALSE)
  maxQ         <- apply(Q_mat, 1, max)
  rmO_idx      <- which(maxQ > q_thresh)
  rmO          <- snp_names[rmO_idx]

  rm_any <- sort(unique(c(rmCD, rmO)))

  message(sprintf(
    "[MRBMA_diagnostics] Influential SNPs (Cook's D): %d | Outlier SNPs (Q-stat): %d | Union: %d",
    length(rmCD), length(rmO), length(rm_any)
  ))

  list(
    rmCD            = rmCD,
    rmO             = rmO,
    rm_any          = rm_any,
    rmCD_idx        = rmCD_idx,
    rmO_idx         = rmO_idx,
    cD_mat          = cD_mat,
    Q_mat           = Q_mat,
    models_examined = models_diag
  )
}
