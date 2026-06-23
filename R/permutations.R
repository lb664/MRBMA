# 135 characters ######################################################################################################################
#' @title Run the permutation procedure for empirical p-values (MRBMA)
#' @description Implements the permutation-based null distribution procedure
#' of \insertCite{Levin2021;textual}{MRBMA} to compute empirical p-values for the Marginal Posterior Probabilities of Inclusion (MPPIs)
#' obtained from MRBMA. The observed SNP-outcome association vector \code{betaY} is permuted \code{nrepeat} times and
#' \code{MRBMA_exhaustive} or \code{MRBMA_SSS} is re-run for each permutation replicate with the same hyperparameters stored in the
#' \code{output} object, producing a null MPPI matrix that is passed to \code{MRBMA_pvalues}
#'
#' @param output      An object of class \code{MRBF} or \code{mvMR_SSS}. The
#'                    method (exhaustive or SSS) and all hyperparameters
#'                    (\code{sigma}, \code{prior_prob} and, for SSS,
#'                    \code{kmin}, \code{kmax}, \code{max_iter}) are read
#'                    directly from this object
#' @param nrepeat     Positive integer. Number of permutation replicates
#'                    (default \code{100000}). Use \code{100} to assess
#'                    runtime; use at least \code{10000} for preliminary
#'                    results; \code{100000} is recommended for publication
#'                    as in \insertCite{Levin2021;textual}{MRBMA}
#' @param save_matrix Logical. If \code{TRUE}, the permutation matrix is saved
#'                    to disk as an \code{.RData} file (default \code{TRUE}).
#'                    Strongly recommended given the long runtime; the matrix
#'                    can be loaded and passed to \code{MRBMA_pvalues} in a
#'                    later session
#' @param file_name   Character scalar. File path and base name for the saved
#'                    \code{.RData} file when \code{save_matrix = TRUE}
#'                    (default \code{"permutation_MRBMA"}). The \code{.RData}
#'                    extension is appended automatically if absent
#' @param verbose     Logical. If \code{TRUE}, prints progress to the console
#'                    at every 10\% of permutation replicates (default
#'                    \code{TRUE})
#'
#' @details
#' The function automatically detects whether the input object is of class
#' \code{MRBF} (exhaustive search) or \code{mvMR_SSS} (stochastic search) and
#' calls the appropriate function for each permutation replicate. Failed
#' permutations (e.g., due to near-singular summary-level statistic matrices in
#' edge cases) are skipped with a warning and their row in the output matrix
#' remains zero.
#'
#' Throughout, \eqn{d} denotes the number of exposures (columns of
#' \code{betaX}), \eqn{n} the number of SNPs (rows of \code{betaX}), and
#' \eqn{k = |\gamma|} the model size.
#'
#' Runtime scales linearly with \code{nrepeat} and is approximately
#' \eqn{n_{\text{repeat}} \times t_{\text{single}}}, where
#' \eqn{t_{\text{single}}} is the runtime of a single MRBMA run. For large
#' analyses, consider running this function on a remote server and setting
#' \code{save_matrix = TRUE}
#'
#' @export
#'
#' @return A numeric matrix of dimension \code{nrepeat x d_Exposure}, where each
#' row contains the Marginal Posterior Probabilities of Inclusion (MPPIs) from one permutation
#' replicate. This matrix is the direct input for \code{MRBMA_pvalues}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{MRBMA_pvalues}}, \code{\link{MRBMA_SSS}},
#' \code{\link{MRBMA_exhaustive}}
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
#' # Run 100 permutations to test runtime (use nrepeat >= 10000 for real analyses)
#' perm_mat <- MRBMA_permutations(bma_out, nrepeat = 100,
#'                                 save_matrix = FALSE, verbose = TRUE)
#'
#' }
MRBMA_permutations <- function(output,
                                nrepeat     = 100000,
                                save_matrix = TRUE,
                                file_name   = "permutation_MRBMA",
                                verbose     = TRUE) {

  # ---- Validate class ----
  if (!is(output, "MRBF") && !is(output, "mvMR_SSS"))
    stop("'output' must be of class 'MRBF' or 'mvMR_SSS'.")

  # ---- Extract shared information ----
  betaX      <- output@betaX
  betaY      <- output@betaY
  sigma      <- output@sigma
  prior_prob <- output@prior_prob
  rf         <- output@Exposure

  n_SNPs     <- nrow(betaX)
  d_Exposure <- ncol(betaX)

  snp_names <- rownames(betaX)
  if (is.null(snp_names)) snp_names <- as.character(seq_len(n_SNPs))

  # ---- Detect method and extract SSS parameters if needed ----
  use_SSS <- is(output, "mvMR_SSS")
  if (use_SSS) {
    kmin     <- output@kmin
    kmax     <- output@kmax
    max_iter <- output@max_iter
  }

  # ---- Guards ----
  if (!is.numeric(nrepeat) || nrepeat < 1)
    stop("nrepeat must be a positive integer.")
  nrepeat <- as.integer(nrepeat)

  if (nrepeat < 100)
    warning(sprintf(
      "nrepeat = %d is very small. For stable p-values, use at least 10,000.",
      nrepeat
    ))

  message(sprintf(
    "[MRBMA_permutations] Running %d permutations (%s method). This may take a while.",
    nrepeat,
    if (use_SSS) "SSS" else "exhaustive"
  ))

  # ---- Permutation loop ----
  permute_bma <- matrix(0.0, nrow = nrepeat, ncol = d_Exposure)
  progress_at <- unique(round(seq(0, nrepeat, length.out = 11)))[-1]

  for (i in seq_len(nrepeat)) {

    if (verbose && i %in% progress_at)
      message(sprintf("  [MRBMA_permutations] %d / %d (%.0f%%)",
                      i, nrepeat, 100 * i / nrepeat))

    perm_idx   <- sample.int(n_SNPs)
    betaY_perm <- as.matrix(betaY[perm_idx, , drop = FALSE])

    perm_input <- new("mvMRInput",
      betaX    = betaX,
      betaY    = betaY_perm,
      snps     = snp_names,
      exposure = rf,
      outcome  = "permutation"
    )

    perm_out <- tryCatch({
      if (use_SSS) {
        MRBMA_SSS(
          perm_input,
          kmin       = kmin,
          kmax       = kmax,
          max_iter   = max_iter,
          sigma      = sigma,
          prior_prob = prior_prob,
          verbose    = FALSE
        )
      } else {
        MRBMA_exhaustive(
          perm_input,
          sigma      = sigma,
          prior_prob = prior_prob
        )
      }
    }, error = function(e) {
      warning(sprintf("Permutation %d failed: %s -- skipping.", i, conditionMessage(e)))
      return(NULL)
    })

    if (is.null(perm_out)) next

    permute_bma[i, ] <- perm_out@pp_marginal
  }

  # ---- Save matrix ----
  if (save_matrix) {
    fname <- if (grepl("\\.(RData|rda)$", file_name, ignore.case = TRUE)) {
      file_name
    } else {
      paste0(file_name, ".RData")
    }
    save(permute_bma, file = fname)
    message("[MRBMA_permutations] Matrix saved to: ", fname)
  }

  message("[MRBMA_permutations] Done.")
  return(permute_bma)
}


# 135 characters ######################################################################################################################
#' @title Compute empirical p-values from the permutation matrix
#' @description Calculates empirical p-values for each exposure by
#' comparing the observed Marginal Posterior Probability of Inclusion (MPPI) against the null distribution generated by
#' \code{MRBMA_permutations}, then applies Benjamini-Hochberg false discovery rate (FDR) correction to adjust for testing across all
#' exposures simultaneously, following the procedure of \insertCite{Levin2021;textual}{MRBMA}
#'
#' @param output      An object of class \code{MRBF} or \code{mvMR_SSS}
#'                    containing the observed MPPIs from the original (unpermuted)
#'                    MRBMA run
#' @param permute_bma Numeric matrix of dimension \code{nrepeat x d_Exposure}
#'                    containing the permuted MPPIs, as returned by
#'                    \code{MRBMA_permutations}. The number of columns must
#'                    equal the number of exposures in \code{output}
#'
#' @details
#' Throughout, \eqn{d} denotes the number of exposures (columns of
#' \code{betaX}).
#'
#' The empirical p-value for exposure \eqn{j} is computed with add-one
#' smoothing to avoid zero p-values:
#' \deqn{p_j = \frac{\sum_{r=1}^{R} \mathbf{1}[\tilde{\pi}_j^{(r)} > \hat{\pi}_j] + 1}{R + 1}}
#' where \eqn{\hat{\pi}_j} is the observed MPPI, \eqn{\tilde{\pi}_j^{(r)}} is
#' the permuted MPPI at replicate \eqn{r}, and \eqn{R} is \code{nrepeat}.
#' Benjamini-Hochberg FDR correction \insertCite{BenjaminiHochberg1995;textual}{MRBMA} is then applied across all \eqn{d} exposures
#' using \code{p.adjust(..., method = "BH")}
#'
#' @export
#'
#' @return A \code{data.frame} with \code{d_Exposure} rows sorted by
#' \code{P-value} (ascending) and four columns:
#' \describe{
#'   \item{\code{Exposure}}{Name of the exposure}
#'   \item{\code{Marginal Posterior Probability of Inclusion (MPPI)}}{Observed marginal
#'         inclusion probability from the original MRBMA run}
#'   \item{\code{P-value}}{Empirical p-value (add-one smoothed)}
#'   \item{\code{FDR}}{Benjamini-Hochberg FDR-adjusted p-value}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{MRBMA_permutations}}, \code{\link{report_MRBMA}}
#'
#' @importFrom stats p.adjust
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
#' perm_mat <- MRBMA_permutations(bma_out, nrepeat = 100, save_matrix = FALSE)
#' MRBMA_pvalues(bma_out, perm_mat)
#' }
MRBMA_pvalues <- function(output, permute_bma) {

  # ---- Validate class ----
  if (!is(output, "MRBF") && !is(output, "mvMR_SSS"))
    stop("'output' must be of class 'MRBF' or 'mvMR_SSS'.")
  if (!is.matrix(permute_bma))
    stop("'permute_bma' must be a matrix (returned by MRBMA_permutations()).")

  mip_obs <- output@pp_marginal
  rf      <- output@Exposure

  # ---- Dimension check ----
  if (ncol(permute_bma) != length(mip_obs))
    stop(sprintf(
      "permute_bma has %d columns but output has %d exposures.",
      ncol(permute_bma), length(mip_obs)
    ))

  n_perm <- nrow(permute_bma)
  if (n_perm < 100)
    warning(sprintf(
      "permute_bma has only %d rows. P-value resolution is 1/%d. Use more permutations.",
      n_perm, n_perm + 1
    ))

  # ---- Empirical p-values (add-one smoothing) ----
  p_val <- vapply(seq_along(mip_obs), function(i) {
    (sum(permute_bma[, i] > mip_obs[i]) + 1) / (n_perm + 1)
  }, numeric(1))

  # ---- BH-FDR correction ----
  p_adj <- p.adjust(p_val, method = "BH")

  res <- data.frame(
    "Exposure"                                = rf,
    "Marginal Posterior Probability of Inclusion (MPPI)" = mip_obs,
    "P-value"                                 = p_val,
    "FDR"                                     = p_adj,
    stringsAsFactors                          = FALSE,
    check.names                               = FALSE
  )

  res[order(res[["P-value"]]), ]
}
