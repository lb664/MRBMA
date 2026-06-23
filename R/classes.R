# 135 characters ######################################################################################################################
#' @title Input class for MRBMA analyses
#' @description S4 class that stores all summary-level data required to run a
#' multivariable Mendelian randomization Bayesian model averaging (MRBMA) analysis. It serves as the single entry point for both
#' \code{MRBMA_exhaustive} and \code{MRBMA_SSS}
#'
#' @slot betaX        Numeric matrix of SNP-exposure associations of dimension
#'                    \code{n_SNPs x d_Exposure}, typically obtained after
#'                    inverse-variance weighting (IVW)
#' @slot betaY        Numeric matrix of SNP-outcome associations of dimension
#'                    \code{n_SNPs x 1}. Must have the same number of rows as
#'                    \code{betaX}
#' @slot betaXse      Numeric matrix of standard errors for \code{betaX}
#'                    (optional; same dimension as \code{betaX}). Set to an
#'                    empty matrix if not available
#' @slot betaYse      Numeric matrix of standard errors for \code{betaY}
#'                    (optional; same dimension as \code{betaY}). Set to an
#'                    empty matrix if not available
#' @slot exposure     Character vector of exposure names, of length
#'                    \code{d_Exposure}. Must match \code{colnames(betaX)}
#' @slot outcome      Character scalar naming the outcome of interest (e.g.,
#'                    \code{"AMD"})
#' @slot snps         Character vector of SNP identifiers (e.g., rs-numbers),
#'                    of length \code{n_SNPs}. Must match \code{rownames(betaX)}
#' @slot effect_allele Character vector of effect alleles, one per SNP
#'                    (optional). Set to \code{character(0)} if not available
#' @slot other_allele Character vector of non-effect alleles, one per SNP
#'                    (optional). Set to \code{character(0)} if not available
#' @slot eaf          Numeric vector of effect allele frequencies, one per SNP
#'                    (optional). Set to \code{numeric(0)} if not available
#' @slot correlation  Numeric matrix of pairwise LD correlations between SNPs,
#'                    of dimension \code{n_SNPs x n_SNPs} (optional; supply an
#'                    identity matrix or empty matrix if SNPs are independent)
#'
#' @details
#' Throughout, \eqn{d} denotes the number of exposures (columns of
#' \code{betaX}) and \eqn{n} the number of SNPs (rows of \code{betaX}).
#'
#' Create an \code{mvMRInput} object using \code{new("mvMRInput", ...)}. Only
#' \code{betaX}, \code{betaY}, \code{snps}, \code{exposure} and \code{outcome}
#' are required for \code{MRBMA_exhaustive} and \code{MRBMA_SSS}. All other
#' slots are stored for completeness and downstream use (e.g., LD adjustment).
#' For details on the MRBMA model see \insertCite{Zuber2020;textual}{MRBMA}
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#'
#' # Minimal construction from pre-processed IVW summary data:
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
#' }
setClass(
  "mvMRInput",
  representation(
    betaX         = "matrix",
    betaY         = "matrix",
    betaXse       = "matrix",
    betaYse       = "matrix",
    exposure      = "character",
    outcome       = "character",
    snps          = "character",
    effect_allele = "character",
    other_allele  = "character",
    eaf           = "numeric",
    correlation   = "matrix"
  )
)

# 135 characters ######################################################################################################################
#' @title Output class for exhaustive Bayesian model averaging (MRBMA)
#' @description S4 class returned by \code{MRBMA_exhaustive}. It stores all
#' model-level and exposure-level summaries produced by an exhaustive evaluation of all \eqn{2^d - 1} sub-models, where \eqn{d} is the
#' number of exposures
#'
#' @slot Exposure           Character vector of exposure names (length
#'                          \code{d_Exposure})
#' @slot Outcome            Character scalar naming the outcome
#' @slot BMAve_Estimate     Numeric vector of BMA-averaged causal estimates,
#'                          one per exposure (length \code{d_Exposure})
#' @slot BestModel_Estimate Numeric vector of causal estimates from the model
#'                          with the highest posterior probability (length
#'                          \code{d_Exposure}; zero for exposures absent from the
#'                          best model)
#' @slot BestModel          Character scalar: comma-separated exposure indices
#'                          (column indices into betaX) of the best
#'                          model (e.g., \code{"2,5,7"})
#' @slot tupel              Character vector of all visited model keys, each
#'                          being comma-separated exposure indices (column
#'                          positions in betaX)
#' @slot pp                 Numeric vector of posterior probabilities, one per
#'                          model in \code{tupel} (sums to 1)
#' @slot pp_marginal        Numeric vector of Marginal Posterior Probabilities of Inclusion (MPPIs),
#'                          one per exposure (length \code{d_Exposure})
#' @slot betaX              SNP-exposure association matrix stored for
#'                          downstream reporting (same as input \code{betaX})
#' @slot betaY              SNP-outcome association matrix stored for
#'                          downstream reporting (same as input \code{betaY})
#' @slot sigma              Numeric scalar: prior standard deviation used in
#'                          the analysis
#' @slot prior_prob         Numeric scalar: prior inclusion probability used
#'                          in the analysis
#'
#' @details
#' The \code{sigma} and \code{prior_prob} slots are stored for permutation replay
#' in \code{MRBMA_permutations}. For the underlying BMA model see
#' \insertCite{Zuber2020;textual}{MRBMA}.
#'
#' @references
#' \insertAllCited{}
#'
setClass(
  "MRBF",
  representation(
    Exposure           = "character",
    Outcome            = "character",
    BMAve_Estimate     = "numeric",
    BestModel_Estimate = "numeric",
    BestModel          = "character",
    tupel              = "character",
    pp                 = "numeric",
    pp_marginal        = "numeric",
    betaX              = "matrix",
    betaY              = "matrix",
    sigma              = "numeric",
    prior_prob         = "numeric"
  )
)

# 135 characters ######################################################################################################################
#' @title Output class for stochastic search Bayesian model averaging (MRBMA)
#' @description S4 class returned by \code{MRBMA_SSS}. Contains the same
#' model-level and exposure-level summaries as \code{MRBF}, plus the search parameters needed to replay the stochastic search
#' identically in permutation runs
#'
#' @slot Exposure           Character vector of exposure names (length
#'                          \code{d_Exposure})
#' @slot Outcome            Character scalar naming the outcome
#' @slot BMAve_Estimate     Numeric vector of BMA-averaged causal estimates
#'                          (length \code{d_Exposure})
#' @slot BestModel_Estimate Numeric vector of causal estimates from the
#'                          highest-posterior model (length \code{d_Exposure})
#' @slot BestModel          Character scalar: comma-separated exposure indices
#'                          (column indices into betaX) of the best
#'                          model's exposures
#' @slot tupel              Character vector of all visited model keys
#' @slot pp                 Numeric vector of posterior probabilities for all
#'                          visited models
#' @slot pp_marginal        Numeric vector of Marginal Posterior Probabilities of Inclusion (MPPIs)
#'                          (length \code{d_Exposure})
#' @slot betaX              SNP-exposure association matrix (stored for
#'                          downstream use)
#' @slot betaY              SNP-outcome association matrix (stored for
#'                          downstream use)
#' @slot kmin               Numeric scalar: minimum model size explored during
#'                          the search
#' @slot kmax               Numeric scalar: maximum model size explored during
#'                          the search
#' @slot max_iter           Numeric scalar: number of SSS iterations performed
#' @slot sigma              Numeric scalar: prior standard deviation used
#' @slot prior_prob         Numeric scalar: prior inclusion probability used
#'
#' @details
#' Throughout, \eqn{d} denotes the number of exposures, \eqn{n} the number of
#' SNPs, \eqn{k = |\gamma|} the model size, and \eqn{k_{\min} \leq k \leq k_{\max}}
#' the search range.
#'
#' The \code{kmin}, \code{kmax}, \code{max_iter}, \code{sigma} and
#' \code{prior_prob} slots are stored so that \code{MRBMA_permutations} can
#' re-run \code{MRBMA_SSS} with identical settings for each permutation
#' replicate, without requiring the user to pass parameters twice. For model
#' details see \insertCite{Zuber2020;textual}{MRBMA}
#'
#' @references
#' \insertAllCited{}
#'
setClass(
  "mvMR_SSS",
  representation(
    Exposure           = "character",
    Outcome            = "character",
    BMAve_Estimate     = "numeric",
    BestModel_Estimate = "numeric",
    BestModel          = "character",
    tupel              = "character",
    pp                 = "numeric",
    pp_marginal        = "numeric",
    betaX              = "matrix",
    betaY              = "matrix",
    kmin               = "numeric",
    kmax               = "numeric",
    max_iter           = "numeric",
    sigma              = "numeric",
    prior_prob         = "numeric"
  )
)
