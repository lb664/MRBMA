#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// =============================================================================
// ---- Internal helpers (not exported) ----------------------------------------
// =============================================================================

// Build comma-separated string key from sorted column indices into betaX
static std::string make_key(const std::vector<int>& gamma_idx) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < gamma_idx.size(); ++i) {
    if (i > 0) oss << ',';
    oss << gamma_idx[i];
  }
  return oss.str();
}


// log10 Bayes factor — internal, works on 0-indexed gamma
static double logBF_internal(
    const arma::vec&  XtY,
    const arma::mat&  XtX,
    double            YtY,
    const arma::vec&  sigma_vec,
    const arma::uvec& gamma0,
    int               n
) {
  const arma::vec XtY_g   = XtY.elem(gamma0);
  const arma::mat XtX_g   = XtX.submat(gamma0, gamma0);
  const arma::vec sigma_g = sigma_vec.elem(gamma0);

  const arma::mat invnu    = arma::diagmat(arma::pow(sigma_g, -2.0));
  const arma::mat invOmega = invnu + XtX_g;

  arma::vec B;
  bool ok = arma::solve(B, invOmega, XtY_g, arma::solve_opts::fast);
  if (!ok) return R_NegInf;

  double log_det_nat, sign_val;
  arma::log_det(log_det_nat, sign_val, invOmega);
  if (sign_val <= 0.0) return R_NegInf;
  const double log10_det = log_det_nat / std::log(10.0);

  const double sum_log10_sigma = arma::sum(arma::log10(sigma_g));
  const double quad  = arma::as_scalar(B.t() * invOmega * B);
  const double denom = YtY - quad;
  if (denom <= 0.0) return R_NegInf;

  const double log10_ratio = std::log10(denom) - std::log10(YtY);
  return -0.5 * log10_det - sum_log10_sigma - (n / 2.0) * log10_ratio;
}

// Causal estimates — internal, works on 0-indexed gamma
static arma::vec beta_internal(
    const arma::vec&  XtY,
    const arma::mat&  XtX,
    const arma::vec&  sigma_vec,
    const arma::uvec& gamma0,
    int               d_Exposure
) {
  const arma::vec XtY_g   = XtY.elem(gamma0);
  const arma::mat XtX_g   = XtX.submat(gamma0, gamma0);
  const arma::vec sigma_g = sigma_vec.elem(gamma0);

  const arma::mat invnu    = arma::diagmat(arma::pow(sigma_g, -2.0));
  const arma::mat invOmega = invnu + XtX_g;

  arma::vec B;
  arma::solve(B, invOmega, XtY_g, arma::solve_opts::fast);

  arma::vec theta(d_Exposure, arma::fill::zeros);
  theta.elem(gamma0) = B;
  return theta;
}

// log10 binomial prior: log10(p^k * (1-p)^(d-k))
static double log_prior(int model_size, int d_Exposure, double prior_prob) {
  return model_size          * std::log10(prior_prob) +
         (d_Exposure - model_size) * std::log10(1.0 - prior_prob);
}

// Generate all combinations of {1,...,d_Exposure} of exactly size k
static std::vector<std::vector<int>>
gen_combos(int d_Exposure, int k) {
  std::vector<std::vector<int>> result;
  std::vector<int> combo;
  combo.reserve(k);

  std::function<void(int)> recurse = [&](int start) {
    if (static_cast<int>(combo.size()) == k) {
      result.push_back(combo);
      return;
    }
    int remaining = k - static_cast<int>(combo.size());
    for (int i = start; i <= d_Exposure - remaining + 1; ++i) {
      combo.push_back(i);
      recurse(i + 1);
      combo.pop_back();
    }
  };
  recurse(1);
  return result;
}

// Generate SSS neighbourhood: delete / swap / add / keep-current moves
static std::vector<std::vector<int>>
create_neighbourhood(const std::vector<int>& cur, int d_Exposure, int kmin, int kmax) {
  const int cur_size = static_cast<int>(cur.size());

  std::vector<int> possible;
  possible.reserve(d_Exposure - cur_size);
  for (int i = 1; i <= d_Exposure; ++i) {
    if (std::find(cur.begin(), cur.end(), i) == cur.end())
      possible.push_back(i);
  }

  std::vector<std::vector<int>> nbhd;

  // Delete moves
  if (cur_size > kmin) {
    for (int i = 0; i < cur_size; ++i) {
      std::vector<int> cfg = cur;
      cfg.erase(cfg.begin() + i);
      nbhd.push_back(cfg);
    }
  }

  // Swap moves
  for (int i = 0; i < cur_size; ++i) {
    std::vector<int> base = cur;
    base.erase(base.begin() + i);
    for (int k : possible) {
      std::vector<int> swapped = base;
      swapped.push_back(k);
      std::sort(swapped.begin(), swapped.end());
      nbhd.push_back(swapped);
    }
  }

  // Add moves
  if (cur_size < kmax) {
    for (int k : possible) {
      std::vector<int> cfg = cur;
      cfg.push_back(k);
      std::sort(cfg.begin(), cfg.end());
      nbhd.push_back(cfg);
    }
  }

  // Keep current
  nbhd.push_back(cur);

  // Deduplicate
  std::unordered_map<std::string, bool> seen;
  std::vector<std::vector<int>> unique_nbhd;
  unique_nbhd.reserve(nbhd.size());
  for (const auto& cfg : nbhd) {
    std::string key = make_key(cfg);
    if (seen.emplace(key, true).second)
      unique_nbhd.push_back(cfg);
  }
  return unique_nbhd;
}


// =============================================================================
// ---- Exported functions -----------------------------------------------------
// =============================================================================

///////// logBF_summary_cpp //////////
//' @title Compute log10 Bayes Factor for a single model
//' @description Evaluates the log10 Bayes Factor (BF) for a model defined by the
//' indicator vector \code{gamma_idx} using pre-computed summary-level statistics \eqn{X^\top X}, \eqn{X^\top Y} and \eqn{Y^\top Y}. This
//' avoids re-multiplying the full data matrices for every model evaluation and is the innermost computational kernel called by
//' \code{exhaustive_bf_cpp} and \code{sss_core_cpp}.
//' Throughout, \eqn{d} denotes the number of exposures, \eqn{n} the number of SNPs.
//'
//' @param XtY       Numeric vector of length \eqn{d}. Pre-computed
//'                  \eqn{X^\top Y = \texttt{betaX}^\top \texttt{betaY}}
//' @param XtX       Numeric matrix of dimension \eqn{d \times d}. Pre-computed
//'                  \eqn{X^\top X = \texttt{betaX}^\top \texttt{betaX}}
//' @param YtY       Positive numeric scalar. Pre-computed
//'                  \eqn{Y^\top Y = \texttt{betaY}^\top \texttt{betaY}}
//' @param sigma_vec Numeric vector of length \eqn{d}. Prior standard
//'                  deviations for causal effects, one per exposure
//' @param gamma_idx Integer vector of column indices into betaX identifying the exposures in model \eqn{\gamma}
//' @param n         Positive integer. Number of SNPs (instruments) \eqn{n}
//'
//' @return Scalar double: log10 Bayes Factor (BF) for model \eqn{\gamma}. Returns
//' \code{-Inf} if the system is (near-)singular
//'
//' @export
// [[Rcpp::export]]
double logBF_summary_cpp(
    const arma::vec&  XtY,
    const arma::mat&  XtX,
    double            YtY,
    const arma::vec&  sigma_vec,
    const arma::uvec& gamma_idx,
    int               n
) {
  const arma::uvec gamma0 = gamma_idx - 1;
  return logBF_internal(XtY, XtX, YtY, sigma_vec, gamma0, n);
}


///////// beta_summary_cpp //////////
//' @title Compute causal estimates for a single model (summary-data version)
//' @description Returns the posterior mean causal estimate vector for the
//' model \eqn{\gamma} defined by \code{gamma_idx} (column indices into betaX).
//' Estimates are computed from pre-computed summary-level statistics to avoid
//' recomputing matrix products. Exposures not in \eqn{\gamma} receive an estimate of zero.
//' Throughout, \eqn{d} denotes the number of exposures.
//'
//' @param XtY       Numeric vector of length \eqn{d}. Pre-computed
//'                  \eqn{X^\top Y}
//' @param XtX       Numeric matrix of dimension \eqn{d \times d}. Pre-computed
//'                  \eqn{X^\top X}
//' @param sigma_vec Numeric vector of length \eqn{d}. Prior standard
//'                  deviations for causal effects
//' @param gamma_idx Integer vector of column indices into betaX for model \eqn{\gamma}
//' @param d_Exposure Positive integer. Total number of exposures \eqn{d}
//'
//' @return Numeric vector of length \eqn{d}. Posterior mean causal estimates;
//' elements corresponding to exposures absent from \eqn{\gamma} are zero
//'
//' @export
// [[Rcpp::export]]
arma::vec beta_summary_cpp(
    const arma::vec&  XtY,
    const arma::mat&  XtX,
    const arma::vec&  sigma_vec,
    const arma::uvec& gamma_idx,
    int               d_Exposure
) {
  const arma::uvec gamma0 = gamma_idx - 1;
  return beta_internal(XtY, XtX, sigma_vec, gamma0, d_Exposure);
}


///////// cooksD_cpp //////////
//' @title Compute Cook's distance for a given MRBMA model
//' @description Computes Cook's distance of \insertCite{Cook1977;textual}{MRBMA}
//' for each SNP in model \eqn{\gamma} using a ridge-regularised hat matrix. The threshold is the 50th percentile of an \eqn{F(k, n -
//' k)} distribution, where \eqn{k = |\gamma|}
//'
//' @param y         Numeric vector of length \eqn{n}.
//'                  SNP-outcome association estimates (\code{betaY})
//' @param x         Numeric matrix of dimension \eqn{n \times |\gamma|}.
//'                  Submatrix of \code{betaX} restricted to exposures in \eqn{\gamma}
//' @param sigma_vec Numeric vector of length \eqn{|\gamma|}. Prior standard
//'                  deviations for the exposures in \eqn{\gamma}
//'
//' @return A named list:
//' \describe{
//'   \item{\code{cooksD}}{Numeric vector of length \eqn{n}:
//'         Cook's distance for each SNP}
//'   \item{\code{cooksD_thresh}}{Scalar double: suggested flagging threshold
//'         (\eqn{F_{0.5}(k, n-k)})}
//' }
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List cooksD_cpp(
    const arma::vec& y,
    const arma::mat& x,
    const arma::vec& sigma_vec
) {
  const int n = static_cast<int>(y.n_elem);  // number of SNPs
  const int k = static_cast<int>(x.n_cols);  // model size

  if (n <= k)
    Rcpp::stop("n_SNPs must exceed model size (n > ncol(x)).");

  const arma::mat sigma_inv = arma::diagmat(arma::pow(sigma_vec, -2.0));
  const arma::mat M         = x.t() * x + sigma_inv;
  const arma::mat H         = x * arma::solve(M, x.t());
  const arma::vec h_i       = H.diag();
  const arma::vec e         = y - H * y;
  const double    s_sq      = arma::as_scalar(e.t() * e) / (n - k);

  const arma::vec cd = arma::pow(e, 2.0) / (s_sq * k) %
                       (h_i / arma::pow(1.0 - h_i, 2.0));

  Rcpp::Function qf("qf");
  const double thresh = Rcpp::as<double>(
    qf(0.5, Rcpp::Named("df1") = k, Rcpp::Named("df2") = n - k)
  );

  return Rcpp::List::create(
    Rcpp::Named("cooksD")        = cd,
    Rcpp::Named("cooksD_thresh") = thresh
  );
}


///////// exhaustive_bf_cpp //////////
//' @title Exhaustive BF search -- evaluate all \eqn{2^d - 1} sub-models in C++
//' @description Enumerates all combinations of exposures of size
//' \eqn{k = 1, \ldots, d} using a recursive combinatorial generator and computes the log10 Bayes Factor (BF), log10 binomial prior, and
//' causal estimate vector for each model. All matrix operations are performed via RcppArmadillo. This function is the computational
//' backend of \code{MRBMA_exhaustive}. For the BF derivation see \insertCite{Zuber2020;textual}{MRBMA}.
//' Throughout, \eqn{d} denotes the number of exposures, \eqn{n} the number of SNPs, and \eqn{k = |\gamma|} the model size.
//'
//' @param XtY       Numeric vector of length \eqn{d}. Pre-computed
//'                  \eqn{X^\top Y}
//' @param XtX       Numeric matrix of dimension \eqn{d \times d}. Pre-computed
//'                  \eqn{X^\top X}
//' @param YtY       Positive numeric scalar. Pre-computed \eqn{Y^\top Y}
//' @param sigma_vec Numeric vector of length \eqn{d}. Prior standard
//'                  deviations
//' @param d_Exposure Positive integer. Number of exposures \eqn{d}
//' @param n_SNPs    Positive integer. Number of SNPs (instruments) \eqn{n}
//' @param prior_prob Numeric scalar in \eqn{(0, 1)}. Prior inclusion
//'                  probability per exposure
//'
//' @return A named list:
//' \describe{
//'   \item{\code{keys}}{Character vector of model keys (comma-separated
//'         column indices into betaX), length \eqn{2^d - 1}}
//'   \item{\code{logBF}}{Numeric vector of log10 Bayes Factors (BFs), same length}
//'   \item{\code{logprior}}{Numeric vector of log10 binomial prior values,
//'         same length}
//'   \item{\code{theta_mat}}{Numeric matrix of dimension
//'         \eqn{d \times (2^d - 1)}: causal estimate vectors for all models}
//' }
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List exhaustive_bf_cpp(
    const arma::vec& XtY,
    const arma::mat& XtX,
    double           YtY,
    const arma::vec& sigma_vec,
    int              d_Exposure,
    int              n_SNPs,
    double           prior_prob
) {
  std::vector<std::string> all_keys;
  std::vector<double>      all_logBF;
  std::vector<double>      all_logprior;
  std::vector<arma::vec>   all_theta;

  for (int k = 1; k <= d_Exposure; ++k) {
    const auto combos = gen_combos(d_Exposure, k);
    const double lp   = log_prior(k, d_Exposure, prior_prob);

    for (const auto& combo : combos) {
      arma::uvec gamma0(combo.size());
      for (std::size_t i = 0; i < combo.size(); ++i)
        gamma0[i] = combo[i] - 1;

      all_keys.push_back(make_key(combo));
      all_logBF.push_back(logBF_internal(XtY, XtX, YtY, sigma_vec, gamma0, n_SNPs));
      all_logprior.push_back(lp);
      all_theta.push_back(beta_internal(XtY, XtX, sigma_vec, gamma0, d_Exposure));
    }
  }

  const int n_models = static_cast<int>(all_keys.size());
  arma::mat theta_mat(d_Exposure, n_models);
  for (int j = 0; j < n_models; ++j)
    theta_mat.col(j) = all_theta[j];

  Rcpp::CharacterVector keys_out(n_models);
  arma::vec logBF_vec(n_models), logprior_vec(n_models);
  for (int j = 0; j < n_models; ++j) {
    keys_out[j]     = all_keys[j];
    logBF_vec[j]    = all_logBF[j];
    logprior_vec[j] = all_logprior[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("keys")      = keys_out,
    Rcpp::Named("logBF")     = logBF_vec,
    Rcpp::Named("logprior")  = logprior_vec,
    Rcpp::Named("theta_mat") = theta_mat
  );
}


///////// sss_core_cpp //////////
//' @title Shotgun Stochastic Search (SSS) core in C++ for MRBMA
//' @description Implements the shotgun stochastic search (SSS) algorithm of
//' \insertCite{Hans2007;textual}{MRBMA} as adapted for multivariable MR by \insertCite{Zuber2020;textual}{MRBMA}, entirely in C++,
//' using a \code{std::unordered_map} for O(1)-average memoisation of all visited models. This function is the computational backend of
//' \code{MRBMA_SSS}
//'
//' @param XtY       Numeric vector of length \eqn{d}. Pre-computed
//'                  \eqn{X^\top Y}
//' @param XtX       Numeric matrix of dimension \eqn{d \times d}. Pre-computed
//'                  \eqn{X^\top X}
//' @param YtY       Positive numeric scalar. Pre-computed \eqn{Y^\top Y}
//' @param sigma_vec Numeric vector of length \eqn{d}. Prior standard
//'                  deviations
//' @param d_Exposure      Positive integer. Number of exposures \eqn{d}
//' @param n_SNPs    Positive integer. Number of SNPs (instruments) \eqn{n}
//' @param prior_prob Numeric scalar in \eqn{(0, 1)}. Prior inclusion
//'                  probability per exposure
//' @param kmin      Positive integer \eqn{k_{\min}}. Minimum model size for the
//'                  deterministic phase and the SSS neighbourhood constraint
//' @param kmax      Positive integer \eqn{k_{\max} \leq d}. Maximum model size
//'                  for the SSS neighbourhood constraint
//' @param max_iter  Positive integer. Number of SSS iterations. Ignored if
//'                  \code{kmax == kmin}
//' @param verbose   Logical. If \code{TRUE}, prints progress to \code{Rcout}
//'                  every 1,000 iterations (default \code{false})
//'
//' @details
//' Throughout, \eqn{d} denotes the number of exposures, \eqn{n} the number of
//' SNPs, and \eqn{k = |\gamma|} the current model size (\eqn{k_{\min} \leq k \leq k_{\max}}).
//'
//' Deterministic phase: all models of size \eqn{1, \ldots, k_{\min}} are
//' enumerated and evaluated. Stochastic phase: at each iteration a model is
//' drawn from the current neighbourhood weighted by posterior evidence
//' (\eqn{\propto 10^{\text{logBF} + \text{logprior}}}); its neighbourhood
//' is created by add, delete and swap moves; unvisited neighbours are
//' evaluated and stored in the hash tables; the process repeats. Log-sum-exp
//' rescaling prevents \eqn{10^x} overflow during weighted sampling
//'
//' @return A named list:
//' \describe{
//'   \item{\code{keys}}{Character vector of all visited model keys}
//'   \item{\code{logBF}}{Numeric vector of log10 Bayes Factors (BFs) for visited
//'         models}
//'   \item{\code{logprior}}{Numeric vector of log10 prior values for visited
//'         models}
//'   \item{\code{theta_mat}}{Numeric matrix of dimension \eqn{d \times M},
//'         where \eqn{M} is the total number of visited models: causal
//'         estimate vectors}
//' }
//'
//' @references
//' \insertAllCited{}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List sss_core_cpp(
    const arma::vec& XtY,
    const arma::mat& XtX,
    double           YtY,
    const arma::vec& sigma_vec,
    int              d_Exposure,
    int              n_SNPs,
    double           prior_prob,
    int              kmin,
    int              kmax,
    int              max_iter,
    bool             verbose = false
) {
  std::unordered_map<std::string, double>    hashlogBF;
  std::unordered_map<std::string, double>    hashlogprior;
  std::unordered_map<std::string, arma::vec> hashTheta;
  hashlogBF.reserve(4096);
  hashlogprior.reserve(4096);
  hashTheta.reserve(4096);

  auto evaluate_config = [&](const std::vector<int>& cfg) {
    std::string key = make_key(cfg);
    if (hashlogBF.count(key)) return;

    arma::uvec gamma0(cfg.size());
    for (std::size_t i = 0; i < cfg.size(); ++i)
      gamma0[i] = cfg[i] - 1;

    hashlogBF[key]    = logBF_internal(XtY, XtX, YtY, sigma_vec, gamma0, n_SNPs);
    hashlogprior[key] = log_prior(static_cast<int>(cfg.size()), d_Exposure, prior_prob);
    hashTheta[key]    = beta_internal(XtY, XtX, sigma_vec, gamma0, d_Exposure);
  };

  // ---- Deterministic phase ----
  std::vector<std::vector<int>> current_nbhd;
  for (int k = 1; k <= kmin; ++k) {
    auto combos = gen_combos(d_Exposure, k);
    for (const auto& cfg : combos) evaluate_config(cfg);
    if (k == kmin) current_nbhd = combos;
  }

  // ---- Stochastic phase ----
  if (kmax > kmin) {
    Rcpp::RNGScope rng_scope;

    for (int iter = 0; iter < max_iter; ++iter) {

      const int nbhd_sz = static_cast<int>(current_nbhd.size());
      arma::vec log_ev(nbhd_sz);
      for (int j = 0; j < nbhd_sz; ++j) {
        const std::string key = make_key(current_nbhd[j]);
        log_ev[j] = hashlogBF.at(key) + hashlogprior.at(key);
      }

      const double max_ev = log_ev.max();
      arma::vec ev = log_ev;
      if (max_ev > 308.0) ev -= (max_ev - 308.0 + 1.0);
      ev = arma::pow(arma::vec(nbhd_sz, arma::fill::value(10.0)), ev);
      const double    sum_ev = arma::sum(ev);
      const arma::vec probs  = ev / sum_ev;

      const double u = R::runif(0.0, 1.0);
      double cumsum  = 0.0;
      int chosen     = nbhd_sz - 1;
      for (int j = 0; j < nbhd_sz; ++j) {
        cumsum += probs[j];
        if (u <= cumsum) { chosen = j; break; }
      }

      auto new_nbhd = create_neighbourhood(current_nbhd[chosen], d_Exposure, kmin, kmax);
      for (const auto& cfg : new_nbhd) evaluate_config(cfg);
      current_nbhd = std::move(new_nbhd);

      if (verbose && (iter + 1) % 1000 == 0) {
        Rcpp::Rcout << "[MRBMA_SSS] iter " << iter + 1
                    << " / " << max_iter
                    << " | models visited: " << hashlogBF.size()
                    << "\n";
      }
    }
  }

  // ---- Pack results ----
  const int n_visited = static_cast<int>(hashlogBF.size());
  Rcpp::CharacterVector all_keys(n_visited);
  arma::vec all_logBF(n_visited), all_logprior(n_visited);
  arma::mat theta_mat(d_Exposure, n_visited);

  int j = 0;
  for (const auto& kv : hashlogBF) {
    all_keys[j]      = kv.first;
    all_logBF[j]     = kv.second;
    all_logprior[j]  = hashlogprior.at(kv.first);
    theta_mat.col(j) = hashTheta.at(kv.first);
    ++j;
  }

  return Rcpp::List::create(
    Rcpp::Named("keys")      = all_keys,
    Rcpp::Named("logBF")     = all_logBF,
    Rcpp::Named("logprior")  = all_logprior,
    Rcpp::Named("theta_mat") = theta_mat
  );
}
