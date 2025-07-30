#' bayesrank_2step: Two-Step Bias-Aware Bayesian Rank Aggregation
#'
#' This function implements a two-stage version of the BayeSRank algorithm.
#' The first stage runs a preliminary MCMC chain with diffuse priors to
#' empirically tune the prior bounds and hyperparameters. These estimates
#' are then used in the second stage, which runs a full posterior inference
#' using the calibrated priors. This two-step procedure improves stability and
#' reduces manual tuning.
#'
#' @param rank_matrix An \eqn{n \times n} matrix of observed rankings, where
#'        \code{rank_matrix[i, j]} denotes the rank assigned by ranker \code{j}
#'        to item \code{i}. Typically, the diagonal represents self-rankings.
#' @param M_itrns Integer. Number of MCMC iterations per stage (pilot and main).
#' @param m_burn Integer. Number of burn-in iterations to discard for posterior summaries.
#' @param init_W Optional \eqn{n \times n} matrix. If provided, used as the initial latent
#'        utility matrix \code{W}. If \code{NULL}, an initial matrix is generated and reordered
#'        according to the observed rankings.
#' @param step1_minbeta,step1_maxbeta Numeric scalars. Truncation bounds for bias parameters \code{beta_j}
#'        in the pilot step.
#' @param step1_Var_beta0 Numeric scalar. Initial variance for \code{beta_j} in the pilot step.
#' @param step1_Var_epsilon0 Numeric scalar. Initial variance of error terms \code{epsilon_j} in the pilot step.
#' @param step1_a,step1_b,step1_c,step1_d Numeric scalars. Hyperparameters for inverse-gamma priors in the pilot step.
#' @param seed Optional integer. If set, ensures reproducibility of the MCMC runs.
#' @param verbose Logical. If \code{TRUE}, prints progress messages for each step.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{preliminary}}{Posterior output from the pilot step with diffuse priors.}
#'   \item{\code{final}}{Posterior output from the main chain using tuned priors.}
#'   \item{\code{tuned_hyper}}{List of calibrated hyperparameters passed to the second step.}
#' }
#'
#' @details
#' The first step of \code{bayesrank_2step()} estimates posterior summaries using
#' weak priors. The output is then used to:
#' \itemize{
#'   \item Estimate the posterior mean and standard deviation of \code{beta_j},
#'         which define a 3-sigma truncation range for the second step;
#'   \item Estimate the shape and rate of inverse-gamma priors on variance components
#'         via moment-matching.
#' }
#'
#' These empirical values are then passed into a second run of \code{\link{BayeSRank}}.
#'
#' @seealso \code{\link{BayeSRank}} for the underlying single-chain rank aggregation algorithm.
#'
#' @examples
#' set.seed(42)
#' n <- 5
#' rank_matrix <- matrix(sample(1:n, n*n, replace = TRUE), n, n)
#'
#' out <- bayesrank_2step(
#'   rank_matrix = rank_matrix,
#'   M_itrns = 1000,
#'   m_burn = 200,
#'   verbose = TRUE,
#'   seed = 123
#' )
#'
#' str(out$final$post_mean_mus)
#' out$tuned_hyper$minBeta
#'
#' @export
bayesrank_2step <- function(rank_matrix,
                            M_itrns      = 5000,
                            m_burn       = 1000,
                            init_W       = NULL,
                            # ↓ optional pilot‑stage hyper‑priors
                            step1_minbeta      = -100,
                            step1_maxbeta      =  100,
                            step1_Var_beta0    = 1,
                            step1_Var_epsilon0 = 1,
                            step1_a = 0.1, step1_b = 0.1,
                            step1_c = 0.1, step1_d = 0.1,
                            # ↓ random‑seed helper so runs are reproducible
                            seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  ## --------------------------------------------------------------------------
  ## 0.  Build or accept an initial W  (un‑ordered latent scores)
  ## --------------------------------------------------------------------------
  n_teams <- nrow(rank_matrix)

  if (is.null(init_W)) {
    init_W <- matrix(rnorm(n_teams^2, mean = 0, sd = 1), nrow = n_teams)
  }

  ## Re–order each column of W so that larger latent utility ⇒ higher rank
  sorted_W <- apply(init_W, 2, sort, decreasing = TRUE)
  initW <- mapply(function(sorted_col, rank_col) sorted_col[rank_col],
                  as.data.frame(sorted_W),
                  as.data.frame(rank_matrix))
  initW <- matrix(unlist(initW), nrow = n_teams)

  ## --------------------------------------------------------------------------
  ## 1.  Preliminary chain with weak / diffuse priors
  ## --------------------------------------------------------------------------
  prelim <- BayeSRank(n = n_teams,
                      r = rank_matrix,
                      M = M_itrns,
                      burnin = m_burn,
                      W = initW,
                      mu0         = rnorm(n_teams, 0, 1),
                      mubeta0     = rnorm(1, 0, 1),
                      Var_beta0   = step1_Var_beta0,
                      Var_epsilon0= step1_Var_epsilon0,
                      minbeta     = step1_minbeta,
                      maxbeta     = step1_maxbeta,
                      a = step1_a, b = step1_b,
                      c = step1_c, d = step1_d)

  ## Empirical tuning from pilot output
  Beta.bar <- mean(prelim$Beta)
  Beta.sd  <-  sd(prelim$Beta)
  minBeta_est <- Beta.bar - 3 * Beta.sd
  maxBeta_est <- Beta.bar + 3 * Beta.sd
  varBeta_est    <- median(prelim$Var_beta)
  varEpsilon_est <- median(prelim$Var_epsilon)

  ## Moment‑match Γ(shape,rate) priors for 1/σ²
  if (!exists("find_shape_rate", mode = "function")) {
    find_shape_rate <- function(mean, sd) {
      shape <- (mean^2) / (sd^2)
      rate  <-  mean    / (sd^2)
      list(shape = shape, rate = rate)
    }
  }
  IG_Beta_est    <- find_shape_rate(varBeta_est,    sd(prelim$Var_beta))
  IG_Epsilon_est <- find_shape_rate(varEpsilon_est, sd(prelim$Var_epsilon))


  ## --------------------------------------------------------------------------
  ## 2.  Main chain with tuned hyper‑priors
  ## --------------------------------------------------------------------------
  res <- BayeSRank(n = n_teams,
                   r = rank_matrix,
                   M = M_itrns,
                   burnin = m_burn,
                   W = initW,
                   mu0         = rnorm(n_teams, 0, 1),
                   mubeta0     = rnorm(1, 0, 1),
                   Var_beta0   = varBeta_est,
                   Var_epsilon0= varEpsilon_est,
                   minbeta     = minBeta_est,
                   maxbeta     = maxBeta_est,
                   a = IG_Beta_est$shape, b = IG_Beta_est$rate,
                   c = IG_Epsilon_est$shape, d = IG_Epsilon_est$rate)


    return(res)
}

