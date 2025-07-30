#' BayeSRank: Bias-Aware Bayesian Aggregation in Peer and Self Ranking
#'
#' Implements the core Gibbs sampling algorithm for the BayeSRank method proposed in
#' Ruan et al. (XXXX). This method estimates latent performance scores and biases
#' from noisy ranking data, while preserving ordinal information and accounting for
#' potential self-evaluation bias and varying levels of ranking quality.
#'
#' This version performs a single-chain run using user-specified priors and does not include
#' the two-step empirical Bayes procedure. For automatic hyperparameter tuning, see
#' [bayesrank_2step()].
#'
#' @param n Integer. Number of items/teams to be ranked (equal to the number of rankers).
#' @param r An \eqn{n \times n} matrix of observed rankings where \code{r[i, j]} denotes the rank
#'        assigned by ranker \eqn{j} to item \eqn{i}. Diagonal entries represent
#'        self-evaluations.
#' @param M Integer. Number of MCMC iterations.
#' @param burnin Integer. Number of burn-in samples to discard when summarizing posteriors.
#' @param W Numeric \eqn{n \times n} matrix. Initial latent score matrix used in data augmentation.
#' @param mu0 Numeric vector of length \code{n}. Initial values for latent team-level performance \eqn{\mu_i}.
#' @param mubeta0 Numeric scalar. Initial value for the mean of evaluator bias parameters.
#' @param beta0 Not used. Placeholder retained for compatibility.
#' @param Var_beta0 Numeric scalar. Initial variance for evaluator bias \eqn{\beta_j}.
#' @param Var_epsilon0 Numeric scalar or vector of length \code{n}. Initial error variance.
#' @param minbeta,maxbeta Numeric scalars. Lower and upper bounds for truncated normal prior
#'        on evaluator bias parameters.
#' @param a,b,c,d Numeric scalars. Hyperparameters for the inverse-gamma priors:
#'        \itemize{
#'          \item \code{a}, \code{b}: shape and rate for \eqn{\text{Var}(\beta_j)}
#'          \item \code{c}, \code{d}: shape and rate for \eqn{\text{Var}(\epsilon_j)}
#'        }
#'
#' @return A list with posterior summaries and evaluation metrics:
#' \describe{
#'   \item{\code{Beta}}{Posterior samples of bias terms \eqn{\beta_j} after burn-in.}
#'   \item{\code{Var_beta}}{Posterior samples of bias variance.}
#'   \item{\code{Var_epsilon}}{Posterior samples of error variance.}
#'   \item{\code{post_mean_mus}}{Posterior means of latent team performance \eqn{\mu_i}.}
#'   \item{\code{post_mean_mubeta}}{Posterior mean of bias prior mean \eqn{\mu_\beta}.}
#'   \item{\code{post_median_Varbeta}}{Posterior median of \eqn{\text{Var}(\beta)}.}
#'   \item{\code{corr_spearman}}{Spearman correlation between estimated ranks and true ranks.}
#'   \item{\code{top1rate}}{Proportion of simulations correctly identifying the top team.}
#'   \item{\code{top3rate}}{Proportion correctly identifying the top-3 teams.}
#' }
#'
#' @seealso \code{\link{bayesrank_2step}} for a two-step version with automatic prior tuning.
#'
#' @references
#' Ruan, H., Wang, K., & Wang, X. (XXXX). *Bias meets Bayes:  A Bayesian Aggregation Method for Peer and Self Ranking*
##' @examples
#' # Simulate toy ranking data
#' set.seed(123)
#' n <- 5
#' true_rank <- 1:n
#' rank_matrix <- matrix(sample(1:n, n * n, replace = TRUE), n, n)
#'
#' initW <- matrix(rnorm(n^2), n, n)
#' sorted_W <- apply(initW, 2, sort, decreasing = TRUE)
#' initW <- mapply(function(sorted_col, rank_col) sorted_col[rank_col],
#'                 as.data.frame(sorted_W),
#'                 as.data.frame(rank_matrix))
#' initW <- matrix(unlist(initW), nrow = n)
#'
#' res <- BayeSRank(
#'   n = n,
#'   r = rank_matrix,
#'   M = 1000,
#'   burnin = 200,
#'   W = initW,
#'   mu0 = rnorm(n),
#'   mubeta0 = rnorm(1),
#'   beta0 = NULL,
#'   Var_beta0 = 1,
#'   Var_epsilon0 = 1,
#'   minbeta = -10,
#'   maxbeta = 10,
#'   a = 0.1, b = 0.1,
#'   c = 0.1, d = 0.1
#' )
#' str(res$post_mean_mus)
#'
#' @export

BayeSRank <- function(n,
                      r,
                      M,
                      burnin,
                      W,             # W=initial Omega matrix
                      mu0,
                      mubeta0,
                      beta0,
                      Var_beta0,
                      Var_epsilon0,
                      minbeta,
                      maxbeta,
                      a,b,c,d
){

  # Dimensions of the rank matrix
  N <- n  # Number of teams (rows of `r`)
  J <- n  # Number of rankers (columns of `r`)

  # Initialize matrices to store Mu, Beta, and Sigma_s2
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M) # DO NOT USE 0 coz NAs will be produced at 2nd itration

  # Initial values
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0

  beta0 <- rnorm(J, mubeta0, sd=sqrt(Var_beta0) )
  Beta[, 1] <- beta0

  Var_epsilon[, 1] <- Var_epsilon0

  for (m in 2:M) {

    # if (m %% verbose == 0) {
    #   print(m)  # Print progress
    # }

    #   # 1.Update Mu

    for (i in 1:N) { #for each team in ROWS

      indBeta <- Beta[, m-1]
      indBeta[-i] <-0

      g <- sum( (W[i, ] - indBeta) / Var_epsilon[, m-1] ) #numerator g
      h <- sum(1 / Var_epsilon[, m-1]) + 1                #denominator h

      Mu[i, m] <- rnorm(1, mean = (g/h), sd = sqrt(1/h) )
    }


    #   # 2.Update Beta

    for (j in 1:J) { #for each ranker in COLUMNS

      g <- ( (W[j, j] - Mu[j, m]) / Var_epsilon[j, m-1] ) + (Mu_beta[m-1]/Var_beta[m-1]) #numerator g
      h <- (1 / Var_epsilon[j, m-1]) + (1/Var_beta[m-1])              #denominator h

      # Beta[j, m] <- rnorm(1, (g/h), sqrt(1/h))

      Beta[j, m] <- truncnorm::rtruncnorm(1, a = minbeta, b = maxbeta,
                                          mean = (g/h),
                                          sd = sqrt(1/h))
    }


    #   # 3.Update mu_beta

    # g <- sum(Beta[,m]/Var_beta[m-1])                                      #numerator g
    # h <- (1 / Var_beta[m-1]) + (1/tau)^2 #denominator h
    #
    # Mu_beta[m] <- rnorm(1, (g/h), sqrt(1/h))

    Mu_beta[m] <- truncnorm::rtruncnorm(1, a = minbeta, b = maxbeta,
                                        mean = mean(Beta[, m]),
                                        sd = sqrt(Var_beta[m-1]/J) )

    #   # 4.Update Var_beta
    # Var_beta[m] <- rinvgamma_trunc(1,
    #                                shape= a + (J/2),
    #                                scale = b + sum( (Beta[, m]-Mu_beta[m])^2 )/2,
    #                                lower = sig2_lower,
    #                                upper = sig2_upper)

    Var_beta[m] <-  1/rgamma(1,
                               shape = a+ (J/2),
                               rate = b+ sum( (Beta[, m]-Mu_beta[,m])^2 )/2 )


    # Var_beta[m] <-  rinvgamma(1,
    #                            alpha = a+ (J/2),
    #                            beta =  b+ sum( (Beta[, m]-Mu_beta[,m])^2 )/2 )
    #

    #   # 5.Update Var_epsilon

    for (j in 1:J) {

      indBeta <- Beta[, m]
      indBeta[-j] <-0

      # Var_epsilon[j, m] <- rinvgamma_trunc(1,
      #                                      shape= c+ (J/2),
      #                                      scale= d + sum( (W[, j] - Mu[, m] - indBeta)^2 )/2,
      #                                      lower= sig2_lower,
      #                                      upper= sig2_upper)

      Var_epsilon[j, m] <- 1/rgamma(1,
                                      shape = c+ (J/2),
                                      rate = d+ sum( (W[, j] - Mu[,m] - indBeta)^2 )/2 )

      # Var_epsilon[j, m] <- rinvgamma(1,
      #                                 alpha = c+ (J/2),
      #                                 beta = d+ sum( (W[, j] - Mu[,m] - indBeta)^2 )/2 )
      #
      }


    #   # 6. Update W

    for (j in 1:J) {
      for (i in 1:N) {

        if (r[i, j] == 1) {
          upper <- Inf
          lower <- W[which(r[, j] == 2), j]

        } else if (r[i,j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -Inf

        } else {

          lower <- W[which(r[, j] == r[i,j] + 1), j]
          upper <- W[which(r[, j] == r[i,j] - 1), j]

        }
        W[i,j] <- truncnorm::rtruncnorm(1, a = lower, b = upper,
                                        mean = Mu[i, m] + Beta[j]*(i==j),
                                        sd = sqrt( Var_epsilon[j, m] ))
      }
    }
  }

  drop = 1:burnin
  post.mean.mus <- apply(Mu[,-drop], 1, mean)
  cor <- cor(rank(-post.mean.mus), true_rank, use = "everything", method="spearman")

  post.mean.betas <- apply(Beta[,-drop], 1, mean)

  post.mean.mu_beta <- mean(Mu_beta[,-drop])
  post.median.Var_beta <- median(Var_beta[,-drop])


  ssq <- calculate_ssq_diff(ranked_list=rank(-post.mean.mus), true_list=true_rank)
  top1rate <- calculate_top1(ranked_list=rank(-post.mean.mus), true_list=true_rank)
  top3rate <- calculate_top3(ranked_list=rank(-post.mean.mus), true_list=true_rank)

  out <- list(# mu = Mu,
    Beta = Beta[,-drop],
    #mu_beta = Mu_beta,
    Var_beta = Var_beta[,-drop],
    Var_epsilon = Var_epsilon[,-drop],
    # W=W,
    post_mean_mus = post.mean.mus,
    # post_mean_betas = post.mean.betas,
    post_mean_mubeta = post.mean.mu_beta,
    post_median_Varbeta = post.median.Var_beta,
    corr_spearman = cor,
    top1rate=top1rate,
    top3rate=top3rate
    # post_rank = rank(-post.mean.mu),
    # true_rank = true_rank,
    # true_beta = mydat$true_beta
  )

  return(out)
}

