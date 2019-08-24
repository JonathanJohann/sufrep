sample_latent_group <- function(n, k) {
  latent <- sample.int(k, size = n, replace = T)
  latent
}

# alpha[l] ~ Laplace(1)
sample_alpha <- function(latent) {
  k <- length(unique(latent))
  sgn <- sample(c(-1, 1), replace = T, size = k)
  base_alpha <- sgn * rexp(sqrt(2), n = k)
  alpha <- base_alpha[latent]
  alpha
}


sample_x_mean <- function(k, p, shifts = 3) {
  mu <- matrix(0, nrow = k, ncol = p)
  for (i in 1:k) {
    shift <- sample(c(-1, 1), size = shifts, replace = TRUE)
    indices <- sample(p, size = shifts, replace = FALSE)
    mu[i, indices] <- shift
  }
  mu
}


sample_x <- function(latent, x_mean, rho = 0.5) {
  p <- dim(x_mean)[2]
  sigma <- toeplitz(rho^seq(0, (p - 1)))
  x <- t(apply(x_mean[latent, ], 1, function(mul) MASS::mvrnorm(mu = mul, Sigma = sigma)))
  colnames(x) <- paste("X", 1:p, sep = "")
  x
}


sample_beta_latent <- function(k, p) {
  if (p < 3) stop("Function sample_beta_latent required p > 2.")
  base_vec <- rep(c(1, 0, -1), p)
  beta <- replicate(k, sample(base_vec, size = p, replace = F))
  beta <- t(apply(beta, 2, function(x) x / sqrt(sum(x^2))))
  beta
}


sample_beta_global <- function(p) {
  sample_beta_latent(k = 1, p = p)
}


sample_observed_category <- function(latent, n, k, ngl, pl) {
  if (pl <= 0.5) stop("Argument pl must be > 0.5.")
  num_cats <- k * ngl
  all_cats <- seq(num_cats)
  latent_to_obs_cat <- matrix(sample(all_cats, num_cats), k, ngl) # scrambling the map
  own_latent <- rbinom(prob = pl, n = n, size = 1)
  obs_cat <- rep(0, n)
  for (i in 1:n) {
    own_cats <- latent_to_obs_cat[latent[i], ]
    if (own_latent[i]) {
      obs_cat[i] <- sample(own_cats, size = 1)
    } else {
      obs_cat[i] <- sample(setdiff(all_cats, own_cats), size = 1)
    }
  }
  obs_cat
}


linear_global_response <- function(alpha, x, beta) {
  if (dim(beta)[1] > 1) stop("Argument beta must have only one row.")
  n <- length(alpha)
  y <- alpha + x %*% t(beta) + rnorm(n)
  c(y)
}


linear_latent_response <- function(latent, alpha, x, betas) {
  if (dim(betas)[1] == 1) stop("Argument beta must multiple rows.")
  n <- length(alpha)
  y <- alpha + apply(x * betas[latent, ], 1, sum) + rnorm(n)
  y
}

# n: number of observations
# p: number of continuous covariates
# k: number of latent groups
# ngl: number of observable cats per latent group
# pl: probability that observable category belongs to their 'own' latent group
create_data <- function(n, p, k, ngl, pl, type = "global") {

  # Draw the latent group L[i] and observable category G[i]
  latent <- sample_latent_group(n, k)
  obs_cat <- sample_observed_category(latent, n, k, ngl, pl)

  # Draw intercept shift alpha[i]
  alpha <- sample_alpha(latent)

  # Draw continuous covariatex X[i]
  x_mean <- sample_x_mean(k, p)
  x <- sample_x(latent, x_mean)

  # Compute response
  if (type == "global") {
    beta <- sample_beta_global(p)
    y <- linear_global_response(alpha, x, beta)
  } else if (type == "latent") {
    beta <- sample_beta_latent(k, p)
    y <- linear_latent_response(latent, alpha, x, beta)
  }

  g <- factor(obs_cat)

  # Put it all together
  data <- list(x = x, y = y, l = latent, alpha = alpha, g = g)
  data
}
