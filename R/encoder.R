one_hot_encode <- function(num_categ) {
  CM <- diag(num_categ)[,1:(num_categ-1)]
  return(CM)
}

helmert_encode <- function(num_categ) {
  CM <- stats::contr.helmert(num_categ)
  return(CM)
}

deviation_encode <- function(num_categ) {
  CM <- stats::contr.sum(num_categ)
  return(CM)
}

repeated_effect_encode <- function(num_categ) {
  TH <- matrix(0, nrow = num_categ, ncol = num_categ)
  TH <- matrix((num_categ - col(TH)) / num_categ, nrow = num_categ, ncol = num_categ)
  TH[lower.tri(TH)] <- 0
  BH <- matrix(-col(TH) / num_categ, nrow = num_categ, ncol = num_categ)
  BH[upper.tri(BH)] <- 0
  diag(BH) <- 0
  CM <- TH + BH
  CM <- CM[, 1:(num_categ - 1)]
  return(CM)
}

difference_encode <- function(num_categ) {
  CM <- matrix(0, nrow = num_categ, ncol = num_categ)
  CM <- matrix(-1 / (col(CM) + 1), nrow = num_categ, ncol = num_categ)
  CM[lower.tri(CM)] <- 0
  CM <- CM[, 1:(num_categ - 1)]
  CM[row(CM) == (col(CM) + 1)] <- -apply(CM, 2, sum)
  return(CM)
}

simple_effect_encode <- function(num_categ) {
  CM <- matrix(-1 / num_categ, nrow = num_categ, ncol = num_categ)
  CM <- CM + diag(num_categ)
  CM <- CM[, 1:(num_categ - 1)]
  return(CM)
}

fisher_encode <- function(X, G, Y) {
  CM <- aggregate(Y, list(G), mean)
  colnames(CM) <- c("label", "X1")
  CM$X1 <- rank(CM$X1)
  ordering <- data.frame(unique(G))
  colnames(ordering) <- "ORD"
  CM <- CM[order(ordering$ORD), ]
  CM <- as.matrix(CM$X1)
  return(CM)
}


means_encode <- function(X, G) {
  CM <- as.matrix(aggregate(X, list(X[, G]), mean))
  return(CM)
}


low_rank_encode <- function(X, G, num_components) {
  if (num_components > dim(X)[2]) {
    stop("Argument num_components cannot be larger than number of X columns.")
  }
  CM <- means_encode(X, G)
  decomp <- svd(CM)
  CM <- as.matrix(decomp$u[, 1:num_components])
  return(CM)
}


sparse_low_rank_encode <- function(X, G, num_components) {
  if (num_components > dim(X)[2]) {
    stop("Argument num_components cannot be larger than number of X columns.")
  }
  CM <- means_encode(X, G)
  decomp <- sparsepca::spca(CM, verbose = FALSE)
  U <- decomp$loadings[, 1:num_components]
  CM <-  CM %*% U
  return(CM)
}


permutation_encode <- function(num_categ, num_perms) {
  CM <- replicate(num_perms, sample(num_categ, size = num_categ, replace = FALSE))
  return(CM)
}


mnl_encode <- function(X, G, num_folds) {

  # Fit glmnet using stratified K-fold cv
  df <- data.frame(X=X, G=G)
  cv_index <- caret::createFolds(factor(training$G), num_folds)
  tc <- caret::trainControl(index = cv_index, method = 'cv', number = num_folds)
  glmnet_fit <- caret::train(G ~ ., data = df, method = "glmnet",
                             trControl = tc, family= "multinomial")
  model <- glmnet_fit$finalModel
  sparse_coef <- coef(model, s=model$lambdaOpt)

  # Retrieve coefficients
  dense_coef <- lapply(sparse_coef, as.matrix)
  thetahat <- as.matrix(t(as.data.frame(dense_coef)))
  rownames(thetahat) <- NULL
  colnames(thetahat) <- NULL
  return(CM)
}


validate_X <- function(X) {
  if (!all(sapply(X, is.numeric))) {
    stop("Argument X contains columns that are not numeric.")
  }
}


validate_G <- function(G) {
  if (!is.factor(G)) {
    stop("Argument G must be factor.")
  }
}

validate_options <- function(method, Y, num_components, num_permutations, num_folds) {
  if ((method != "fisher") && (!is.null(Y))) {
    stop("Only method 'fisher' requires Y.")
  }
  if (!(method %in% c("low_rank", "sparse_low_rank")) && (!is.null(Y))) {
    stop("Only methods 'low_rank' and 'sparse_low_rank' requires num_components.")
  }
  if ((method != "multi_permutation") && (!is.null(Y))) {
    stop("Only method 'multi_permutation' requires num_permutations.")
  }
  # etc.
}


#' @export
make_make_encoder <- function(method, X, G,
                         prefix = "E",
                         Y = NULL,
                         num_components = NULL,
                         num_permutations = NULL,
                         num_folds = 3) {

  # Type validation
  validate_X(X)
  validate_G(G)

  # Argument validation
  validate_options(method, Y, num_components, num_permutations, num_folds)

  # Compute encoding
  num_categ <- length(unique(G))
  CM <- switch(method,
    one_hot = one_hot_encode(num_categ),
    helmert = helmert_encode(num_categ),
    deviation = deviation_encode(num_categ),
    repeated_effect = repeated_effect_encode(num_categ),
    difference = difference_encode(num_categ),
    simple_effect = simple_effect_encode(num_categ),
    fisher = fisher_encode(X, G, Y),
    means = means_encode(X, G),
    low_rank = low_rank_encode(X, G, num_components),
    sparse_low_rank = sparse_low_rank_encode(X, G, num_components),
    permutation = permutation_encode(num_categ, num_permutations = 1),
    multi_permutation = permutation_encode(num_categ, num_permutations),
    mnl = mnl_encode(X, G, num_folds)
  )

  # Create encoding function
  encoding_fun <- function(X, G) {
    validate_X(X)
    validate_G(G)

    # Augment original matrix
    X_aug <- cbind(X, CM[G,])

    # Maintain X type
    if (is.data.frame(X)) {
      X_aug <- as.data.frame(X_aug)
    }

    # Maintain row and column names, if appropriate
    rownames(X_aug) <- rownames(X)
    if (!is.null(colnames(X))) {
      colnames(X_aug) <- c(colnames(X), paste(prefix, 1:dim(CM)[2], sep = ""))
    }

    return(X_aug)
  }

  return(encoding_fun)
}
