one_hot_encode <- function(num_categ) {
  CM <- data.frame(diag(num_categ))
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
  colnames(CM) <- paste("E", 1:(num_categ - 1), sep = "")
  return(CM)
}

fisher_encode <- function(X, G, Y) {
  CM <- aggregate(X[, Y], list(X[, G]), mean)
  colnames(CM) <- c("label", "X1")
  CM$X1 <- rank(CM$X1)
  ordering <- data.frame(unique(X[, G]))
  colnames(ordering) <- "ORD"
  CM <- CM[order(ordering$ORD), ]
  CM <- data.frame(CM$X1)
  return(CM)
}


means_encode <- function(X, G) {
  CM <- aggregate(X, list(X[, G]), mean)
  return(CM)
}

low_rank_encode <- function(X, G, num_components) {
  CM <- means_encode(X, G)
  CM <- as.data.frame(CM)
  for (i in 1:dim(CM)[2]) {
    CM[, i] <- as.numeric(CM[, i])
  }
  CM <- as.matrix(CM)
  decomp <- tryCatch({
    svd(CM)
  },
  error = function(e) {
    return(svd(CM[, colSums(!is.finite(CM)) == 0]))
  }
  )

  np <- min(num_components, dim(decomp$u)[2])
  CM <- data.frame(decomp$u[, 1:np])
  return(CM)
}

sparse_low_rank_encode <- function(X, G, num_components) {
  CM <- means_encode(X, G)
  CM <- as.matrix(CM)
  for (i in 1:dim(CM)[2]) {
    CM[, i] <- as.numeric(CM[, i])
  }
  CM <- as.matrix(CM)
  decomp <- tryCatch({
    sparsepca::spca(CM, verbose = FALSE)
  },
  error = function(e) {
    return(sparsepca::spca(CM[, colSums(!is.finite(CM)) == 0]))
  })

  np <- min(num_components, dim(decomp$loadings)[2])
  U <- decomp$loadings[, 1:np]
  CM <- as.matrix(CM)
  CM <- tryCatch({
    CM %*% U
  },
  error = function(e) {
    return(CM[, colSums(!is.finite(CM)) == 0] %*% U)
  }
  )
  CM <- data.frame(CM)
  return(CM)
}


permutation_encode <- function(num_categ, num_permutations = 1) {
  CM <- data.frame(sample(num_categ, size = num_categ, replace = FALSE))
  if (num_permutations > 1) {
    for (i in 2:num_permutations) {
      CM <- cbind(CM, data.frame(sample(num_categ, size = num_categ, replace = FALSE)))
    }
  }
  colnames(CM) <- paste("E", 1:num_permutations, sep = "")
  return(CM)
}


mnl_encode <- function(X, G, num_folds = 3) {

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
  thetahat <- data.frame(dense_coef)
  colnames(thetahat) <- names(dense_coef)
  rownames(thetahat) <- NULL
  return(CM)
}


validate_X <- function(X) {
  if (!all(sapply(X, is.numeric))) {
    stop("Argument X contains columns that are not numeric.")
  }
}


validate_G < function(X) {

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
encoder <- function(method, X, G,
                    prefix = "E",
                    Y = NULL,
                    num_components = NULL,
                    num_permutations = NULL,
                    num_folds = 3) {

  validate_X(X)
  validate_G(G)
  validate_options(method, Y, num_components, num_permutations, num_folds)

  id <- data.frame(G=unique(G))
  num_categ <- dim(id)[1]
  n <- dim(X)[1]
  p <- dim(X)[2]

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
  colnames(CM) <- paste(prefix, 1:dim(CM)[2], sep = "")
  map <- data.frame(cbind(id, CM))


  f <- function(X) {
    X_enc <- map[X[, G], ]
    rownames(X_enc) <- NULL
    X <- cbind(X, X_enc)
    X <- X[, -which(colnames(X) %in% c(G))]
    return(X)
  }

  return(f)
}
