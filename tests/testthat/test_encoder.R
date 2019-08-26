library(sufrep)
library(testthat)


test_that("all encoders have consistent encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y
  for (method in methods) {
    enc <- make_encoder("permutation", X, G, Y = Y)
    X_enc <- enc(X, G)
    expect_true(has_consistent_encoding(X, G, X_enc))
  }
})


test_that("classical encodings have same R-squared", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y
  rsquared <- rep(0, length(classical_methods))
  for (i in seq_along(classical_methods)) {
    method <- classical_methods[i]
    enc <- make_encoder(method, X, G)
    X_enc <- enc(X, G)
    ols <- lm(Y ~ X_enc)
    rsquared[i] <- summary(ols)$r.squared
    expect_equal(rsquared[1], rsquared[i], tol = 1e-5)
  }
})


test_that("means encoding correctly takes means", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y
  enc <- make_encoder(method = "means", X, G)
  X_enc <- enc(X, G)
  G_levels <- unique(G)
  for (g in seq_along(G_levels)) {
    idx <- G == g
    G_mean1 <- apply(X_enc[idx, (p + 1):dim(X_enc)[2], drop = F], 2, mean)
    G_mean2 <- apply(X[idx,], 2, mean)
    expect_equal(G_mean1, G_mean2, tol = 1e-8)
  }
})


test_that("low rank encodings have correct number of columns", {
  data <- make_data()
  X <- data$X
  G <- data$G
  p <- dim(X)[2]
  for (method in c("low_rank", "sparse_low_rank")) {
    for (q in seq(p)) {
      X_enc <- make_encoder(method = method, X, G, num_components = q)(X, G)
      expect_equal(dim(X_enc)[2], p + q)
    }
  }
})


test_that("low rank encodings break if number of components is too large", {
  data <- make_data()
  X <- data$X
  G <- data$G
  p <- dim(X)[2]
  for (method in c("low_rank", "sparse_low_rank")) {
      expect_error(make_encoder(method = method, X, G, num_components = p + 1))
  }
})


test_that("multi_permutation method has correct number of columns", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y
  p <- dim(X)[2]
  for (q in seq(p)) {
    X_enc <- make_encoder(method = "multi_permutation", X, G, num_permutations = q)(X, G)
    expect_equal(dim(X_enc)[2], p + q)
  }
})


test_that("continuous variable names are kept untouched", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y
  p <- dim(X)[2]

  colnames(X) <- letters[1:p]
  Xd <- as.data.frame(X)
  Xm <- as.matrix(X)

  for (method in methods) {
    X_enc_d <- make_encoder(method, X = Xd, G = G, Y = Y, num_permutations = q)(Xd, G)
    X_enc_m <- make_encoder(method, X = Xm, G = G, Y = Y, num_permutations = q)(Xm, G)
    expect_equal(colnames(X_enc_d)[1:p], letters[1:p])
    expect_equal(colnames(X_enc_m)[1:p], letters[1:p])
  }
})


test_that("output type matches input X type", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y

  for (method in methods) {
    X_enc_d <- make_encoder(method, X = Xd, G = G, Y = Y, num_permutations = q)(Xd, G)
    X_enc_m <- make_encoder(method, X = Xm, G = G, Y = Y, num_permutations = q)(Xm, G)
    expect_equal(class(X_enc_m), "matrix")
    expect_equal(class(X_enc_d), "data.frame")
  }
})


test_that("raises an error if test data has an unseen level", {
  data <- make_data()
  X <- data$X
  G <- data$G
  Y <- data$Y

  Gnew <- sample(c(levels(G), "newlevelzz"), replace = T, size = dim(X)[1])
  for (method in methods) {
    enc = make_encoder(method, X = X, G = G, Y = Y, num_permutations = q)
    expect_error(enc(Xd, Gnew))
  }
})



