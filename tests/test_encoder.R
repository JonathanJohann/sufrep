library(sufrep)
library(testthat)


# Checks that all the encodings for one category are the same,
#  regardless of the method used.
has_consistent_encoding <- function(X, G, X_enc) {
  G_num <- as.integer(G)
  p <- dim(X)[2]
  p_enc <- dim(X_enc)[2]
  is_consistent <- sapply(unique(G_num), function(g) {
    idx <- G_num == g
    encoded <- X_enc[idx, (p + 1):p_enc, drop = F]
    apply(encoded, 2, function(x) {
      all(x == x[1])
    })
  })
  all(is_consistent)
}


make_data <- function(n = 200, p = 4, m = 6) {
  cats <- apply(expand.grid(letters, letters), 1, function(s) paste0(s, collapse = ""))[1:m]
  G <- factor(sample(cats, replace = T, size = n))
  G_num <- as.integer(G)
  Y <- G_num + rnorm(n)
  X <- apply(matrix(runif(n * p), n, p), 2, function(x) x + G_num)
  list(X = X, G = G, Y = Y)
}

test_that("permutation encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("permutation", X, G)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})

test_that("multi_permutation encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("multi_permutation", X, G, num_permutations = 11)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})


test_that("means encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("means", X, G)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})


test_that("low_rank encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("low_rank", X, G, num_components = 3)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})


test_that("sparse_low_rank encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("sparse_low_rank", X, G, num_components = 3)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})


test_that("mnl encoding", {
  data <- make_data()
  X <- data$X
  G <- data$G
  enc <- make_encoder("mnl", X, G)
  X_enc <- enc(X, G)
  expect_true(has_consistent_encoding(X, G, X_enc))
})
