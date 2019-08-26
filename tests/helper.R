methods <- c("one_hot", "helmert", "deviation",
    "repeated_effect", "difference",
    "simple_effect", "fisher", "means",
    "low_rank", "sparse_low_rank",
    "permutation", "multi_permutation", "mnl")

classical_methods =  c("one_hot", "helmert", "deviation",
                       "repeated_effect", "difference", "simple_effect")

# Creates dummy data.
make_data <- function(n = 200, p = 4, m = 6) {
  cats <- apply(expand.grid(letters, letters), 1, function(s) paste0(s, collapse = ""))[1:m]
  G <- factor(sample(cats, replace = T, size = n))
  G_num <- as.integer(G)
  Y <- G_num + rnorm(n)
  X <- apply(matrix(runif(n * p), n, p), 2, function(x) x + G_num)
  list(X = X, G = G, Y = Y)
}
