library(sufrep)
library(grf)

time_seed <- function() {
  as.integer((as.numeric(Sys.time()) * 1e+07) %% 1e+07)
}
encoding_methods <- c("means")
method <- sample(encoding_methods, 1)
n <- 1000
p <- 10
k <- 5 #sample(c(2, 5))
ngl <- 5 #sample(c(5, 10, 25, 50, 100), 1)
pl <- .9
type <- "global"
filename <- paste0(method, "_", time_seed(), ".csv", collapse = "")

for (i in seq(100)) {
  print(i)

  data <- create_data(n, p, k, ngl = ngl, pl = pl, type = type)

  enc_method <- make_encoder(method, X = data$x, G = data$g)
  enc_onehot <- make_encoder("one_hot", X = data$x, G = data$g)

  x_enc <- enc_method(data$x, data$g)
  x_onehot <- enc_onehot(data$x, data$g)

  forest_enc <- regression_forest(x_enc, data$y)
  forest_onehot <- regression_forest(x_onehot, data$y)

  mse_enc <- mean(predict(forest_enc)$debiased.error, na.rm = TRUE)
  mse_onehot <- mean(predict(forest_onehot)$debiased.error, na.rm = TRUE)

  config <- cbind(n, p, k, ngl, pl, method, other = NA, mse_enc, mse_onehot)
  write.table(config, file = filename, append = T, col.names = F, sep = ",")
  print(config)

}
