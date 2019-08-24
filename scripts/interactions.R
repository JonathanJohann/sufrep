
rm(list = ls())

library(tidyverse)
library(xgboost)
library(glmnet)
library(sufrep)
library(sparsepca)

source("dgp.R")
source("utils.R")

config <- expand.grid(
  setup = "interactions",
  num_groups = c(2, 10),
  n = c(10000),
  p = c(20),
  rhos = c(0.25),
  num_trees = c(2000),
  noise_ratio = c(0.5),
  iterations = 1,
  signal_fraction = c(0.5),
  comp_percent = c(0.9),
  num_levels = c(100, 500),
  num_interactions = c(1)
)

model <- "regression_forest"

encodings <- c(
  "MNL", "one_hot", "multi_permutation", "permutation",
  "difference", "helmert", "fisher", "deviation",
  "low_rank", "sparse_low_rank", "repeated_effect",
  "means", "simple_effect"
)

seed <- time_seed()

for (i in 1:dim(config)[1]) {
  cfg <- config[i, ]
  sim_data <- make_dataset(
    setup = as.character(cfg$setup),
    n = cfg$n, p = cfg$p,
    rho = cfg$rhos,
    size_l = cfg$num_groups,
    size_a = cfg$num_levels,
    lambda = cfg$signal_fraction,
    snr = cfg$noise_ratio,
    comp_percent = cfg$comp_percent,
    seed = seed,
    num_interactions = cfg$num_interactions,
    test = TRUE
  )

  train <- data.frame(sim_data$TRAIN %>% dplyr::mutate_at("A", fix_factors))
  test <- data.frame(sim_data$TEST %>% dplyr::mutate_at("A", fix_factors))

  result <- evaluate_method(rbind(train, test), categorical = "A", response = "Y", k = 4, model = model)
  result[, "seed"] <- seed
  filename <- paste(model, "_", cfg$setup, "_n", cfg$n,
    "_p", cfg$p,
    "_rhos", cfg$rhos,
    "_snr", cfg$noise_ratio,
    "_comp", cfg$comp_percent,
    "_ngroups", cfg$num_groups,
    "_nlevels", cfg$num_levels,
    "_", time_seed(),
    ".csv",
    sep = ""
  )
  write.csv(result, file = filename)
}
