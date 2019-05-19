
rm(list=ls())

library(tidyverse)
library(xgboost)
library(glmnet)
library(sufrep)
library(sparsepca)

source("dgp.R")
source("utils.R")

config <- expand.grid(
  setup="interactions",
  num_groups=c(2,10),
  n=c(5000),
  p=c(20),
  rhos=c(0.25),
  num_trees=c(2000),
  noise_ratio=c(0.5),
  iterations=1,
  signal_fraction=c(0.5),
  comp_percent=c(0.9),
  num_levels=c(500,100),
  num_interactions=c(1)
)

shuffled.indices <- sample(dim(config)[1]-1,
                           size=dim(config)[1]-1,
                           replace=FALSE)

config <- config[c(1,(shuffled.indices+1)),]

encodings <- c("MNL","one_hot","multi_permutation","permutation",
               "difference","helmert","fisher","deviation",
               "low_rank","sparse_low_rank","repeated_effect",
               "means","simple_effect")

walk(seq(dim(config)[1]), function(i) {

  cfg <- config[i,]

  output <- rerun(.n=cfg$iterations, {

    seed <- time_seed()
    sim_data <- make_dataset(setup=as.character(cfg$setup),
                             n=cfg$n,p=cfg$p,
                             rho=cfg$rhos,
                             size_l=cfg$num_groups,
                             size_a=cfg$num_levels,
                             lambda=cfg$signal_fraction,
                             snr=cfg$noise_ratio,
                             comp_percent=cfg$comp_percent,
                             seed=seed,
                             num_interactions = cfg$num_interactions,
                             test=TRUE)

    mses <- matrix(0,nrow=1,ncol=13)
    print("Starting...")
    train = data.frame(sim_data$TRAIN %>% dplyr::mutate_at("A",fix_factors))
    test = data.frame(sim_data$TEST %>% dplyr::mutate_at("A",fix_factors))
    test = test[-which(test$A %in% c(which(table(train$A)<4))),]
    train = train[which(train$A==which(table(train$A)<4)),]
    for(ii in 1:13){
      if(encodings[ii]=="fisher"){
        mses[ii] <- evaluate(method=encodings[ii],train=train,test=test,categorical = "A",response="Y",Y="Y")
      } else if(encodings[ii] %in% c("low_rank","sparse_low_rank")) {
        mses[ii] <- evaluate(method=encodings[ii],train=train,test=test,categorical = "A",response="Y",num_components=3)
      } else {
        mses[ii] <- evaluate(method=encodings[ii],train=train,test=test,categorical = "A",response="Y")
      }
      print(ii)
    }

    colnames(mses) <- encodings
    result <- cbind(cfg,mses)
    result[,"seed"] <- seed
    as_tibble(result)
  }) %>% bind_rows()

  # Write to csv
  filename <- paste(cfg$setup,"_n", cfg$n,
                    "_p",cfg$p,
                    "_rhos",cfg$rhos,
                    "_snr", cfg$noise_ratio,
                    "_comp", cfg$comp_percent,
                    "_ngroups", cfg$num_groups,
                    "_nlevels", cfg$num_levels,
                    "_", time_seed(),
                    ".csv",
                    sep="")
  write_csv(output, filename)
})

