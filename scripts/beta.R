
rm(list=ls())

library(tidyverse)
library(xgboost)
library(glmnet)
library(sufrep)
library(sparsepca)

source("new_simulations.R")
source("utils.R")

model = "regression_forest"

config <- expand.grid(
  k=c(5,10),
  Nl=1000,
  Ng=c(50,100),
  Pl=0.9,
  p=20,
  type=c("global","latent")
)


seed <- time_seed()

for(i in 1:8){
  cfg <- config[i,]
  df <- simulation(p=cfg$p,k=cfg$k,
                   nl=cfg$Nl,ng=cfg$Ng,
                   pl=cfg$Pl,type=cfg$type,seed=seed)
  result = evaluate_method(df,categorical="G",response="Y",k=4,model=model)
  result[,"seed"] <- seed
  result[,"model"] <- model
  result <- cbind(result,cfg)
  filename <- paste(model,"_",cfg$type,"type_",cfg$k,"k_",
                    cfg$Nl,"Nl_",cfg$Ng,"Ng_",
                    cfg$Pl,"Pl_",cfg$p,"p_",
                    time_seed(),".csv",sep="")
  write.csv(result,file = filename)
}
