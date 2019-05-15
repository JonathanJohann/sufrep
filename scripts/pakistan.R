

rm(list=ls())

library(tidyverse)
library(grf)
library(caret)
library(glmnet)
library(sufrep)

source("dgp.R")
source("utils.R")

model = "regression_forest"
folds = 4
seeds = 1
y_var = "Education.score"
target_categorical = "City"

dataset = sufrep::pakistan
dataset$Y = dataset[,y_var]
dataset = dataset[,-which(colnames(dataset) %in% c(y_var))]

start_time = Sys.time()
eval.df = evaluate_method(dataset,response="Y",categorical = target_categorical,k=folds,seeds=seeds,model=model)
filename <- paste(model,"_PAKISTAN_",time_seed(),".rds",sep="")
saveRDS(eval.df,file=filename)
end_time = Sys.time()
print(end_time - start_time)
