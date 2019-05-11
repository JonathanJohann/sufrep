

library(tidyverse)
library(grf)
library(caret)
library(glmnet)
library(xgboost)
library(sufrep)


filename <- paste("xgb_ames_45seeds_10fold_",time_seed(),".rds",sep="")
folds = 4
seeds = 1
y_var = "SalePrice"
target_categorical = "Neighborhood"

dataset = data("ames",package="sufrep")
dataset$Y = dataset[,y_var]
dataset = dataset[,-which(colnames(dataset) %in% c(y_var))]


start_time = Sys.time()
eval.df = evaluate_method(df=dataset,response = "Y",categorical=target_categorical,k=folds,seeds=seeds)
saveRDS(eval.df,file=filename)
end_time = Sys.time()
print(end_time - start_time)
