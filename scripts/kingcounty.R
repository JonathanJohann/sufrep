library(tidyverse)
library(grf)
library(levels)
library(here)
library(glmnet)
library(sufrep)

direc <- here::here("Datasets")
king_file <- here::here("Datasets", "kc_house_data.csv")

kingcounty <- readr::read_csv(king_file)
kingcounty <- kingcounty %>%
  na.omit() %>%
  dplyr::select(-date) %>%
  dplyr::select(-id)
kingcounty <- kingcounty %>% mutate_at("zipcode",funs(as.factor))



filename <- paste("NEW_RF_king_10folds",time_seed(),".rds",sep="")
folds = 4
seeds = 1
y_var = "price"
target_categorical = "zipcode"
dataset = data.frame(kingcounty)
dataset$Y = dataset[,y_var]
dataset = dataset[,-which(colnames(dataset) %in% c(y_var))]
dataset[,target_categorical] <- as.numeric(as.character(dataset[,target_categorical]))

start_time = Sys.time()
eval.df = evaluate_method(dataset,response="Y",categorical = target_categorical,k=folds,seeds=seeds)
saveRDS(eval.df,file=filename)
end_time = Sys.time()
print(end_time - start_time)

