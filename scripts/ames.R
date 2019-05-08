library(tidyverse)
library(grf)
library(levels)
library(here)
library(caret)
library(glmnet)
library(xgboost)
library(sufrep)
ames_file <- "C:\\Users\\jonat\\Desktop\\levels\\datasets\\ames.csv"#here::here("Datasets", "ames.csv")

ames <- readr::read_csv(ames_file)
ames_columns <- ames %>% Filter(f = is.character) %>% names
ames[which(is.na(ames$Alley)),"Alley"] <- 0
ames[which(is.na(ames$LotFrontage)),"LotFrontage"] <- 0
ames <- ames %>% mutate_at(ames_columns,funs(as.factor)) %>%
  mutate_at(ames_columns,funs(as.integer))

ames <- ames %>% mutate_at("Neighborhood",funs(as.factor)) %>% dplyr::select(-PoolQC) %>% dplyr::select(-GarageQual) %>% dplyr::select(-GarageYrBlt)
ames <- data.frame(ames)
ames[is.na(ames)] <- 0
ames <- as.tibble(ames)

ames <- ames[-which(ames$Neighborhood==2),]

seed <- time_seed()
filename <- paste("xgb_ames_45seeds_10fold_",seed,".rds",sep="")
folds = 4
seeds = 1
y_var = "SalePrice"
target_categorical = "Neighborhood"
#dataset = ames


dataset = data.frame(ames)
dataset$Y = dataset[,y_var]
dataset = dataset[,-which(colnames(dataset) %in% c(y_var))]
dataset[,target_categorical] <- as.numeric(as.character(dataset[,target_categorical]))



start_time = Sys.time()
eval.df = evaluate_method(df=dataset,response = "Y",categorical=target_categorical,k=folds,seeds=seeds)
saveRDS(eval.df,file=filename)
end_time = Sys.time()
print(end_time - start_time)
