



library(tidyverse)
library(grf)
library(here)
library(levels)
library(caret)
library(glmnet)
library(sufrep)
direc <- here::here("Datasets", "pakistan.csv")

pakistan <- readr::read_csv("C:\\Users\\jonat\\Desktop\\levels\\datasets\\pakistan.csv")
pakistan <- pakistan[-which(pakistan$City %in% c('Tor Ghar')),]
pakistan <- pakistan %>% na.omit()

relevant_columns <- c('Education score','Toilet','Province','Population','School infrastructure score','Total number of schools','Primary Schools with single teacher','Primary Schools with single classroom',
                      'Pakistan Economic Growth','Number of secondary schools','Electricity','No Facility',
                      'City','Global Terrorism Index - Pakistan','Complete Primary Schools','Building condition satisfactory',
                      'Drone attacks in Pakistan','Drinking water','Boundary wall','Bomb Blasts Occurred','% Complete Primary Schools','% Boys Enrolled')

pakistan <- pakistan[,relevant_columns]
percent_boys <- pakistan %>% pull(names(pakistan)[1])
pakistan[,1] <- sapply(gsub("%","",percent_boys),as.numeric)


pakistan_columns <- pakistan %>% Filter(f = is.character) %>% names
pakistan <- pakistan %>% mutate_at(pakistan_columns,funs(as.factor)) %>%
  mutate_at(pakistan_columns,funs(as.integer))


pakistan <- pakistan %>% mutate_at("City",funs(as.factor))

filename <- paste("xgb_pakistan_45seeds_10folds",time_seed(),".rds",sep="")
folds = 4
seeds = 1
y_var = "Education.score"
target_categorical = "City"
dataset = data.frame(pakistan)
dataset$Y = dataset[,y_var]
dataset = dataset[,-which(colnames(dataset) %in% c(y_var))]
dataset[,target_categorical] <- as.numeric(as.character(dataset[,target_categorical]))

start_time = Sys.time()
eval.df = evaluate_method(dataset,response="Y",categorical = target_categorical,k=folds,seeds=seeds)
saveRDS(eval.df,file=filename)
end_time = Sys.time()
print(end_time - start_time)
