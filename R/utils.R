


#' @export
evaluate <- function(train,test,
                     categorical,response,
                     method="one_hot",
                     model = "regression_forest",...){
  
  remove_response <- which(colnames(train) %in% c(response))
  
  if(method=="fisher"){
    train.X <- train
  } else {
    train.X <- train[,-remove_response]
  }

  
  enc <- encoder(train.X,G=categorical,method=method,...)
  train <- enc(train)
  test <- enc(test)
  
  return(ifelse(model=="regression_forest",
                get_forest_mse(train,test),
                get_xgboost_mse(train,test)))
}






#' @export
get_noise_scale <- function(noiseless, snr) {
  sqrt(var(noiseless)/snr)
}


#' @export
time_seed <- function() {
  as.integer((as.numeric(Sys.time()) * 1e+07)%%1e+07)
}



#' @export
scale_vals <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

#' @export
category_stratify <- function(categories,num_folds=4){
  indices <- 1:length(categories)
  size_per_fold <- as.integer(length(categories)/num_folds)
  folds <- list()
  unique_categories <- unique(categories)
  tmpcat <- categories
  for(j in 1:num_folds){
    vals <- c()
    for(i in 1:length(unique_categories)){
      available <- which(tmpcat==unique_categories[i])
      total_available <- which(categories==unique_categories[i])
      num_select <- as.integer(length(total_available)/num_folds)
      if(length(available)<=num_select){
        added_indices <- available
      }
      else{
        added_indices <- sample(available,num_select,replace=FALSE)
      }
      vals <- c(vals,added_indices)  
    }
    folds[[j]] <- vals
    tmpcat[vals] <- rep(-99999,length(vals))
    indices[vals] <- rep(-99999,length(vals))
  }
  return(folds)
  
}


#' @export
multiple_permutations <- function(df, group, nums = 4, drop = TRUE) {
  output <- df
  observable_groups <- df %>% 
                            dplyr::mutate_at(group,as.character) %>% 
                            dplyr::mutate_at(group,as.numeric) %>% 
                            dplyr::pull(group)
  for (i in 1:nums) {
    set.seed(time_seed())
    tmp <- random_mapping(observable_groups,1)
    tmp <- data.frame(tmp) %>% magrittr::set_colnames(paste("MRP",i,sep=""))
    output <- cbind(output, tmp)
  }
  if (drop == TRUE) {
    output <- output %>% dplyr::select(-group)
  }
  output <- output %>% dplyr::mutate_if(is.factor, as.integer)
  return(output)
}

#' @export
random_mapping <- function(X, levels_per_group) {
  unique_groups <- unique(X)
  g <- length(unique_groups)
  n <- length(X)
  mappings <- matrix(sample(g * levels_per_group, replace = FALSE, 
                            size = (g * levels_per_group)), nrow = g, ncol = levels_per_group)
  mapped_X <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    ind <- sample(1:levels_per_group, size = 1)
    mapped_X[i] <- mappings[which(unique_groups == X[i]), ind]
  }
  return(factor(mapped_X))
}


#' @export
get_regression_forest_prediction <- function(train_data, test_data, ...) {
  train_X <- train_data %>% dplyr::select(-Y)
  forest <- grf::regression_forest(X = train_X, Y = train_data$Y, ...)
  yhat <- predict(forest, 
                  newdata = test_data %>% dplyr::select(-Y), 
                  estimate.variance = FALSE)$predictions
  return(yhat)
}



#' @export
get_xgboost_mse <- function(train,test,...){
  train_Y <- train %>% dplyr::pull(Y)
  train_X <- as.matrix(train %>% dplyr::select(-Y))
  test_Y <- test %>% dplyr::pull(Y)
  test_X <- as.matrix(test %>% dplyr::select(-Y))
  
  xgb_grid_1 = expand.grid(nrounds = c(20,50,100),  # this is n_estimators in the python code above
                           max_depth = c(3,6,9,12),#, 15, 20, 25),
                           colsample_bytree = c(0.5,0.7,0.9),#seq(0.5, 0.9, length.out = 5),
                           ## The values below are default values in the sklearn-api. 
                           eta = c(0.1,0.3,0.5),
                           gamma=c(0,0.1),
                           min_child_weight = c(1,5,10),
                           subsample = c(0.5,0.75,1.0)
  )
  
  # pack the training control parameters
  xgb_trcontrol_1 = trainControl(
    method = "cv",
    number = 3,  
    allowParallel = TRUE
  )
  
  
  # train the model for each parameter combination in the grid, 
  #   using CV to evaluate
  xgb_train_1 = train(x=train_X,
                      y=train_Y,
                      trControl = xgb_trcontrol_1,
                      tuneGrid = xgb_grid_1,
                      method = "xgbTree"
  )
  
  
  predictions <- predict(xgb_train_1,test_X)
  mse <- mean((test_Y-predictions)^2)
  return(mse)
}


#' @export
get_forest_mse <- function(train_data, test_data, ...) {
  yhat <- get_regression_forest_prediction(train_data = train_data, test_data = test_data, ...)
  mse <- mean((test_data$Y - yhat)^2)
  return(mse)
}

#' @export
fix_factors <- function(x) {
  # if you directly convert factors to integers, you won't actually convert the factor '1' to 1 when '1' is not the first observed
  # factor in your data.
  return(as.numeric(as.character(x)))
}

