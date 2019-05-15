


evaluate <- function(train,test,
                     categorical,response,
                     method="one_hot",
                     model = "regression_forest",...){
  
  remove_response <- which(colnames(train) %in% c(response))
  
  if(method %in% c("fisher","sparse_low_rank","low_rank")){
    train.X <- train
  } else {
    train.X <- train[,-remove_response]
  }

  
  enc <- encoder(X=train.X,G=categorical,method=method,...)
  train <- enc(train)
  test <- enc(test)
  
  return(ifelse(model=="regression_forest",
                get_forest_mse(train,test),
                get_xgboost_mse(train,test)))
}


evaluate_method <- function(df,categorical,response,k=10,seeds=1,model="regression_forest",stratify=TRUE){
  output <- c()
  seed_val <- time_seed()
  set.seed(seed_val)
  randomized_df <- df[sample(nrow(df)), ]
  methods <- c("one_hot","multi_permutation","means","low_rank",
               "sparse_low_rank","MNL","permutation","difference",
               "deviation","repeated_effect","helmert","fisher",
              "simple_effect")
  for (i in 1:seeds) {
    
    total_indices <- 1:nrow(df)
    
    categ <- c(randomized_df %>% dplyr::mutate_at(categorical,as.character) 
               %>% dplyr::mutate_at(categorical,as.numeric) 
               %>% dplyr::select(categorical))
    
    fold_cat <- category_stratify(categ[[1]],num_folds=k)
    
    for (j in 1:k) {
      testIndexes <- fold_cat[[j]]
      testData <- randomized_df[testIndexes, ]
      trainData <- randomized_df[-testIndexes, ]
      mses <- c()
      for (q in 1:13){
        mse <- evaluate(method=methods[q],train=trainData,test=testData,
                        categorical=categorical,response=response,
                        model=model,Y=response)
        print(mse)
        mses <- c(mses,mse)
        print(methods[q])
      }
      
      new_row <- c(nrow(df), model, seed_val, j, mses)
      output <- rbind(output, new_row)
      print(paste("Done with -- ",j,sep=""))
    }
    colnames(output) <- c("file", "model", "seed", "fold", "one_hot", "multi_perm", 
                          "add_means", "add_svd", "add_spca", "add_pax_weight", "perm", 
                          "difference", "deviation", "repeated", "helmert", "fisher","simple_effect")  
    saveRDS(output, file = paste("Evaluation_", k, "_seed", i, "_", seed_val, ".rds", sep = ""))
  }
  return(data.frame(output))
  
}






get_noise_scale <- function(noiseless, snr) {
  sqrt(var(noiseless)/snr)
}



time_seed <- function() {
  as.integer((as.numeric(Sys.time()) * 1e+07)%%1e+07)
}




scale_vals <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


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



get_regression_forest_prediction <- function(train_data, test_data, ...) {
  train_X <- train_data %>% dplyr::select(-Y)
  forest <- grf::regression_forest(X = train_X, Y = train_data$Y, ...)
  yhat <- predict(forest, 
                  newdata = test_data %>% dplyr::select(-Y), 
                  estimate.variance = FALSE)$predictions
  return(yhat)
}




get_xgboost_mse <- function(train,test,...){
  train_Y <- train %>% dplyr::pull(Y)
  train_X <- as.matrix(train %>% dplyr::select(-Y))
  test_Y <- test %>% dplyr::pull(Y)
  test_X <- as.matrix(test %>% dplyr::select(-Y))
  
  xgb_grid_1 = expand.grid(nrounds = c(20,50,100),  
                           max_depth = c(3,6,9,12),
                           colsample_bytree = c(0.5,0.7,0.9),
                           eta = c(0.1,0.3,0.5),
                           gamma=c(0,0.1),
                           min_child_weight = c(1,5,10),
                           subsample = c(0.5,0.75,1.0)
  )
  
  xgb_trcontrol_1 = trainControl(
    method = "cv",
    number = 3,  
    allowParallel = TRUE
  )

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



get_forest_mse <- function(train_data, test_data, ...) {
  yhat <- tryCatch({get_regression_forest_prediction(train_data = train_data, test_data = test_data, ...)},
                   error=function(e){return(rbind(train_data,test_data))})
  print(yhat)
  mse <- mean((test_data$Y - yhat)^2)
  return(mse)
}


fix_factors <- function(x) {
  return(as.numeric(as.character(x)))
}

