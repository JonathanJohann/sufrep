library(sufrep)
library(grf)
library(caret)
library(xgboost)
type <- "ames"
filename <- paste0(method, "_", time_seed(), ".csv", collapse = "")

num_comp <- c(5,10,15)

#==========================
df <- sufrep::ames
data <- list(x = df[,-c(13,78)],
             y = df$SalePrice,
             g = factor(df$Neighborhood))

#==========================
df <- sufrep::pakistan
data <- list(x = df[,-c(1,12)],
             y=df[,1],
             g=factor(df[,12]))


#==========================
df <- sufrep::kingcounty
data <- list(x= df[,-c(1,15)],
             y=df[,1],
             g=factor(df[,15]))

set.seed(123123)
start <- Sys.time()
methods <-c("mnl")#multi_permutation","means","low_rank","sparse_low_rank","mnl")
results <- c()
for (iz in 1:length(methods)) {
  print(i)
  method <- methods[iz]
  try({
    folds1 <- createFolds(factor(data$g),k = 4,returnTrain = T)
    #data <- create_data(n, p, k, ngl = ngl, pl = pl, type = type)
    mses <- c()
    for(i in 1:4){
      train <- list(x = apply(data$x[folds1[[i]],],
                              2,as.numeric),
                    g = data$g[folds1[[i]]],
                    y = data$y[folds1[[i]]])
      
      test <- list(x = apply(data$x[-folds1[[i]],],
                             2,as.numeric),
                   g = data$g[-folds1[[i]]],
                   y = data$y[-folds1[[i]]])
      
      if(method %in% c("fisher")){
        enc_method <- make_encoder(method, X = train$x, G = train$g,Y=train$y)
        x_enc <- enc_method(train$x, train$g)
        x_test_enc <- enc_method(test$x, test$g)
      }
      else if(method %in% c("multi_permutation")){
        enc_method <- make_encoder(method, X = train$x, G = train$g,num_permutations = p)
        x_enc <- enc_method(train$x, train$g)
        x_test_enc <- enc_method(test$x, test$g)
      }
      else if(method %in% c("low_rank","sparse_low_rank")){
        folds <- createFolds(factor(train$g),k = 3,returnTrain = T)
        cv_mses <- c()
        for(ii in 1:length(num_comp)){
          cmp <- num_comp[ii]
          fold_mses <- c()
          for(iii in 1:3){
            #make encoder function based on training subset
            enc_method_CV <- make_encoder(method, X = train$x[folds[[iii]],], G = train$g[folds[[iii]]],num_components = cmp)
            
            #CV train & test
            x_enc_CV <- enc_method_CV(train$x[folds[[iii]],], train$g[folds[[iii]]])
            
            x_test_enc_CV <- enc_method_CV(train$x[-folds[[iii]],], train$g[-folds[[iii]]])
            
            forest_enc_CV <- regression_forest(x_enc_CV, train$y[folds[[iii]]])
            mse_enc_CV <- mean((predict(forest_enc_CV,x_test_enc_CV)$predictions - train$y[-folds[[iii]]])^2, na.rm = TRUE)
            fold_mses <- c(fold_mses,mse_enc_CV)
            print(paste0("CV ... ",iii))
          }
          cv_mses <- c(cv_mses,mean(fold_mses))
          print(paste0("Checking - ",ii))
        }
        enc_method <- make_encoder(method, X = train$x, G = train$g,num_components = which.min(cv_mses))
        x_enc <- enc_method(train$x, train$g)
        x_test_enc <- enc_method(test$x, test$g)
      } else if(method %in% c("mnl")){
        enc_method <- tryCatch({make_encoder(method, X = train$x, G = train$g)},
                               error=function(e){make_encoder(method, X = train$x, G = train$g)})
        x_enc <- enc_method(train$x, train$g)
        x_test_enc <- enc_method(test$x, test$g)
      }
      else{
        enc_method <- make_encoder(method, X = train$x, G = train$g)
        x_enc <- enc_method(train$x, train$g)
        x_test_enc <- enc_method(test$x, test$g)
      }
      enc_onehot <- make_encoder("one_hot", X = train$x, G = train$g)
      x_onehot <- enc_onehot(train$x, train$g)
      x_test_onehot <- enc_onehot(test$x,test$g)
      
      #print(get_xgboost_mse(cbind(x_onehot,data.frame(Y=train$y)),cbind(x_test_onehot,data.frame(Y=test$y))))
      forest_enc <- regression_forest(x_enc, train$y)
      forest_onehot <- regression_forest(x_onehot, train$y)
      
      mse_enc <- mean((predict(forest_enc,x_test_enc)$predictions - test$y)^2, na.rm = TRUE)
      mse_onehot <- mean((predict(forest_onehot,x_test_onehot)$predictions - test$y)^2, na.rm = TRUE)
      
      config <- cbind(method, type, other = NA, mse_enc, mse_onehot)
      write.table(config, file = filename, append = T, col.names = F, sep = ",")
      mses <- c(mses, mse_enc/mse_onehot)
    }
    
    results <- rbind(results,mean(mses))
  })
  end <- Sys.time()
  print(end-start)
}
print(results)
