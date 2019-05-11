




one_hot_encode <- function(k){
  CM <- data.frame(diag(k))
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

helmert_encode <- function(k){
  CM <- stats::contr.helmert(k)
  colnames(CM) <- paste("E", 1:(k-1), sep = "")

  return(CM)
}

deviation_encode <- function(k){
  CM <- stats::contr.sum(k)
  colnames(CM) <- paste("E", 1:(k-1), sep = "")
  return(CM)
}

repeated_effect_encode <- function(k){
  TH <- matrix(0, nrow = k, ncol = k)
  TH <- matrix((k - col(TH))/k, nrow = k, ncol = k)
  TH[lower.tri(TH)] <- 0
  BH <- matrix(-col(TH)/k, nrow = k, ncol = k)
  BH[upper.tri(BH)] <- 0
  diag(BH) <- 0
  CM <- TH + BH
  CM <- CM[, 1:(k - 1)]
  colnames(CM) <- paste("E", 1:dim(CM)[2], sep = "")
  return(CM)
}

difference_encode <- function(k){
  CM <- matrix(0, nrow = k, ncol = k)
  CM <- matrix(-1/(col(CM) + 1), nrow = k, ncol = k)
  CM[lower.tri(CM)] <- 0
  CM <- CM[, 1:(k - 1)]
  CM[row(CM) == (col(CM) + 1)] <- -apply(CM, 2, sum)
  colnames(CM) <- paste("E", 1:(k-1), sep = "")
  return(CM)
}

simple_effect_encode <- function(k){
  CM <- matrix(-1/k, nrow = k, ncol = k)
  CM <- CM + diag(k)
  CM <- CM[, 1:(k - 1)]
  colnames(CM) <- paste("E", 1:(k-1), sep = "")
  return(CM)
}

fisher_encode <- function(X,G,Y){

  CM <- aggregate(X[,Y],list(X[,G]),mean)
  colnames(CM) <- c("label","X1")
  CM$X1 <- rank(CM$X1)
  ordering <- data.frame(unique(X[,G]))
  colnames(ordering) <- "ORD"
  CM <- CM[order(ordering$ORD),]
  CM <- data.frame(CM$X1)
  colnames(CM) <- paste("E1", sep = "")
  return(CM)
}

means_encode <- function(X,G){

  CM <- aggregate(X,list(X[,G]),mean)
  colnames(CM) <- paste("E",1:dim(CM)[2],sep="")
  ordering <-
  return(CM)
}

low_rank_encode <- function(X,G,Y=NULL,num_components=NULL,folds=3,cv_vals=c(5,10,15),cross_validate=TRUE,model="regression_forest"){

  if((cross_validate==TRUE)|(is.null(num_components))){
    mses <- c()
    randomized_df <- X[sample(nrow(X)), ]
    rownames(randomized_df) <- NULL
    fold_cat <- category_stratify(randomized_df[,G],num_folds=folds)
    for(i in 1:length(cv_vals)){
      mse <- c()
      print(paste("--CV FOLD--",i,sep=""))
      for(j in 1:folds){
        testIndexes <- fold_cat[[j]]
        testData <- randomized_df[testIndexes, ]
        trainData <- randomized_df[-testIndexes, ]
        id <- data.frame(unique(trainData[, G]))
        colnames(id) <- G
        k <- dim(id)[1]
        train.X2 <- trainData[,-which(colnames(trainData) %in% c(Y))]
        CM <- low_rank_encode(train.X2,G=G,num_components=cv_vals[i],cross_validate=FALSE,model=model)

        map <- data.frame(cbind(id, CM))


        enc <- function(X) {

          X_enc <- map[X[,G],]
          rownames(X_enc) <- NULL
          X <- cbind(X,X_enc)
          X <- X[,-which(colnames(X) %in% c(G))]
          return(X)
        }
        trainData <- enc(trainData)
        testData <- enc(testData)
        if(model=="regression_forest"){
          print(dim(trainData))
          print(dim(testData))
          mse <- c(mse,get_forest_mse(trainData,testData))
        }
        else{
          mse <- c(mse,get_xgboost_mse(trainData,testData))
        }
      }
      mses <- c(mses,mean(mse))
    }
    mx <- which(mses==min(mses))[1]
    num_components <- cv_vals[mx]
  }
  CM <- means_encode(X,G)
  decomp <- tryCatch({svd(CM)},
                     error=function(e){
                       return(svd(CM[,colSums(!is.finite(CM))==0]))
                     })

  np <- min(num_components,dim(decomp$u)[2])
  CM <- data.frame(decomp$u[,1:np])
  colnames(CM) <- paste("E",1:np,sep="")
  return(CM)
}

sparse_low_rank_encode <- function(X,G,Y=NULL,num_components=NULL,folds=3,cv_vals=c(5,10,15),cross_validate=TRUE,model="regression_forest"){

  if((cross_validate==TRUE)|(is.null(num_components))){
    mses <- c()
    randomized_df <- X[sample(nrow(X)), ]
    rownames(randomized_df) <- NULL
    fold_cat <- category_stratify(randomized_df[,G],num_folds=folds)
    for(i in 1:length(cv_vals)){
      print(paste("--CV FOLD--",i,sep=""))
      mse <- c()
      for(j in 1:folds){
        testIndexes <- fold_cat[[j]]
        testData <- randomized_df[testIndexes, ]
        trainData <- randomized_df[-testIndexes, ]
        id <- data.frame(unique(trainData[, G]))
        colnames(id) <- G
        k <- dim(id)[1]
        train.X2 <- trainData[,-which(colnames(trainData) %in% c(Y))]
        CM <- sparse_low_rank_encode(train.X2,G=G,num_components=cv_vals[i],cross_validate=FALSE,model=model)

        map <- data.frame(cbind(id, CM))


        enc <- function(X) {

          X_enc <- map[X[,G],]
          rownames(X_enc) <- NULL
          X <- cbind(X,X_enc)
          X <- X[,-which(colnames(X) %in% c(G))]
          return(X)
        }
        trainData <- enc(trainData)
        testData <- enc(testData)
        if(model=="regression_forest"){
          mse <- c(mse,get_forest_mse(trainData,testData))
        }
        else{
          mse <- c(mse,get_xgboost_mse(trainData,testData))
        }
      }
      mses <- c(mses,mean(mse))
    }
    mx <- which(mses==min(mses))[1]
    num_components <- cv_vals[mx]
  }
  CM <- means_encode(X,G)
  decomp <- tryCatch({sparsepca::spca(CM,verbose=FALSE)},
                     error=function(e){
                       return(sparsepca::spca(CM[,colSums(!is.finite(CM))==0]))
                     })

  np <- min(num_components,dim(decomp$loadings)[2])
  U <- decomp$loadings[,1:np]
  CM <- as.matrix(CM)
  CM <- tryCatch({CM %*% U},
                 error=function(e) {
                   return(CM[,colSums(!is.finite(CM))==0] %*% U)
                 })
  CM <- data.frame(CM)
  colnames(CM) <- paste("E",1:np,sep="")
  return(CM)
}


permutation_encode <- function(k,num_permutations=1){
  set.seed(time_seed())
  CM <- data.frame(E1 = sample(k,size=k,replace=FALSE))
  if(num_permutations>1){
    mm <- time_seed()
    for(i in 2:num_permutations){
      set.seed(time_seed()*i %% mm)
      CM <- cbind(CM,data.frame(sample(k,size=k,replace=FALSE)))
    }
  }
  colnames(CM) <- paste("E",1:num_permutations,sep="")
  return(CM)
}

mnl_encode <- function(X,G,k,folds=3){

  categorical_names <- which(colnames(X) %in% c(G))
  train.Y <- as.matrix(X[,categorical_names])
  train.X <- as.matrix(X[,-categorical_names])
  fold_cat <- category_stratify(train.Y,num_folds=folds)
  labels <- matrix(0,ncol=length(train.Y))
  for(j in 1:folds){
    labels[,fold_cat[[j]]] <- j
  }
  remainder <- which(labels==0)
  labels[,remainder] <- sample(c(1:folds),size=length(remainder),replace=TRUE)

  fit <- cv.glmnet(x=train.X,y=train.Y,foldid=labels,family="multinomial")
  coef_vals <- coef(fit,s="lambda.1se")
  CM <- data.frame(matrix(0,ncol=length(colnames(train.X))))
  colnames(CM) <- colnames(train.X)
  all_cols <- colnames(CM)
  for(i in 1:k){
    tmp <- tryCatch({as.data.frame(as.matrix(t(coef_vals[[i]])))},error=function(e){return(CM[1,])})
    tmp_cols <- colnames(tmp)
    needed <- all_cols[-which(tmp_cols %in% all_cols)]
    if(length(needed)>0){
      for(l in 1:length(needed)){
          tmp[,needed[l]] <- 0
      }
    }
    CM <- rbind(CM,tmp[,all_cols])
  }

  CM <- CM[-1,]
  rownames(CM) <- NULL
  colnames(CM) <- paste("E",1:dim(CM)[2],sep="")
  return(CM)
}

#' @export
encoder <- function(X, G = "A", Y=NULL, num_components = NULL, num_folds=4, method = "one_hot",...) {


    id <- data.frame(unique(X[, G]))
    colnames(id) <- G
    k <- dim(id)[1]

    CM <- switch(method,
                 one_hot = one_hot_encode(k),
                 helmert = helmert_encode(k),
                 deviation = deviation_encode(k),
                 repeated_effect = repeated_effect_encode(k),
                 difference = difference_encode(k),
                 simple_effect = simple_effect_encode(k),
                 fisher = fisher_encode(X,G,Y),
                 means = means_encode(X,G),
                 low_rank = low_rank_encode(X,G,Y=Y,cross_validate=TRUE),
                 sparse_low_rank = sparse_low_rank_encode(X,G,Y=Y,cross_validate=TRUE),
                 permutation = permutation_encode(k,num_permutations=1),
                 multi_permutation = permutation_encode(k,num_permutations = (dim(X)[2]-1)),
                 MNL = mnl_encode(X,G,k))

    map <- data.frame(cbind(id, CM))


    f <- function(X) {

        X_enc <- map[X[,G],]
        rownames(X_enc) <- NULL
        X <- cbind(X,X_enc)
        X <- X[,-which(colnames(X) %in% c(G))]
        return(X)
    }

    return(f)
}









