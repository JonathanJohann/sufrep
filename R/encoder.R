




one_hot_encode <- function(k){
  CM <- data.frame(diag(k))
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

helmert_encode <- function(X,G,k){
  CM <- stats::model.matrix(~G,data=X,contrasts = list(G = "contr.helmert"))
  CM <- CM[,-1] #Drop the intercept term
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

deviation_encode <- function(X,G,k){
  CM <- stats::model.matrix(~G,data=X,contrasts = list(G = "contr.sum"))
  CM <- CM[,-1] #Drop the intercept term
  colnames(CM) <- paste("E", 1:k, sep = "")
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
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

difference_encode <- function(k){
  CM <- matrix(0, nrow = k, ncol = k)
  CM <- matrix(-1/(col(CM) + 1), nrow = k, ncol = k)
  CM[lower.tri(CM)] <- 0
  CM <- CM[, 1:(k - 1)]
  CM[row(CM) == (col(CM) + 1)] <- -apply(CM, 2, sum)
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

simple_effect_encode <- function(k){
  CM <- matrix(-1/k, nrow = k, ncol = k)
  CM <- CM + diag(k)
  CM <- CM[, 1:(k - 1)]
  colnames(CM) <- paste("E", 1:k, sep = "")
  return(CM)
}

fisher_encode <- function(X,G,Y){
  
  CM <- aggregate(X[,Y],list(X[,G]),mean)
  colnames(CM) <- c("label","X1")
  CM$X1 <- rank(CM$X1)
  ordering <- data.frame(unique(X[,G]))
  colnames(ordering) <- "ORD"
  CM <- CM[order(ordering$ORD),]
  #colnames(CM) <- paste("E1", sep = "")
  return(CM$X1)
}

means_encode <- function(X,G){
  
  CM <- aggregate(X,list(X[,G]),mean)
  colnames(CM) <- paste("E",1:dim(CM)[2],sep="")
  ordering <- 
  return(CM)
}

low_rank_encode <- function(X,G,num_components){
  
  CM <- means_encode(X,G)
  decomp <- tryCatch({svd(CM)},
                     error=function(e){
                       return(svd(CM[,colSums(!is.finite(CM))==0]))
                     })
  np <- min(num_components,dim(decomp$u)[2])
  CM <- data.frame(decomp$u[,1:np])
  colnames(CM) <- paste("E",1:np,sep="")
}

sparse_low_rank_encode <- function(X,G,num_components){
  
  CM <- means_encode(X,G)
  decomp <- tryCatch({sparsepca::spca(CM)},
                     error=function(e){
                       return(sparsepca::spca(CM[,colSums(!is.finite(CM))==0]))
                     })
  np <- min(num_components,dim(decomp$loadings)[2])
  U <- decomp$loadings[,1:np]
  CM <- tryCatch({CM %*% U},
                 error=function(e) {
                   return(CM[,colSums(!is.finite(M))==0] %*% U)
                 })
  colnames(CM) <- paste("E",1:np,sep="")
  return(CM)
}

mnl_encode <- function(X,G){
  
  categorical_names <- which(colnames(X) %in% c(G))
  Y_mnl <- as.matrix(X[,categorical_names])
  X_mnl <- as.matrix(X[,-categorical_names])
  
}

#' @export 
encoder <- function(X, G = "A", Y=NULL, num_components = NULL, num_folds=4, method = "one_hot",...) {
    
    id <- data.frame(unique(X[, G]))
    k <- dim(id)[1]
    
    CM <- switch(method,
                 one_hot = one_hot_encode(k),
                 helmert = helmert_encode(X,G,k),
                 deviation = deviation_encode(X,G,k),
                 repeated_effect = repeated_effect_encode(k),
                 difference = difference_encode(k),
                 simple_effect = simple_effect_encode(k),
                 fisher = fisher_encode(X,G,Y),
                 means = means_encode(X,G),
                 low_rank = low_rank_encode(X,G,num_components),
                 sparse_low_rank = sparse_low_rank_encode(X,G,num_components),
                 MNL = mnl_encode(X,G))
  
    map <- cbind(id, CM)
    
    
    f <- function(X, drop = TRUE) {
        X <- merge(x = X, y = map, by = G)
        X <- X[,-which(colnames(X) %in% c(G))]
    }
    
    return(f)
}
