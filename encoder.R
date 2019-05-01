



encoder <- function(X,G="A",method="one hot"){
  
  id <- data.frame(unique(X[,G]))
  k <- dim(id)[1]
  
  if(method=="one hot"){
    
    CM <- data.frame(diag(k))
    colnames(CM) <- paste("E",1:k,sep="")
    
  }
  else if(method=="helmert"){
    
    CM <- matrix(0, nrow = k, ncol = k)
    CM <- matrix(-1/(k - col(CM) + 1), nrow = k, ncol = k)
    CM[upper.tri(CM)] <- 0
    CM <- CM + diag(k)
    CM <- CM[, 1:(k - 1)]
  
  }
  else if(method=="deviation"){
    CM <- diag(k)
    CM[k, ] <- CM[k, ] - 1
    CM <- CM[, 1:(k - 1)]
  }
  else if(method=="repeated effect"){
    
    TH <- matrix(0, nrow = k, ncol = k)
    TH <- matrix((k - col(TH))/k, nrow = k, ncol = k)
    TH[lower.tri(TH)] <- 0
    BH <- matrix(-col(TH)/k, nrow = k, ncol = k)
    BH[upper.tri(BH)] <- 0
    diag(BH) <- 0
    CM <- TH + BH
    CM <- CM[, 1:(k - 1)]
    
  }
  else if(method=="difference"){
    CM <- matrix(0, nrow = k, ncol = k)
    CM <- matrix(-1/(col(CM) + 1), nrow = k, ncol = k)
    CM[lower.tri(CM)] <- 0
    
    CM <- CM[, 1:(k - 1)]
    CM[row(CM) == (col(CM) + 1)] <- -apply(CM, 2, sum)
    
  }
  else if(method=="simple effect"){
    
    CM <- matrix(-1/k, nrow = k, ncol = k)
    CM <- CM + diag(k)
    CM <- CM[, 1:(k - 1)]
  }
  map <- cbind(id,CM)
  
  
  f <- function(X,drop=TRUE){
    X <- left_join(x=X,y=map,by=G)
    if(drop){X %>% dplyr::select(-G)}
  }
  
  return(f)
}



X_temp <- sim_df

z1 <- encoder(X_temp,G="A",method="helmert")


X_temp2 <- z1(X_temp)
