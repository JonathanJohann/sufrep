






#' @export 
encoder <- function(X, G = "A", Y=NULL, num_components = NULL, num_folds=4, method = "one hot",...) {
    
    id <- data.frame(unique(X[, G]))
    k <- dim(id)[1]
    
    if (method == "one hot") {
        
        CM <- data.frame(diag(k))
        colnames(CM) <- paste("E", 1:k, sep = "")
        
    } else if (method == "helmert") {
        
        CM <- matrix(0, nrow = k, ncol = k)
        CM <- matrix(-1/(k - col(CM) + 1), nrow = k, ncol = k)
        CM[upper.tri(CM)] <- 0
        CM <- CM + diag(k)
        CM <- CM[, 1:(k - 1)]
        colnames(CM) <- paste("E", 1:k, sep = "")
        
    } else if (method == "deviation") {
      
        CM <- diag(k)
        CM[k, ] <- CM[k, ] - 1
        CM <- CM[, 1:(k - 1)]
        colnames(CM) <- paste("E", 1:k, sep = "")
    
    } else if (method == "repeated effect") {
        
        TH <- matrix(0, nrow = k, ncol = k)
        TH <- matrix((k - col(TH))/k, nrow = k, ncol = k)
        TH[lower.tri(TH)] <- 0
        BH <- matrix(-col(TH)/k, nrow = k, ncol = k)
        BH[upper.tri(BH)] <- 0
        diag(BH) <- 0
        CM <- TH + BH
        CM <- CM[, 1:(k - 1)]
        colnames(CM) <- paste("E", 1:k, sep = "")
        
    } else if (method == "difference") {
    
        CM <- matrix(0, nrow = k, ncol = k)
        CM <- matrix(-1/(col(CM) + 1), nrow = k, ncol = k)
        CM[lower.tri(CM)] <- 0
        CM <- CM[, 1:(k - 1)]
        CM[row(CM) == (col(CM) + 1)] <- -apply(CM, 2, sum)
        colnames(CM) <- paste("E", 1:k, sep = "")
        
    } else if (method == "simple effect") {
        
        CM <- matrix(-1/k, nrow = k, ncol = k)
        CM <- CM + diag(k)
        CM <- CM[, 1:(k - 1)]
        colnames(CM) <- paste("E", 1:k, sep = "")
        
    } else if (method == "fisher") {
      
        CM <- X %>% 
                dplyr::group_by(!!rlang::sym(G)) %>% 
                dplyr::summarise_at(Y,mean)
        CM <- data.frame(CM[order(id),] %>% 
                           dplyr::select(-!!rlang::sym(G)))
        colnames(CM) <- paste("E1", sep = "")
      
    } else if (method == "means") {
      
      CM <- X %>%
              dplyr::group_by(!!rlang::sym(G)) %>%
              dplyr::summarise_all(mean)
      CM <- data.frame(CM[order(id),] %>%
                         dplyr::select(-!!rlang::sym(G)))
      colnames(CM) <- paste("E",1:dim(CM)[2],sep="")
      
    } else if (method == "low rank") {
      
      CM <- X %>%
        dplyr::group_by(!!rlang::sym(G)) %>%
        dplyr::summarise_all(mean)
      CM <- as.matrix(CM[order(id),] %>%
                         dplyr::select(-!!rlang::sym(G)))
      decomp <- tryCatch({svd(CM)},
                         error=function(e){
                           return(svd(CM[,colSums(!is.finite(CM))==0]))
                           })
      np <- min(num_components,dim(decomp$u)[2])
      CM <- data.frame(decomp$u[,1:np])
      colnames(CM) <- paste("E",1:np,sep="")
      
    } else if (method == "sparse low rank") {
      
      CM <- X %>%
        dplyr::group_by(!!rlang::sym(G)) %>%
        dplyr::summarise_all(mean)
      CM <- as.matrix(CM[order(id),] %>%
                        dplyr::select(-!!rlang::sym(G)))
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
      
    } else if (method == "MNL") {
      
      Y_mnl <- as.matrix(X %>% dplyr::pull(!!rlang::sym(G)))
      X_mnl <- as.matrix(X %>% dplyr::select(-!!rlang::sym(G)))
      
      return("Not done.")
    }
    map <- cbind(id, CM)
    
    
    f <- function(X, drop = TRUE) {
        X <- left_join(x = X, y = map, by = G)
        if (drop) {
            X %>% dplyr::select(-G)
        }
    }
    
    return(f)
}

