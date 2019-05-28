

# Simulation Steps

# Step 1: Choose k, Nl, Ng, Pl

# Step 2: Sample k mu's where mu in Unif([-1,0,1])^p

# Step 3: Sample Nl points from N(u_{l=i},I_{pxp}) \forall i in {1,...,k}

# Step 4: Sample beta in Unif([-1,0,1])^p, alpha_{l=i} in Unif([-10,...,10]),
#         ||alpha|| = 1, and ||beta||=1

# Step 5: Permute {1,...,|G|} where |G| = k * Ng and
#         G_{l=i} = {gj in G | (l-1)*Ng + 1 < j <= l*Ng}
#         where i in {1,...,k}

# Step 6: Sample G ~ { Unif(Gl),   prob = Pl,
#                      Unif(G\Gl), prob = 1 - Pl}

# Step 7: Yi = alpha_{l=l'} + X_{l=l'} \beta + eps, eps in N(0,2)


k = 10
Nl = 1000
Ng = 50
Pl = 0.9
p = 20
sample_mu <- function(k,p){
  mu = matrix(sample(c(-1,0,1),
                     size=k*p,
                     replace=TRUE),
              nrow=k,ncol=p)
  return(mu)
}

sample_points <- function(mu,nl){
  k <- dim(mu)[1]
  p <- dim(mu)[2]
  X <- c()
  for(i in 1:k){
    Xl <- MASS::mvrnorm(nl,mu=mu[i,],Sigma=diag(p))
    X <- rbind(X,cbind(Xl,rep(i,nl)))
  }
  colnames(X) <- c(paste("X",1:p,sep=""),"L")
  return(X)
}

sample_beta_global <- function(p){
  beta = matrix(sample(c(-1,0,1),
                       size=p,
                       replace=TRUE),
                nrow=1,ncol=p)
  beta = beta/sum(beta^2)
  return(beta)
}

sample_beta_latent <- function(k,p){
  betas = c()
  for(i in 1:k){
    beta = matrix(sample(c(-1,0,1),
                         size=p,
                         replace=TRUE),
                  nrow=1,ncol=p)
    beta = beta/sum(beta^2)
    betas <- rbind(betas,beta)
  }
  return(betas)
}

sample_alpha <- function(k){
  alpha = matrix(stats::rnorm(k),
                 nrow=k,ncol=1)
  alpha = alpha/sum(alpha^2)
  return(alpha)
}

permute_categories <- function(k,ng){
  map <- matrix(sample(c(1:(k*ng)),
                       replace=FALSE,
                       size=k*ng),
                nrow=k,ncol=ng)
  return(map)
}

sample_observed <- function(X,cat,pl){
  L = dim(X)[2]
  n = dim(X)[1]
  t = stats::rbinom(n=n,size=1,prob=pl)
  G = matrix(0,nrow=n,ncol=1)
  for(i in 1:n){
    G[i,] <- ifelse(t[i]==1,
                     sample(cat[X[i,L],],size=1),
                     sample(cat[-X[i,L],],size=1))
  }
  colnames(G) <- "G"
  return(G)
}

create_response <- function(X,alpha,beta,type="global"){
  n = dim(X)[1]
  p = dim(X)[2]
  k = max(dim(alpha))
  eps = rnorm(n,sd=sqrt(2))
  if(type=="global"){
    Y <- X[,-p] %*% t(beta) + alpha[X[,p]] + eps
  }
  else if(type=="latent"){
    Y <- matrix(0,nrow=n,ncol=1)
    for(i in 1:k){
      l <- which(X[,"L"]==i)
      bl <- as.matrix(beta[i,],nrow=p,ncol=1)
      Y[l,] <- X[l,-p] %*% bl  + alpha[i,] + eps[l]
    }
  }
  colnames(Y) <- "Y"
  return(Y)
}







simulation <- function(p,k,nl,ng,pl,type="global",seed=NULL){
  if(is.null(seed)){seed = time_seed()}
  set.seed(seed)
  u = sample_mu(k,p)
  X = sample_points(u,nl)
  a = sample_alpha(k)
  if(type=="global"){
    b = sample_beta_global(p)
  } else if(type=="latent"){
    b = sample_beta_latent(k,p)
  }
  g = permute_categories(k,ng)
  G = sample_observed(X=X,cat=g,pl=pl)
  Y = create_response(X=X,alpha=a,beta=b,type=type)
  output = as.data.frame(cbind(X[,-(p+1)],G,Y))
  return(output)
}


