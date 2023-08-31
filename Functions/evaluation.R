#########################
## model evaluation
##############################
roc = function(est, pop, thres){
  if(dim(pop)[1]!=dim(est)[1]|dim(pop)[2]!=dim(est)[2]){
    stop('the dimension of estimated precision matrix must match the population one')
  }
  if(!isSymmetric(est)){
    stop('estimated precision matrix must be symmetric')
  }
  if(!isSymmetric(pop)){
    stop('true precision matrix must be symmetric')
  }
  p = dim(est)[1]
  est[est<thres] = 0 
  v.est <- c(est[upper.tri(est)])
  v.pop <- c(pop[upper.tri(pop)])

  TP = TN = FP = FN = 0
  
  TP = sum(v.est & v.pop)
  TN = sum(!v.est & !v.pop)
  FP = sum(v.est & !v.pop)
  FN = sum(!v.est & v.pop)
  
  TPR = TP / (TP + FN)
  FPR = FP / (TN + FP)
  precision = TP / (TP + FP)
  recall = TP / (TP + FN)
  
  return(c(TPR,FPR))
}

##############
## roc_curve
###################
roc_curve = function(est, pop, n_grid = 10){
  est_vec = as.vector(est[lower.tri(est)])
  threshold = seq(from = min(est_vec), to = max(est_vec), length.out = n_grid)-10^-5
  sapply(X = threshold, FUN = roc, est = est, pop = pop)
}

################
## roc_curve_bglasso
#################
# roc_curve_bglasso = function(X, p, fun, gamma_n, Omega, ridge=FALSE, rho=1){
#   n = dim(X)[1]
#   m = dim(X)[2]/p
#   if (fun=="cmcopula"){
#     #compute the normal scores
#     a = numeric()
#     for (j in 1:p) {
#       index.j = (((j-1)*m+1):(j*m))
#       a_j = cmcopula(X[,index.j], fun = 'cm', grid.method = 'halton', truncation = TRUE)
#       a = cbind(a, a_j)
#     }
#     #compute the covariance matrix of the normal scores
#   }
#   if (fun=='copula'){
#     a = npn(X, npn.func = "truncation")
#   }
#   if (fun=='gaussian'){
#     a = X
#   }
#   #a = a/sd(a[, 1])
#   Scov = cor(a) 
#   if(ridge){
#     scale = eigen(Scov)$values[1]
#     ignore = rho*scale
#     Scov = Scov + ignore*diag(dim(Scov)[1])
#     Scov = cov2cor(Scov)
#   }
#   roc_bglasso = numeric()
#   for(i in 1:length(gamma_n)){
#     est_prec_bglasso = bglasso(S = Scov, lambda = gamma_n[i], n = n, p = p, tol.Theta=1e-2, tol.beta=1e-3, 
#                                tol.beta_k=1e-7, maxit.Theta=10, maxit.beta=5, maxit.beta_k=1500, 
#                                Theta=NULL, Theta.inv=NULL, method=1)$Theta
#     est_adj = screenS(est_prec_bglasso, p = p, lambda = 0)$A
#     fptp = roc(est = est_adj, pop = Omega, thres = 0)
#     roc_bglasso = cbind(roc_bglasso, fptp)
#   }
#   return(roc_bglasso)
# }


roc_curve_bglasso = function(S, n, p, gamma_n, Omega, tol.Theta=1e-2, tol.beta=1e-3, 
                             tol.beta_k=1e-7, maxit.Theta=10, maxit.beta=5, maxit.beta_k=1500, 
                             Theta=NULL, Theta.inv=NULL, method=1){
  roc_bglasso = numeric()
  adj_matrix_list = list()
  for(k in 1:length(gamma_n)){
    A = bglasso(S = S, lambda = gamma_n[k], n = n, p = p, tol.Theta=tol.Theta, tol.beta=tol.beta, 
                               tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta, 
                               maxit.beta_k=maxit.beta_k, Theta=Theta, Theta.inv=Theta.inv, method=method)$alpha
    # A = matrix(0, nrow=p, ncol=p)
    # M = dim(S)[1]/p
    # 
    # for(i in 1:(p-1)){
    #   index.i = ((i-1)*M+1):(i*M)
    #   for(j in (i+1):p){
    #     index.j = ((j-1)*M+1):(j*M)
    #     
    #     if(norm(est_prec_bglasso[index.i, index.j], type="F") > 10^-1){
    #       A[i,j] = 1
    #     }      
    #   }
    # }
    # A = A + t(A) + diag(p)
    adj_matrix_list[[k]] = A
    fptp = roc(est = A, pop = Omega, thres = 0)
    roc_bglasso = cbind(roc_bglasso, fptp)
  }
  return(list(roc_bglasso = roc_bglasso, adj_matrix_list = adj_matrix_list))
}


# a function to calculate TP, FP, TN, FN
# and return precision, recall, TPR and FPR.
# Input:
#   G.true, the true p*p adjacency matrix
#   G.mat, the estimated p*p adjacency matrix
#   type,
#     AND: when two nodes both recognize each other as neighbors
#     OR: when either of them recognize each other as neighbors
# Output:
#   prec, precision
#   rec, recall
#   TPR and FPR, as their names

prec.rec <- function(G.true, G.mat, type=c("AND","OR")){
  
  p <- nrow(G.true)
  TP <- 0; TN <- 0; FP <- 0; FN <- 0
  
  if(type=="AND"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 & G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 & G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 & G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 & G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }else if(type=="OR"){
    for(i in 1:p){
      for(j in 1:p){
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==1 | G.mat[j,i]==1)) TP <- TP + 1
        if(i!=j & G.true[i,j]==1 & (G.mat[i,j]==0 | G.mat[j,i]==0)) FN <- FN + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==1 | G.mat[j,i]==1)) FP <- FP + 1
        if(i!=j & G.true[i,j]==0 & (G.mat[i,j]==0 | G.mat[j,i]==0)) TN <- TN + 1
      }
    }
  }
  prec <- TP / (TP + FP)
  if(TP+FP==0) prec <- 1
  
  rec <- TP / (TP + FN)
  if(TP+FN==0) rec <- 0
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(list(prec=prec, rec=rec, TPR=TPR, FPR=FPR))
}
