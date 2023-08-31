#####################
## VCGGM 
#####################
## n sample size
## p number of nodes
## m dimension of vectors on each nodes
## X is a vector with dim n*(m*p)

require(MASS)

#########################
## generate gaussian data with
## sample size: n
## mean: mu 
## precision matrix: Omega
######################
rMVNormP = function(n, mu, Omega){
  p = length(mu)
  Z = matrix(rnorm(p*n), p, n)
  U = chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X = backsolve(U, Z) # more efficient and stable than acctually inverting
  X = sweep(X, 1, mu, FUN=`+`)
  
  return(t(X))
}


##########################
## copula transformation
##############################
copula_trans = function(x, type, alpha = 3, std = TRUE){
  ## Paras:
  ##  x: 1-dim array
  ##  m: > 1
  ##  type: {'power', 'exp', 'cdf'}
  ##  alpha: para for power tranformation
  
  if(type == 'power'){
    y = x^alpha
  }
  if(type == 'logit'){
    y = exp(x)/(1+exp(x))
  }
  if(type == 'cdf'){
    y = pnorm(x)
  }
  if(type == 'exp'){
    y = exp(x)
  }
  
  if(std){
    y = y - mean(y)
    y = y * sd(x) /sd(y)
  }
  
  return(y)
}

############################
## generate data
#############################
## Input:
## n: sample size
## p: number of nodes
## m: dimension of blocks
## model: 
##     Model 1: AR2 Theta_{j, j-1} = 0.4*I_m Theta_{j, j-2} = 0.2*I_m
##     Model 3: 2 hubs connected
##     Model 2: Random block-sparse
##     Model 4: Theta_{1, j} = 0.2*I_m
## trans.type: 'copula', 'cmcopula', 'pcmcopula'
## trans.func: 'power', 'exp', 'logit', 'cdf'
## Output:
##      X: n by p*d matrix
############################################
gen_data = function(n, p, m, k=2, model, run.ind, trans.type = NULL, trans.func='exp'){
  
  if (m < 2) {stop('m must be greater or equal to 2')}
  if (p < 2) {stop('p must be greater or equal to 2')}
  if (trans.type == 'pcmcopula'){
    if (m <= k) stop('m must be greater than k')
  }
  
  Omega = matrix(0, nrow = m*p, ncol = m*p)
  Omega0 = matrix(0, nrow = p, ncol = p)
  
  if(model == 'model1'){
    A = matrix(0, m, m)
    A[1,1] = A[2,2] = 1
    for (j in 1:(p-1)) {
      index.j = (((j-1)*m+1):(j*m))
      Omega[index.j,(index.j+m)] = 0.4*A 
      Omega0[j, j+1] = 1
      if(j != p-1) {
        Omega[index.j,index.j+2*m] = 0.2*A
        Omega0[j, j+2] = 1
      } 
    }
    Omega = Omega + t(Omega) + diag(m*p)
    Omega0 = Omega0 + t(Omega0)
  }
  
  if(model == 'model2'){
    delta = 0.3
    set.seed(run.ind)
    for(i in 1:(p/2-1)){
      index.i = (((i-1)*m+1):(i*m))
      for(j in (i+1):(p/2)){
        index.j = (((j-1)*m+1):(j*m))
        c = rbinom(n=1, size = 1, prob = 0.1)
        if(c){
          Omega[index.i, index.j] = 0.3*diag(m)
          Omega0[i, j] = 1
        }
      }
    }
    for(i in (p/2+1):(p-1)){
      index.i = (((i-1)*m+1):(i*m))
      for(j in (i+1):p){
        index.j = (((j-1)*m+1):(j*m))
        c = rbinom(n=1, size = 1, prob = 0.1)
        if(c){
          Omega[index.i, index.j] = 0.3*diag(m)
          Omega0[i, j] = 1
        }
      }
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 1.0, 1.0)
    Omega0 = Omega0 + t(Omega0) + diag(p)
  }
  
  if(model == 'model3'){
    set.seed(run.ind)
    seed = run.ind
    G.true = edge.gen(p=p, sparsity=0.99, hubnumber=2, hubsparsity=0.5)$G
    Omega = Omega.gen(G=G.true, M=m)
    Omega = hub.trans(Omega, p)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 1.0, 1.0)
    Omega0 = G.true
  }
  
  if(model == 'model4'){
    for (j in 1:(p-1)) {
      index.j = (j*m+1):((j+1)*m)
      if(j <= p-1){
        Omega[(1:m),(index.j)] = 0.2*diag(m)
      }
      Omega0[1,(2:p)] = 1
    }
    Omega = Omega + t(Omega)
    diag(Omega) = 0
    ee = min(eigen(Omega, only.values = T)$values)
    diag(Omega) = ifelse(ee < 0, -ee + 0.1, 0.1)
    Omega0 = Omega0 + t(Omega0) + diag(p)
  }
  
  if(model == 'model5'){
    set.seed(run.ind)
    seed = run.ind
    
    G1 = edge.gen(p=p/2, sparsity=0.99, hubnumber=1, hubsparsity=0.5)$G
    G2 = edge.gen(p=p/2, sparsity=0.99, hubnumber=1, hubsparsity=0.5)$G
    G = matrix(0, p, p)
    G[1:(p/2), 1:(p/2)] = G1
    G[(p/2+1):p, (p/2+1):p] = G2
    Omega = Omega.gen.new(G=G, m=m)
    Omega0 = G
  }
  
  Sigma = solve(Omega)
  diag.sigma = matrix(0, m*p, m*p)
  for (j in 1:p){
    index.j = (((j-1)*m+1):(j*m))
    diag.sigma[index.j, index.j] = matpower(Sigma[index.j, index.j], -1/2)
  }
  
  cov.X = diag.sigma %*% Sigma %*% diag.sigma
  X.pre = rmvnorm(n, sigma = cov.X)
  X = X.pre
  
  ## whether take copula gaussian transformation
  if(trans.type == 'copula'){
    for (j in 1:p) {
      index.j = (((j-1)*m+1):(j*m))
      X[,index.j] = apply(X = as.matrix(X.pre[,index.j]), MARGIN = 2, copula_trans, type = trans.func)
    }
  }
  
  ## whether take cm-copula gaussian transformation
  if(trans.type == 'cmcopula'){
    # if (m > 2) {stop('m must be less or equal to 2')}
    for (j in 1:p) {
      index.j = (((j-1)*m+1):(j*m))
      B = diag(m)
      for (j in 1:floor(m/2)) {
        B[(2*j-1), 2*j] = 1
        B[2*j, (2*j-1)] = -1
      }
      B = B/ sqrt(2)
      
      X[,index.j] = X.pre[,index.j]%*%B
      X[,index.j] = apply(X = as.matrix(X[,index.j]), MARGIN = 2, copula_trans, type = trans.func)
      X[,index.j] = X[,index.j]%*%t(B)
    }
  }
  
  if(trans.type == 'pcmcopula'){
    k = 2
    for (j in 1:p) {
      index.j = ((j-1)*m+1):(j*m)
      index.j.1 = ((j-1)*m+1):((j-1)*m+k)
      index.j.2 = ((j-1)*m+k+1):(j*m)
      B = matrix(0, m, k)
      B[1:k,] = diag(k)
      for (j in 1:floor(k/2)) {
        B[(2*j-1), 2*j] = 1
        B[2*j, (2*j-1)] = -1
      }
      B[1:k,] = B[1:k,]/sqrt(2)
      
      # X.tmp = apply(X = X.pre[,index.j]%*%B, MARGIN = 2, copula_trans, 
      #               std = FALSE, type = trans.func)
      # X[,index.j.1] = X.tmp%*%t(B)[,1:k]
      # X[,index.j.2] = X.pre[,index.j.2]
      ## get orthogonal complement of U
      Pv = diag(m) - B %*% t(B)
      X.tmp = apply(X = X.pre[,index.j]%*%B, MARGIN = 2, copula_trans, 
                    std = FALSE, type = trans.func)
      X[,index.j] = X.tmp%*%t(B) + X.pre[,index.j] %*% Pv
    }
  }
  
  if(trans.type == 'pcmcopula1'){
    k = 1
    for (j in 1:p) {
      index.j = ((j-1)*m+1):(j*m)
      index.j.1 = ((j-1)*m+1):((j-1)*m+k)
      index.j.2 = ((j-1)*m+k+1):(j*m)
      B = matrix(0, m, k)
      B[1,1] = B[2,1] = 1/sqrt(2)
      
      ## get orthogonal complement of U
      Pv = diag(m) - B %*% t(B)
      X.tmp = apply(X = X.pre[,index.j]%*%B, MARGIN = 2, copula_trans, 
                    std = FALSE, type = trans.func)
      X[,index.j] = X.tmp%*%t(B) + X.pre[,index.j] %*% Pv
    }
  }
  
  return(list(X = X, Omega = Omega, Omega0 = Omega0, X.pre = X.pre))
}

