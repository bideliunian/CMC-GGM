###
library(BB)
######################
## algorthm to solve block glasso
########################

## bglasso
##########################
# bglasso = function(X, m, fun, gamma_n, thres, max_iter = 6, rho, ridge = FALSE){
#   ## X can be 1) n by p data matrix 2) p by p covriance matrix after transformation
#   if(isSymmetric(X)){
#     p = dim(X)[2]/m
#     S = X
#   }
#   else{
#     n = dim(X)[1]
#     p = dim(X)[2]/m
#     if (fun=="mrank"){
#       #compute the normal scores
#       a = numeric()
#       for (j in 1:p) {
#         a_j = cmcopula(X[,((j-1)*m+1):(j*m)])
#         a = cbind(a, a_j)
#       }
#       #compute the covariance matrix of the normal scores
#     }
#     if (fun=='npn.truncation'){
#       a = npn(X, npn.func = "truncation", verbose = FALSE)
#     }
#     if (fun=='gaussian'){
#       a = X
#     }
#     #a = a/sd(a[, 1])
#     S = cor(a)
#     if(ridge){
#       scale = eigen(S)$values[1]
#       ridge = scale*rho
#       S = S + ridge*diag(m*p) 
#     }
#   }
#   result = matrix(0,p,p)
#   ## step1: initialize
#   Theta = diag(m*p)
#   Sigma = diag(m*p)
#   sigma_j = matrix(0, nrow = m*(p-1), ncol = m)
#   initial_w_j = matrix(0, nrow = m*(p-1), ncol = m)
#   ## step2
#   Sigma_old = 2*Sigma
#   Theta_old = 2*Theta
#   t = 1
#   while(sum((Theta-Theta_old)^2) > thres & t<=max_iter){
#     Theta_old = Theta
#     for (j in 1:p) {
#       index_j = ((j-1)*m+1):(j*m)
#       # (a)
#       # Theta_mj_inv = Sigma[-index_j, -index_j] - sigma_j%*%mppower(Sigma[index_j,index_j],-1,ignore = 10^-7)%*%t(sigma_j)
#       # ## (b) bcd
#       # S_jj = S[index_j,index_j]
#       # s_j = S[-index_j,index_j]
#       # w_j = bcd(initial_w_j, Theta_mj_inv, S_jj, s_j, gamma_n, thres)
#       # 
#       # 
#       # ## (c)
#       # U_j = Theta_mj_inv%*%w_j
#       # Sigma[index_j,index_j] = S[index_j,index_j]
#       # sigma_j = -U_j%*%S[index_j, index_j]
#       # Sigma[-index_j,-index_j] = Theta_mj_inv + U_j%*%S[index_j, index_j]%*%t(U_j)
#       # Sigma[-index_j,index_j] = sigma_j
#       # Sigma[index_j,-index_j] = t(sigma_j)
#       
#       
#        Theta_mj_inv = mppower(Theta[-index_j, -index_j],-1,ignore = 10^-3)
#        #Theta_mj_inv = solve(Theta[-index_j, -index_j])
#        S_jj = S[index_j,index_j]
#        s_j = S[-index_j,index_j]
#        w_j = bcd(initial_w_j, Theta_mj_inv, S[index_j,index_j], S[-index_j,index_j], gamma_n, thres)
#        Theta[-index_j,index_j] = w_j
#        Theta[index_j, -index_j] = t(w_j)
#        Theta[index_j,index_j] = solve(S_jj)+t(w_j)%*%Theta_mj_inv%*%w_j
#       
#     }
#     t = t+1
#     print(t)
#   }
#   #Theta = mppower(Sigma,-1,10^-7)
#   for (i in 2:p) {
#     for (j in 1:(i-1)){
#       mat = Theta[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)]
#       result[i,j] <- result[j,i] <- norm(mat,type="F")
#     }
#   }
#   rm(m, p, S)
#   return(list(result = result, Theta = Theta))
# }





###################
## block coordinate descent 
###########################

# bcd = function(initial, Theta_mj_inv, S_jj, s_j, gamma_n, thres, max_iter=5){
#   ## w_j: m*(p-1) by m
#   ## Theta_mj_inv: m*(p-1) by m*(p-1)
#   ## S_jj: m by m
#   m = dim(S_jj)[1]
#   p = dim(initial)[1]/m + 1
#   w_j = initial
#   old_w_j = initial + sqrt(thres)
#   t = 1
#   while (sqrt(sum((w_j-old_w_j)^2)) > thres & t<=max_iter) {
#     old_w_j = w_j
#     for (k in 1:(p-1)) {
#       index_k = ((k-1)*m+1):(k*m)
#       w_jk = w_j[index_k,]
#       r_jk = t(Theta_mj_inv[-index_k,index_k])%*%w_j[-index_k,]%*%S_jj + s_j[index_k,] ## r_jk: m by m
#       if(sqrt(sum(r_jk^2))<=gamma_n){
#         w_jk = 0
#       }
#       else{
#         #fn = function(x) crossprod(kronecker(Theta_mj_inv[index_k,index_k],S_jj)%*%x + c(r_jk) + gamma_n/sqrt(sum(x^2))*x)
#         #w_jk = optim(rep(0.1,m^2), fn, upper = 10^4*rep(1,m^2), lower = -10^4*rep(1, m^2), method="L-BFGS-B")$par
#         fn <- function(x) return (kronecker(Theta_mj_inv[index_k,index_k],S_jj)%*%x + c(r_jk) + gamma_n/sqrt(sum(x^2))*x)
#         w_jk = nleqslv(rep(0.1,m^2), fn, xscalm = "auto",global = "dbldog",method = "Broyden")$x
#       }
#       w_j[index_k,] = matrix(w_jk, nrow = m, ncol = m)
#     }
#     t = t+1
#   }
#   return(w_j)
# }

# bglasso_1 = function(X, m, fun, gamma_n, thres, max_iter = 10, ridge = FALSE, rho = 1){
#   ## X can be 1) n by p data matrix 2) p by p covriance matrix after transformation
#   if(isSymmetric(X)){
#     p = dim(X)[2]/m
#     S = X
#   }
#   else{
#     n = dim(X)[1]
#     p = dim(X)[2]/m
#     if (fun=="mrank"){
#       #compute the normal scores
#       a = numeric()
#       for (j in 1:p) {
#         a_j = cmcopula(X[,((j-1)*m+1):(j*m)])
#         a = cbind(a, a_j)
#       }
#       #compute the covariance matrix of the normal scores
#     }
#     if (fun=='npn.truncation'){
#       a = npn(X, npn.func = "truncation", verbose = FALSE)
#     }
#     if (fun=='gaussian'){
#       a = X
#     }
#     #a = a/sd(a[, 1])
#     S = cor(a)
#     if(ridge){
#       scale = eigen(S)$values[1]
#       ridge = scale*rho
#       S = S + ridge*diag(m*p)
#     }
#   }
#   result = matrix(0,p,p)
#   ## step1: initialize
#   Theta = diag(m*p)
#   Sigma = diag(m*p)
#   sigma_j = matrix(0, nrow = m*(p-1), ncol = m)
#   initial_w_j = matrix(0, nrow = m*(p-1), ncol = m)
#   ## step2
#   Sigma_old = 2*Sigma
#   Theta_old = 2*Theta
#   t = 1
#   while(sum((Sigma-Sigma_old)^2) > thres & t<=max_iter){
#     Sigma_old = Sigma
#     for (j in 1:p) {
#       index_j = ((j-1)*m+1):(j*m)
#       ## (a)
#       Theta_mj_inv = Sigma[-index_j, -index_j] - sigma_j%*%mppower(Sigma[index_j,index_j],-1,ignore = 10^-7)%*%t(sigma_j)
#       ## (b) bcd
#       S_jj = S[index_j,index_j]
#       s_j = S[-index_j,index_j]
#       w_j = bcd(initial_w_j, Theta_mj_inv, S_jj, s_j, gamma_n, thres)
#
#
#       ## (c)
#       U_j = Theta_mj_inv%*%w_j
#       Sigma[index_j,index_j] = S[index_j,index_j]
#       sigma_j = -U_j%*%S[index_j, index_j]
#       Sigma[-index_j,-index_j] = Theta_mj_inv + U_j%*%S[index_j, index_j]%*%t(U_j)
#       Sigma[-index_j,index_j] = sigma_j
#       Sigma[index_j,-index_j] = t(sigma_j)
#     }
#     t = t+1
#     print(t)
#   }
#   Theta = mppower(Sigma,-1,10^-7)
#   for (i in 2:p) {
#     for (j in 1:(i-1)){
#       mat = Theta[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)]
#       result[i,j] <- result[j,i] <- norm(mat,type="F")
#     }
#   }
#   rm(m, p, S)
#   return(list(result = result, Theta = Theta))
# }


#######################
## solve a M^2 nonlinear equation using BB package
## beta.k is the initial value, cannot be zero
## B, b and lambda are constants
## output is a M^2 matrix solution

nleq <-function(B, b, lambda, beta.k, maxit.beta_k, tol.beta_k){
  M <- sqrt(length(b))
  
  froth <-function(x){
    f <- rep(NA,length(x))
    for(i in 1:length(x)){
      f[i] <- sum(B[i,] * x) + b[i] +lambda*x[i]/sqrt(sum(x^2))
    }
    f
  }
  
  ## if zero initial value, denominator gives a small permutation
  if(sum(abs(beta.k)) == 0){
    beta.k <- matrix(1e-6, nrow=M, ncol=M)
  }
  
  eq <- BBsolve(par=as.vector(beta.k), fn=froth, quiet=T, control=list(maxit=maxit.beta_k, tol=tol.beta_k)) 
  
  ans <- matrix(eq$par, nrow=M, ncol=M)
  return(ans)
}



################################################
## Coordinate Descent Algorithm
######################################################

#### Input:
## p: Number of variables
## S: Sample covariance matrix Mp by Mp
## lambda: tunning parameter
## tol.Theta: tolerace value for matrix Theta
## iter.Theta.max: maximum iterative steps for solving Theta
## Witten's idea in high dimensional setting
## Theta: initial precision matrix, default is 1/diagonal(Sample Cov)
## Theta.inv: initial covariance matrix, default is diagonal(Sample Cov) 
## method: 1(general) 2(simplified)

#### Output
## Theta: Inverse covariance matrix Mp by Mp
## alpha: Indicator matrix for network strucure
## crit: criterion as iteration goes, increasing sequence
## aic, bic: information criterion

bglasso <-function(S, lambda, n, p, tol.Theta=1e-2, tol.beta=1e-3, tol.beta_k=1e-7, 
                   maxit.Theta=10, maxit.beta=10, maxit.beta_k=1500, Theta=NULL, Theta.inv=NULL, method=1){
  
  M <-dim(S)[1]/p ## Block dimension
  
  ## Initialization, retrieve the block diagonal information from sample covariance matrix
  Sd <-list()
  Sd.inv <-list()
  for(i in 1:p){
    index.i <-((i-1)*M+1):(i*M)
    Sd[[i]] <-S[index.i, index.i]
    Sd.inv[[i]] <-solve(S[index.i, index.i])
  }
  if(is.null(Theta)){
    Theta <- diag(diag(1/S))
  }
  if(is.null(Theta.inv)){
    Theta.inv <- diag(diag(S))
  }
  
  dif.Theta <-1
  iter.Theta <-1
  alpha <-matrix(0, p, p) ## Indicator matrix
  crit <-NULL
  
  ## Current step
  while((dif.Theta >= tol.Theta) & (iter.Theta <= maxit.Theta)){ 
    
    Theta.old <-Theta
    
    ## Sweeping each column block from 1 to p
    for(j in 1:p){
      index.j <- ((j-1)*M+1):(j*M)
      s.jj <- S[index.j, index.j] ## M by M
      s.nj <- S[-index.j, index.j]  ## (p-1)M by M
      
      ## Formula (9): Solve Theta.j.inv given Theta.inv
      Theta.j.inv <- Theta.inv[-index.j, -index.j]- Theta.inv[-index.j, index.j] %*% solve(Theta.inv[index.j, index.j]) %*% t(Theta.inv[-index.j, index.j])
      beta <- Theta[-index.j, index.j]		

      ##Solve beta=w_j: (p-1) M by M blocks given a fixed Theta.j.inv by inner Coordinate Descent Algorithm 
      dif.beta <-1
      iter.beta <-1
      while((dif.beta >= tol.beta) & (iter.beta <= maxit.beta)){     
        
        beta.old <-beta  
        
        for(k in 1:(p-1)){        
          index.k <-((k-1)*M+1):(k*M)        
          A.kk <-Theta.j.inv[index.k, index.k]
          r.k <-t(Theta.j.inv[-index.k, index.k]) %*% beta[-index.k,] %*% s.jj + s.nj[index.k,]
          B <-t(s.jj) %x% A.kk
          b <-as.vector(r.k)
          
          if(method==c(1)){
            if(norm(r.k, type="F") <= lambda){
              alpha[-j, j][k] = alpha[j, -j][k] <-0
              beta[index.k, ] <-matrix(0,M,M)
            } 
            
            else{
              alpha[-j, j][k] = alpha[j, -j][k] <-1		 
              
              ## initialization for beta.k
              if((iter.Theta == 1) & (iter.beta == 1)){
                beta.k <- -solve(A.kk) %*% (r.k+lambda) %*% solve(s.jj)
              }
              else{
                beta.k <-beta[index.k,] ## warm start using beta from previous step
              }
              
              beta.k <-nleq(B, b, lambda, beta.k, maxit.beta_k, tol.beta_k)       
              beta[index.k,] <-beta.k      
            }
          }
          ## end of method 1
          
          if(method== c(2)){
            if((iter.Theta ==1)){
              if(norm(r.k, type="F") <= lambda){
                alpha[-j, j][k] = alpha[j, -j][k] <-0
                beta[index.k, ] <-matrix(0,M,M)
              } 
              else{
                alpha[-j, j][k] = alpha[j, -j][k] <-1
                if((iter.beta ==1)){
                  beta.k <- -solve(A.kk) %*% (r.k+lambda) %*% solve(s.jj)
                }
                else{
                  beta.k <-beta[index.k,]
                }
                
                beta.k <-nleq(B, b, lambda, beta.k, maxit.beta_k, tol.beta_k)        
                beta[index.k,] <-beta.k
              }
            }
            ## do method 1 on the first iterative step for Theta
            
            else{
              beta.k <-beta[index.k,]
              
              if((sum(abs(beta.k)) == 0)){
                alpha[-j, j][k] = alpha[j, -j][k] <-0 
                beta[index.k, ] <-matrix(0,M,M)
              }
              else{
                if(norm(r.k, type="F") <= lambda){
                  alpha[-j, j][k] = alpha[j, -j][k] <-0 
                  beta[index.k, ] <-matrix(0,M,M)
                }
                else{
                  alpha[-j, j][k] = alpha[j, -j][k] <-1
                  beta.k <-nleq(B, b, lambda, beta.k, maxit.beta_k, tol.beta_k)        
                  beta[index.k,] <-beta.k
                }
              }
            }
          }
          ## end of method 2          
        }
        ## end of k
        
        if(sum(abs(beta.old)) == 0){
          beta.prev <-matrix(1e-6, nrow=(p-1)*M, ncol=M) ## this appears in the denominator
          dif.beta <-norm(beta-beta.prev,type="F")/norm(beta.prev,type="F")
        }
        else{
          dif.beta <-norm(beta-beta.old,type="F")/norm(beta.old,type="F")
        }
        iter.beta <-iter.beta+1  
        
        #print(c(iter.beta, iter.Theta, j))
      }
      ## end of while for beta      
      
      Theta[-index.j, index.j] <- beta	
      Theta[index.j, -index.j] <- t(beta) ## row block by symmetry	  
      
      ##Solve w_{jj} by Formula (13) given Theta.j.inv
      Theta[index.j , index.j] <-solve(s.jj) + t(beta) %*% Theta.j.inv %*% beta
      
      ## Step 2 (c): Update Theta.inv given Theta.j.inv
      U <-Theta.j.inv %*% beta
      E <-solve(Theta[index.j,index.j]- t(beta) %*% U)
      Theta.inv[-index.j, -index.j] <-Theta.j.inv + U %*% E %*% t(U)
      Theta.inv[-index.j, index.j] <- -U %*% E
      Theta.inv[index.j, -index.j] <- t(Theta.inv[-index.j, index.j])
      Theta.inv[index.j, index.j] <- E
    }
    ## end of j
    
    sumTheta <-0
    for(l in 1:(p-1)){
      index.l <-((l-1)*M+1):(l*M)
      for(m in (l+1):p){
        index.m <-((m-1)*M+1):(m*M)
        sumTheta <-sumTheta + norm(Theta[index.l,index.m], type="F")
      }
    }
    loglik <- log(det(Theta)) - sum(diag(S%*%Theta))
    crit <-c(crit, loglik - 2*lambda*sumTheta)  ## an increasing sequence
    dif.Theta <-norm(Theta-Theta.old, type="F")/norm(Theta.old, type="F")
    iter.Theta <-iter.Theta + 1
    
    ## plot network structure
    #g1 = graph.adjacency(alpha,mode="undirected")
    #plot(g1)
    #print(iter.Theta)
  }
  ## end of while for Theta
  
  aic <- -loglik + 2*(M^2*sum(alpha)/2 + M*(M+1)*p/2)/n
  bic <- -loglik + log(n)*(M^2*sum(alpha)/2 + M*(M+1)*p/2)/n
  
  ans <-list(Theta=Theta, Theta.inv=Theta.inv, alpha=alpha, iter.Theta=iter.Theta-1, crit=crit, aic=aic, bic=bic)
  return(ans)
}


