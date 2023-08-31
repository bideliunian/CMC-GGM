#################
## thresholding methods
##
## X: n by m*p matrix
## m: size of vector in each nodes
## fun: 'cmcopula': cm copula method
##      'copula': gaussian copula method
##      'gaussian': ggm
## rho: control the MP inverse accuracy
##################################
vggm.thr = function(X, p, fun='cmcopula', rho=NULL, ridge=FALSE){
  m = dim(X)[2]/p
  Theta0 = matrix(0,p,p)
  Sigma = matrix(0,m*p,m*p)
  if(isSymmetric(X)){
    Sigma = X
    }
  else{
    n = dim(X)[1]
    Theta0 = matrix(0,p,p)
    Sigma = matrix(0,m*p,m*p)
    if (fun=="cmcopula"){
      #cm transformation
      a = numeric()
      for (j in 1:p) {
        index.j = (((j-1)*m+1):(j*m))
        a_j = cmcopula(X[,index.j], func = 'cm', truncation = TRUE)
        a = cbind(a, a_j)
      }
    }
    if (fun=='copula'){
      a = npn(X, npn.func = "truncation")
    }
    if (fun=='gaussian'){
      a = X
    }
    #a = a/sd(a[, 1])
    Sigma = cov(a)
  }
  
  scale = eigen(Sigma)$values[1]
  ignore = scale*rho
  if(ridge){
    Theta = riginv(Sigma, ignore)
  }
  else{
    #Theta = mppower(Sigma,-1,ignore)
    Theta = solve(Sigma)
  }
  
  for (i in 2:p) {
    for (j in 1:(i-1)){
      mat = Theta[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)]
      Theta0[i,j] <- Theta0[j,i] <- norm(as.matrix(mat),type="F")
    }
  }
  return(list(Theta0 = Theta0, Theta = Theta))
}
