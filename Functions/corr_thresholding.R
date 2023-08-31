######################################
# correlation thresholding
######################
vggm.ct = function(X, m, fun = "mrank", lambda = NULL){
  n = dim(X)[1]
  p = dim(X)[2]/m
  result = matrix(0,p,p)
  S = matrix(0,m*p,m*p)
  bS = matrix(0,p,p)
  fit = list()


  if (fun=="mrank"){
    # cm transformation
    a = numeric()
    for (j in 1:p) {
      a_j = cmcopula(X[,((j-1)*m+1):(j*m)], func = "cm")
      a = cbind(a, a_j)
    }
    #compute the covariance matrix of the normal scores
  }
  if (fun=='npn.truncation'){
    a = npn(X, npn.func = "truncation", verbose = FALSE)
  }
  if (fun=='gaussian'){
    a = X
  }
  S = cor(a)
  for (i in 2:p) {
    for (j in 1:(i-1)){
      mat = S[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)]
      bS[i,j] <- bS[j,i] <- norm(mat,type="F")
    }
  }
  rm(mat)

  S.rank = order(bS,decreasing = TRUE)

  if(is.null(lambda)){
    nlambda = 10
    lambda.min.ratio = 1

    density.max = lambda.min.ratio*p*(p-1)/2
    density.min = 1
    density.all = ceiling(seq(density.min,density.max,length = nlambda))*2

    fit$path = list()
    for(k in 1:nlambda){
      S.block = matrix(0,p,p)
      S.block[S.rank[1:density.all[k]]] = 1
      S.thr = S

      for (i in 2:p) {
        for (j in 1:(i-1)){
          if(S.block[i,j] == 0){
            S.thr[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)] = S.thr[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = 0
          }
        }
      }

      Theta.thr = solve(S.thr)
      Theta.block = matrix(0, p, p)

      for (i in 2:p) {
        for (j in 1:(i-1)){
          mat = Theta.thr[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)]
          Theta.block[i,j] <- Theta.block[j,i] <- norm(mat,type="F")
        }
      }

      fit$path[[k]] = Theta.block
    }
  }

  return(fit$path)
}


