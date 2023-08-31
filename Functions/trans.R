#####################################
## non para normal transformation
#################################
## input:
## x: The n by d data matrix representing n observations in d dimensions
## npn.func: The transformation function used in the npn transformation; "shrinkage","truncation","skeptic"
## npn.thresh: The truncation threshold used in nonparanormal transformation, ONLY applicable when npn.func = "truncation"
## verbose: If verbose = FALSE, tracing information printing is disabled.
npn = function (x, npn.func = "truncation", npn.thresh = NULL, verbose = FALSE) {
  n = nrow(x)
  d = ncol(x)
  if (npn.func == "shrinkage") {
    if (verbose) 
      cat("Conducting the nonparanormal (npn) transformation via shrunkun ECDF....")
    x = qnorm(apply(x, 2, rank)/(n + 1))
    x = x/sd(x[, 1])
    if (verbose) 
      cat("done.\n")
    rm(n, d, verbose)
  }
  if (npn.func == "truncation") {
    if (verbose) 
      cat("Conducting nonparanormal (npn) transformation via truncated ECDF....")
    if (is.null(npn.thresh)) {npn.thresh = 1/(4 * (n^0.25) * sqrt(pi * log(n)))}
    x = qnorm(pmin(pmax(apply(x, 2, rank)/n, npn.thresh), 
                   1 - npn.thresh))
    x = x/sd(x[, 1])
    if (verbose) 
      cat("done.\n")
    rm(n, d, npn.thresh, verbose)
  }
  if (npn.func == "skeptic") {
    if (verbose) 
      cat("Conducting nonparanormal (npn) transformation via skeptic....")
    x = 2 * sin(pi/6 * cor(x, method = "spearman"))
    if (verbose) 
      cat("done.\n")
    rm(n, d, verbose)
  }
  return(x)
}

#################
## cmc transformation
##############################
## x: n by d matrix
## return: n by d matrix
cmcopula = function(x, fun, grid.method = 'random', truncation = TRUE){
  x = as.matrix(x)
  n = nrow(x)
  m = ncol(x)
  ## whether to standardize x
  sig = matpower(cov(x), 1/2)
  mu = colMeans(x)

  if(m==1){
    grid = seq(from = 1/(n+1), to = n/(n+1), by = 1/(n+1))
  }
  else{
    ## three ways to generate regular points: halton; regular; random
    if (grid.method == 'halton'){
      grid = halton(n, m)
      grid = qnorm(grid)
    }
    if (grid.method == 'regular'){
      n_grid = ceiling(n^{1/m})
      a = seq(0,1,length.out = n_grid + 2)[-c(1,n_grid + 2)]
      grid = as.matrix(expand.grids(a,m))
      grid = qnorm(grid)
    }
    if (grid.method == 'random'){
      grid = matrix(rnorm(n*m), n, m)
    }
  }
  
  grid = as.matrix(grid)
  distmat = matrix(0, nrow = n, ncol = n)
  
  for(i in 1:n){
    distmat[i,]=apply((x[i,]-t(grid)),2,Norm)^2 
  }
  
  if(fun == 'cm'){
    assignmentFUN = solve_LSAP(distmat)
    assignmentSOL = cbind(seq_along(assignmentFUN),assignmentFUN)
    x.transformed = as.matrix(grid[assignmentSOL[,2],])
  }
  
  else if(fun == 'sinkhorn'){
    skh = sinkhorn(x, grid, p = 2, lambda = 1)
    x.transformed = n*(skh$plan)%*%grid
  }
  
  if (truncation) {
    npn.thresh = sqrt(2 * log(n))
    x.transformed = pmin(pmax(as.matrix(x.transformed), -npn.thresh), npn.thresh)
  }
  
  return(x.transformed)
}

#################
## compositional cmc transformation/ rank cmc transformation
##############################
## x: n by d matrix
## return: n by d matrix
rankcmcopula = function(x, grid.method = 'random', truncation = TRUE){
  x = as.matrix(x)
  n = nrow(x)
  m = ncol(x)
  ## whether to standardize x
  sig = matpower(cov(x), 1/2)
  mu = colMeans(x)
  
  if(m==1){
    grid = seq(from = 1/(n+1), to = n/(n+1), by = 1/(n+1))
  }
  else{
    ## three ways to generate regular points: halton; regular; random
    if (grid.method == 'halton'){
      grid = halton(n, m)
    }
    if (grid.method == 'regular'){
      n_grid = ceiling(n^{1/m})
      a = seq(0,1,length.out = n_grid + 2)[-c(1,n_grid + 2)]
      grid = as.matrix(expand.grids(a,m))
    }
    if (grid.method == 'random'){
      grid = matrix(runif(n*m), n, m)
    }
  }
  
  grid = as.matrix(grid)
  distmat = matrix(0, nrow = n, ncol = n)
  
  for(i in 1:n){
    distmat[i,]=apply((x[i,]-t(grid)),2,Norm)^2 
  }
  
  assignmentFUN = solve_LSAP(distmat)
  assignmentSOL = cbind(seq_along(assignmentFUN),assignmentFUN)
  x.transformed = as.matrix(grid[assignmentSOL[,2],])
  
  if (truncation) {
    npn.thresh = 1/n
  }
  x.transformed = qnorm(pmin(pmax(x.transformed, npn.thresh), 
                 1 - npn.thresh))
  
  return(x.transformed)
}



#################
## blockwise nonparanormal transformation by projected cmc copula
##############################
## x: n by d matrix
## U: d by k matrix
## return: n by d matrix

pcmcopula = function(x, grid.method = 'random', truncation = FALSE, U){
  x = as.matrix(x)
  n = nrow(x)
  m = ncol(x)
  k = ncol(U)
  
  sig = matpower(cov(x), 1/2)
  mu = colMeans(x)
  
  ## normalize U
  U = U %*% matpower((t(U)%*%U), -1/2)
  
  ## get orthogonal complement of U
  Pv = diag(m) - U %*% t(U)
  
  ## three ways to generate regular points: halton; regular; random
  if (grid.method == 'halton'){
    grid = halton(n, m)
    grid = qnorm(grid)
  }
  if (grid.method == 'regular'){
    n_grid = ceiling(n^{1/m})
    a = seq(0,1,length.out = n_grid + 2)[-c(1,n_grid + 2)]
    grid = as.matrix(expand.grids(a,m))
    grid = qnorm(grid)
  }
  if (grid.method == 'random'){
    grid = matrix(rnorm(n*m), n, m)
  }
  
  if (truncation) {
    npn.thresh = sqrt(2 * log(n))
    grid = pmin(pmax(as.matrix(grid), -npn.thresh), npn.thresh)
  }
  
  distmat = matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      distmat[i,j] = Norm((x[i,]%*% U-grid[j,]%*% U),2)^2  
    }
  }
  
  assignmentFUN = solve_LSAP(distmat)
  assignmentSOL = cbind(seq_along(assignmentFUN),assignmentFUN)
  proj = as.matrix((grid%*% U)[assignmentSOL[,2],])
  
  #data.transformed = proj %*% t(U) + grid %*% Pv
  data.transformed = proj %*% t(U) + x %*% Pv
  
  return(data.transformed)
}

###########################################
## cm transformation by discrete to continuous ot
## only applied to dimension 2
#################################
cmcopula.d2c = function(x, ranks = TRUE, truncation = TRUE){
  x = as.matrix(x)
  n = nrow(x)
  M = ncol(x)
  
  X.OTM = tos.fit(x)
  
  if(ranks){
    R = tos.rank(X.OTM, x)
    ranks = R$Rank
  }
  else{
    ranks = X.OTM$Centroid
  }
  
  npn.thresh = 1/(4 * (n^0.25) * sqrt(pi * log(n)))
  data.transformed = qnorm(pmin(pmax(as.matrix(ranks), npn.thresh), 
                                1 - npn.thresh))
  
  return(data.transformed)
}


##############################
## unified function to apply multiple transformation methods
#####################################
mult.trans = function(X, p, m.list=c(0), fun, U.list=NULL, grid.method = 'random'){
  # Param:
  #   X n by p*m matrix
  #   fun {'cmcopula', 'cmcopula.d2c', 'copula', 'scale', 'pcmcopula'}
  #   U projection direction for pcmcopula, NULL by default
  # Return:
  #   Normal score
  if (!(fun %in% c('scale', 'copula', 'cmcopula', 'pcmcopula','rankcmcopula'))) 
  {stop('could not find transformation method')}
  n = dim(X)[1]
  if (m.list[1] == 0){
    m.list = rep(dim(X)[2]/p, p)
  }
  m = m.list[1]
  # if (m > 3) {grid.method = 'random'}
  # else {grid.method = 'halton'}
  
  m.cumsum = cumsum(m.list)
  m.cumsum =c(0, m.cumsum)
  
  normal.score = numeric()
  
  if (fun=='scale'){
    normal.score = X
    for (j in 1:p) {
      index.j = (m.cumsum[j]+1):m.cumsum[j+1]
      sig.nrt = matpower(cov(X[, index.j]), -1/2)
      normal.score[, index.j] = X[, index.j] %*% sig.nrt
    }
  }
  
  if (fun == 'copula'){
    normal.score = npn(X, npn.func = "truncation")
  }
  
  if (fun == "cmcopula"){
    for (j in 1:p) {
      index.j = (m.cumsum[j]+1):m.cumsum[j+1]
      a_j = cmcopula(X[,index.j], fun = 'cm', grid.method = grid.method, truncation = TRUE)
      normal.score = cbind(normal.score, a_j)
    }
  }
  
  if (fun == "rankcmcopula"){
    for (j in 1:p) {
      index.j = (m.cumsum[j]+1):m.cumsum[j+1]
      a_j = rankcmcopula(X[,index.j], grid.method = grid.method, truncation = TRUE)
      normal.score = cbind(normal.score, a_j)
    }
  }
  
  if (fun == 'pcmcopula'){
    for (j in 1:p) {
      index.j = (m.cumsum[j]+1):m.cumsum[j+1]
      U = as.matrix(U.list[[j]])
      a_j = pcmcopula(X[,index.j], grid.method = grid.method, truncation = TRUE, U = U)
      normal.score = cbind(normal.score, a_j)
    }
  }
  return(normal.score)
}

