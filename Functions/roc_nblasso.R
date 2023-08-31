#################### helper functions ##########################
get.TPFP <- function(G.mat, L, G){
  p <- sqrt(nrow(G.mat))
  tpfp <- matrix(NA, L, 2)
  for (l in 1:L) {
    G.est <- matrix(G.mat[, l], nrow = p)
    tpfp[l, ] <- roc(est = G.est, pop = G, thres = 0)
  }
  return(tpfp)
}

## transform to a list
mat2list <- function(data, M){
  data.list <- list()
  for (j in 1:M) {
    data.list[[j]] <- data[, j - M + M*(1:p)]
  }
  return(data.list)
}

######## nbd group lasso ################
get.roc.nblasso <- function(data.trans, p, G, L){
  M <- ncol(data.trans) / p
  lambda.max <- lambda.sup(data.trans, p)
  data.list <- mat2list(data.trans, M)
  G.est <- multivarNetwork(data.list, select="none", lambda.max=lambda.max,
                           nlambda=L, min.ratio=1e-3)
  
  G.and <- G.est$networks.and
  G.or <- G.est$networks.or
  
  roc.and <- get.TPFP(G.and, L, G)
  roc.or <- get.TPFP(G.or, L, G)
  
  return(cbind(roc.and, roc.or))
}