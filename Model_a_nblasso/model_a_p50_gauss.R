#######################
## Model_a ROC (Model 1)
## p = 50, n = 100
## ggm
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_a_nblasso"
save.path <- "~/work/CMC-GGM/Model_a_nblasso/Results"

## packages
library(pdist)
library(pracma)
library(clue)
library(huge)
library(igraph)
library(mvtnorm)
library(pdist)
library(randtoolbox)
library(gglasso)
library(Matrix)
library(pbmcapply)
library(parallel)

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

## global parameters
n = 200
p = 50
M = 2
L = 20 ## number of grids in roc curve
model = 'model1'
trans.type = 'no_trans'
n.expr = 50

#################### PART 1: helper functions ##########################
get.TPFP <- function(G.mat, L, G){
  p <- sqrt(nrow(G.mat))
  tpfp <- matrix(NA, L, 2)
  for (l in 1:L) {
    G.est <- matrix(G.mat[, l], nrow = p)
    tpfp[l, ] <- roc(est = G.est, pop = G, thres = 0)
  }
  tpfp
}

## transform to a list
mat2list <- function(data, M){
  data.list <- list()
  for (j in 1:M) {
    data.list[[j]] <- data[, j - M + M*(1:p)]
  }
  return(data.list)
}

get.roc.nblasso <- function(data.X, p, G, L, trans){
  M <- ncol(data.X) / p
  data.trans <- mult.trans(X=data.X, p = p, fun = trans)
  data.list <- mat2list(data.trans, M)
  G.est <- multivarNetwork(data.list, select="none",
                           nlambda=L, min.ratio=1e-3)
  
  G.and <- G.est$networks.and
  G.or <- G.est$networks.or
  
  roc.and <- get.TPFP(G.and, L, G)
  roc.or <- get.TPFP(G.or, L, G)
  
  return(cbind(roc.and, roc.or))
}

#################### PART 2: Training with neighborhood lasso #######################

time_start <- proc.time()
df.list <- list()
for (arg in 1:n.expr) {
  set.seed(arg)
  ##########################
  ## PART 1: Data Generation
  ###########################
  data = gen_data(n=n, p=p, m=M, model=model, run.ind = arg, trans.type = trans.type)
  data.X = data$X
  G = data$Omega0
  Omega = data$Omega
  
  roc.cmc <- get.roc.nblasso(data.X, p, G, L, 'cmcopula')
  roc.copula <- get.roc.nblasso(data.X, p, G, L, 'copula')
  roc.linear <- get.roc.nblasso(data.X, p, G, L, 'scale')
  
  df <- data.frame(trans.methods = rep(c("cmc","copula",'linear',"cmc","copula",'linear'),each=L),
                   methods = rep(c("and","or"),each=3*L),
                   TP = c(roc.cmc[, 1], roc.copula[, 1], roc.linear[, 1], 
                          roc.cmc[, 3], roc.copula[, 3], roc.linear[, 3]),
                   FP = c(roc.cmc[, 2], roc.copula[, 2], roc.linear[, 2], 
                          roc.cmc[, 4], roc.copula[, 4], roc.linear[, 4])
  )
  
  df.list[[arg]] <- df
}

time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

save(df.list, file=paste(save.path,"/ROC.a.050.gauss.nbd.Rdata",sep=""))

#################################
# using ADMM in Zhao(2021)
# cl <- makeCluster(20)
# registerDoParallel(cl)
# 
# get.lambda <- function(data, p, L){
#   lambda.max <- lambda.sup(data, p)
#   lambda.min <- 0.01
#   lambdas <- seq(lambda.max, lambda.min, length.out=L)
#   return(lambdas)
# }
# 
# lambdas.cmcopula <- get.lambda(data.cmcopula, p, L)
# roc.cmc.nblasso <- foreach(l = 1:L, .combine="rbind") %dopar% {
#   V <- neighbor.vtv.glasso(data.cmcopula, p, lambdas.cmcopula[l])
#   V.sym <- V + t(V)
#   G.and <- G.or <- V.sym
#   G.and[G.and < 2] <- 0 
#   G.and[G.and > 0] <- 1
#   G.or[G.or > 0] <- 1
#   
#   TFPR.and <- roc(est = G.and, pop = G, thres = 0)
#   TFPR.or <- roc(est = G.or, pop = G, thres = 0)
#   
#   list(TFPR.and, TFPR.or)
# }
# 
# stopCluster(cl)

