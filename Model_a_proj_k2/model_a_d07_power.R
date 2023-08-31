#######################
## Model_a ROC (Model 1)
## p = 50, n = 200, d = 7
## pcmcopula
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_a_proj"
read.path <- "~/work/CMC-GGM/proj_cmc"
save.path <- "~/work/CMC-GGM/Model_a_proj/Results"

## packages
library(pdist)
library(pracma)
library(clue)
library(huge)
library(igraph)
library(BB)
library(mvtnorm)
library(doParallel)
library(pdist)
library(randtoolbox)
library(Matrix)
library(gglasso)

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

## global parameters
n = 200
p = 50
M = 7
k = 2
args = as.numeric(commandArgs(trailingOnly=TRUE))
rand.seed = args
L = 20 ## number of grids in roc curve
model = 'model1'
trans.type = 'pcmcopula'
trans.func = 'power'

##########################
## PART 1: Data Generation
###########################
data = gen_data(n=n, p=p, m=M, model=model, run.ind = rand.seed, trans.type = trans.type, trans.func=trans.func)
data.X = data$X
G = data$Omega0
Omega = data$Omega


###################
## PART 2: Three transformations: cmcopula, copula, linear
######################
U.df = read.csv(file=paste(read.path,"/proj_subspace_model1_exp_d07.csv",sep=""), header = FALSE)
U.list = list()
for (i in 1:p) {
  U.list[[i]] = matrix(as.numeric(U.df[i,]), nrow = M, ncol = k)
}

data.pcmcopula = mult.trans(X=data.X, p = p, fun = 'pcmcopula', U.list = U.list)
data.cmcopula = mult.trans(X=data.X, p = p, fun = 'cmcopula')
data.copula = mult.trans(X=data.X, p = p, fun = 'copula')

## corresponding sample covariance matrix
S.cmcopula = cov(data.cmcopula)
S.pcmcopula = cov(data.pcmcopula)
S.copula = cov(data.copula)


#########################
## PART 3: Thresholding
###########################
est.prec.cmc.thres = vggm.thr(X=S.cmcopula, p = p, rho=10^-3, ridge = TRUE)$Theta0
est.prec.pcmc.thres = vggm.thr(X=S.pcmcopula, p = p, rho=10^-3, ridge = TRUE)$Theta0
est.prec.copula.thres = vggm.thr(X=S.copula, p = p, rho=10^-3, ridge = TRUE)$Theta0

## roc_curve_thresholding
roc.cmc.thres = roc_curve(est.prec.cmc.thres, G, n_grid = L)
roc.pcmc.thres = roc_curve(est.prec.pcmc.thres, G, n_grid = L)
roc.copula.thres = roc_curve(est.prec.copula.thres, G, n_grid = L)


########################
## PART 4: bglasso
#########################


## tuning parameter
lambda.cmc.max = max(S.cmcopula)
lambda.pcmc.max = max(S.pcmcopula)
lambda.copula.max = max(S.copula)

lambda.cmc = exp(seq(log(lambda.cmc.max/2), log(lambda.cmc.max/5), length = L))
lambda.pcmc = exp(seq(log(lambda.pcmc.max/2), log(lambda.pcmc.max/5), length = L))
lambda.copula = exp(seq(log(lambda.copula.max), log(lambda.copula.max/10), length = L))


## first eigenvalue
rho.cmc = eigen(S.cmcopula)$values[1]
rho.pcmc = eigen(S.pcmcopula)$values[1]
rho.copula = eigen(S.copula)$values[1]


## roc_curve_bglasso
# parameters needed to run FGGM 
tol.Theta <-1e-3
tol.beta <-1e-3
maxit.Theta <-10
maxit.beta <-5
tol.beta_k <-1e-7
maxit.beta_k <-1500
method <- 1

time_start <- proc.time()

cl <- makeCluster(20)
registerDoParallel(cl)

## bglasso with cmcopula transformation
roc.cmc.bglasso <- foreach(l = 1:L, .combine="cbind", .packages=c("matrixcalc","igraph","BB")) %dopar% {
  G.bglasso <- bglasso(S = S.cmcopula+10^-3*diag(M*p), lambda = lambda.cmc[l], n = n, p = p,
                       tol.Theta=tol.Theta, tol.beta=tol.beta,
                       tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                       maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}

roc.pcmc.bglasso <- foreach(l = 1:L, .combine="cbind", .packages="matrixcalc") %dopar% {
  G.bglasso <- bglasso(S = S.pcmcopula+10^-3*diag(M*p), lambda = lambda.pcmc[l], n = n, p = p,
                       tol.Theta=tol.Theta, tol.beta=tol.beta,
                       tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                       maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}


roc.copula.bglasso <- foreach(l = 1:L, .combine="cbind", .packages="matrixcalc") %dopar% {
  G.bglasso <- bglasso(S = S.copula+10^-3*diag(M*p), lambda = lambda.copula[l], n = n, p = p,
                       tol.Theta=tol.Theta, tol.beta=tol.beta,
                       tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                       maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}

stopCluster(cl)


######## nbd group lasso ################
get.roc.nblasso <- function(data.trans, p, G, L){
  M <- ncol(data.X) / p
  data.list <- mat2list(data.trans, M)
  lambda.max <- lambda.sup(data.trans, p) / 10
  G.est <- multivarNetwork(data.list, select="none", lambda.max=lambda.max,
                           nlambda=L, min.ratio=5*1e-2)
  
  G.and <- G.est$networks.and
  G.or <- G.est$networks.or
  
  roc.and <- get.TPFP(G.and, L, G)
  
  return(t(roc.and))
}

roc.pcmc.nb <- get.roc.nblasso(data.pcmcopula, p, G, L)
roc.cmc.nb <- get.roc.nblasso(data.cmcopula, p, G, L)
roc.copula.nb <- get.roc.nblasso(data.copula, p, G, L)

time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

###############################
## PART 5: Collect and save data
###############################
tpfp = cbind(roc.pcmc.thres, roc.cmc.thres, roc.copula.thres, 
             roc.pcmc.bglasso, roc.cmc.bglasso, roc.copula.bglasso,
             roc.pcmc.nb, roc.cmc.nb, roc.copula.nb)

df = data.frame(trans.methods = rep(rep(c("pcmc","cmc","copula"), 3), each=L),
                methods = rep(c("thresholding","group-glasso",'nbd-group-lasso'),each=3*L), 
                TP = tpfp[1,], FP = tpfp[2,])

save(df, file=paste(save.path,"/ROC_a_d07_power.seed",rand.seed,".Rdata",sep=""))
