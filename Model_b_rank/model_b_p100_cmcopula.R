#######################
## Model_b ROC (Model 1)
## p = 100, n = 300
## cmcopula
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_b_rank"
save.path <- "~/work/CMC-GGM/Model_b_rank/Results/p100"

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
n = 300
p = 100
M = 2
args = as.numeric(commandArgs(trailingOnly=TRUE))
rand.seed = args
L = 20 ## number of grids in roc curve
model = 'model2'
trans.type = 'cmcopula'
trans.func = 'exp'

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
data.rankcmcopula = mult.trans(X=data.X, p = p, fun = 'rankcmcopula')

## corresponding sample covariance matrix
S = cov(data.rankcmcopula)

#########################
## PART 3: thresholding
###########################
est.prec.thres = vggm.thr(X=S, p = p, rho=10^-1, ridge = TRUE)$Theta0
roc.thres = roc_curve(est.prec.thres, G, n_grid = L)

########################
## PART 4: group glasso
#########################

## tuning parameter
lambda.max = max(S)
lambda = exp(seq(log(lambda.max), log(5*10^-2*lambda.max), length = L))
rho = eigen(S)$values[1]

## group glasso with rank cmcopula transformation
tol.Theta <-1e-2
tol.beta <-1e-3
maxit.Theta <-10
maxit.beta <-5
tol.beta_k <-1e-7
maxit.beta_k <-1500
method <- 1

cl <- makeCluster(20)
registerDoParallel(cl)

roc.bglasso <- foreach(l = 1:L, .combine="cbind", .packages=c("matrixcalc","igraph","BB")) %dopar% {
  G.bglasso <- bglasso(S = S+10^-3*rho*diag(M*p), lambda = lambda[l], n = n, p = p,
                    tol.Theta=tol.Theta, tol.beta=tol.beta,
                    tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                    maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}

stopCluster(cl)

########################
## PART 5: neighborhood group lasso
roc.nb <- get.roc.nblasso(data.rankcmcopula, p, G, L)

###############################
## PART 6: Collect and save data
###############################

tpfp = cbind(roc.thres, roc.bglasso, roc.nb)
df = data.frame(trans.methods = rep('rankcmc',each=L),
                methods = rep(c("thresholding","group-glasso",'nbd-group-lasso'),each=3*L), 
                TP = tpfp[1,], FP = tpfp[2,])

save(df, file=paste(save.path,"/ROC.b.100.cmcopula.seed",rand.seed,".Rdata",sep=""))
