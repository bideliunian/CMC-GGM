#######################
## Model_b ROC (Model 2)
## p = 100, n = 200
## gaussian graphical model
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_b"
save.path <- "~/work/CMC-GGM/Model_b/Results/p100"

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
trans.type = 'no_transform'
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
data.cmcopula = mult.trans(X=data.X, p = p, fun = 'cmcopula')
data.copula = mult.trans(X=data.X, p = p, fun = 'copula')
data.linear = mult.trans(X=data.X, p = p, fun = 'scale')

## corresponding sample covariance matrix
S.cmcopula = cov(data.cmcopula)
S.copula = cov(data.copula)
S.linear = cor(data.linear)


#########################
## PART 3: Thresholding
###########################
est.prec.cmc.thres = vggm.thr(X=S.cmcopula, p = p, rho=10^-2, ridge = TRUE)$Theta0
est.prec.copula.thres = vggm.thr(X=S.copula, p = p, rho=10^-2, ridge = TRUE)$Theta0
est.prec.linear.thres = vggm.thr(X=S.linear, p = p, rho=10^-2, ridge = TRUE)$Theta0

## roc_curve_thresholding
roc.cmc.thres = roc_curve(est.prec.cmc.thres, G, n_grid = L)
roc.copula.thres = roc_curve(est.prec.copula.thres, G, n_grid = L)
roc.linear.thres = roc_curve(est.prec.linear.thres, G, n_grid = L)


########################
## PART 4: bglasso
#########################


## tuning parameter
lambda.cmc.max = max(S.cmcopula)
lambda.copula.max = max(S.copula)
lambda.linear.max = max(S.linear)

lambda.cmc = exp(seq(log(lambda.cmc.max), log(5*10^-2*lambda.cmc.max), length = L))
lambda.copula = exp(seq(log(lambda.copula.max), log(5*10^-2*lambda.copula.max), length = L))
lambda.linear = exp(seq(log(lambda.linear.max), log(5*10^-2*lambda.linear.max), length = L))


## first eigenvalue
rho.cmc = eigen(S.cmcopula)$values[1]
rho.copula = eigen(S.copula)$values[1]
rho.linear = eigen(S.linear)$values[1]


## roc_curve_bglasso
# parameters needed to run FGGM 
tol.Theta <-1e-2
tol.beta <-1e-2
maxit.Theta <-5
maxit.beta <-5
tol.beta_k <-1e-4
maxit.beta_k <-1500
method <- 1

cl <- makeCluster(20)
registerDoParallel(cl)

## bglasso with cmcopula transformation
roc.cmc.bglasso <- foreach(l = 1:L, .combine="cbind", .packages=c("matrixcalc","igraph","BB")) %dopar% {
  G.bglasso <- bglasso(S = S.cmcopula+10^-4*rho.cmc*diag(M*p), lambda = lambda.cmc[l], n = n, p = p,
                    tol.Theta=tol.Theta, tol.beta=tol.beta,
                    tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                    maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}


roc.copula.bglasso <- foreach(l = 1:L, .combine="cbind", .packages="matrixcalc") %dopar% {
  G.bglasso <- bglasso(S = S.copula+10^-4*rho.copula*diag(M*p), lambda = lambda.copula[l], n = n, p = p,
                    tol.Theta=tol.Theta, tol.beta=tol.beta,
                    tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                    maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}

roc.linear.bglasso <- foreach(l = 1:L, .combine="cbind", .packages="matrixcalc") %dopar% {
  G.bglasso <- bglasso(S = S.linear+10^-4*rho.linear*diag(M*p), lambda = lambda.linear[l], n = n, p = p,
                    tol.Theta=tol.Theta, tol.beta=tol.beta,
                    tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                    maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.l <- roc(est = G.bglasso, pop = G, thres = 0)
  TFPR.l
}

stopCluster(cl)


###############################
## PART 5: Collect and save data
###############################
df = data.frame(trans.methods = rep(c("cmc","copula",'linear',"cmc","copula",'linear'),each=L),
                 methods = rep(c("thres","bglasso"),each=3*L),
                 TP = cbind(roc.cmc.thres, roc.copula.thres, roc.linear.thres,
                            roc.cmc.bglasso, roc.copula.bglasso, roc.linear.bglasso)[1,],
                 FP = cbind(roc.cmc.thres, roc.copula.thres, roc.linear.thres,
                           roc.cmc.bglasso, roc.copula.bglasso, roc.linear.bglasso)[2,])

# df$TP[df$FP == 0] = 0
# plot<-ggplot(df, aes(x=FP, y=TP)) +
#   geom_line(aes(color=trans.methods,linetype=methods))
# plot

save(df, file=paste(save.path,"/ROC.b.100.gauss.seed",rand.seed,".Rdata",sep=""))
