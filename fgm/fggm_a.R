library(pdist)
library(pracma)
library(clue)
library(ggplot2)
library(huge)
library(igraph)
library(BB)
library(mvtnorm)
library(randtoolbox)
library(mvtnorm)
library(fda)
library(matrixcalc)

############ pathes
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/fgm"
save.path <- "~/work/CMC-GGM/fgm/Results"


# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)
source(paste(working.path,"gen_fct_data.R", sep="/"))
# Global Parameter Settings

model <- 'model1'
m <- 5 # number of basis used to generate data
M <- 5
p <- 50
n <- 300
tau <- 100 # number of observations
L <- 20 # number of bases used to estimate a_ij hat
trans.type = 'copula'
trans.func = 'exp'

args <- (commandArgs(TRUE))
run.ind <- args
set.seed(args)
time.start <- proc.time()

####################################
#     PART 1: DATA GENERATION      #
####################################

#    Generate Random Functions and Observation h_ijk       
# h_ijk = Fourier basis func(t_k) %*% delta_ij + iid error

# 0. Generate delta
delta.list <- gen_data(n=n, p=p, m=m, model=model, run.ind = rand.seed, 
                       trans.type = trans.type, trans.func=trans.func)
delta <- delta.list$X
Omega <- delta.list$Omega
Omega0 <- delta.list$Omega0

# 1. Observation time
obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta

# 2. Fourier basis function for data generation
b.mat <- obs.fourier.bases(obs.time, m)

# 3. Observations h_ijk
h <- array(0, c(n, p, tau))
for(i in 1:n){
  for(j in 1:p){
    h[i,j,] <- b.mat %*% matrix(delta[i, ((j-1)*m+1) : (j*m)], ncol=1) + rnorm(tau, 0, 0.5)
  }
}

####################################
#     PART 2: GAIN FPC SCORE       #
####################################

# fpca by fda package
fpca.list2 <- list()
fd.list <- list()
for(j in 1:p){
  # bspline.basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=20)
  # fd.object.array <- Data2fd(argvals=obs.time, y=t(h[,j,]), basisobj=bspline.basis)
  fourier.basis <- create.fourier.basis(rangeval=c(0, 1), nbasis=10)
  fd.object.array <- Data2fd(argvals=obs.time, y=t(h[,j,]), basisobj=fourier.basis)
  fd.list[[j]] <- fd.object.array
  fpca.list2[[j]] <- pca.fd(fd.object.array, nharm=M)
  print(sum(fpca.list2[[j]]$varprop))
}

## collect scores
fpca.score <- numeric()
for(j in 1:p){
  fpca.score <- cbind(fpca.score, inprod(fd.list[[j]], fpca.list2[[j]]$harmonics))
}
fpca.score.cen <- scale(fpca.score, center=T, scale=F)

pdf("~/work/CMC-GGM/fgm/test.pdf")
plot(fpca.score.cen[,1], fpca.score.cen[,2])
dev.off()

pdf("~/work/CMC-GGM/fgm/raw_test.pdf")
plot(delta[,1], delta[,2])
dev.off()

pdf("~/work/CMC-GGM/fgm/signal_test.pdf")
plot(obs.time, h[1,1,], type = "l")
lines(obs.time, h[2,1,])
dev.off()
# ############# using common basis across different nodes
# # produce splines functions defined on interval [0,1]
# data.aggregated <- numeric()
# for (k in 1:p){
#   data.aggregated <- rbind(data.aggregated,(h[,k,]))
# }
# 
# # estimate the functions using splines
# bspline.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=M)
# data.fd <- Data2fd(argvals=obs.time, y=t(data.aggregated), basisobj=bspline.basis)
# 
# # compute the  first M functional principal components
# fpca.object <- pca.fd(data.fd, nharm=M)
# fpca.score.agg <- inprod(data.fd, fpca.object$harmonics)
# fpca.score <- matrix(0, nrow=n, ncol=M*p)
# 
# for (k in 1:p){
#   fpca.score[(1:n),((k-1)*M+1):(k*M)] <- fpca.score.agg[((k-1)*n+1):(k*n),]
# }
# 
# print(sum(fpca.object$varprop))
  
####################################
#             PART 3:              #
#    THRESHOLDING AND FGLASSO      #
#        USING FPCA BASIS          #
#        AND GET ROC CURVE         #
####################################

U.df = read.csv(file=paste(read.path,"/proj_subspace_model1_exp_d05.csv",sep=""), header = FALSE)

# three transformations
data.pcmcopula = mult.trans(X=fpca.score.cen, p = p, fun = 'pcmcopula', U.list = U.list)
data.cmcopula = mult.trans(X=fpca.score.cen, p = p, fun = 'cmcopula')
data.copula = mult.trans(X=fpca.score.cen, p = p, fun = 'copula')

pdf("~/work/CMC-GGM/fgm/test_copula.pdf")
plot(data.copula[,1], data.copula[,2])
dev.off() 


# covariance matrices
S.cmcg <- cor(data.cmcg)
S.copula <- cor(data.copula)
S.linear <- cor(data.linear)

## thresholding
est.prec.cmcg.thres <- vggm.thr(X=S.cmcg, p = p, rho=10^-2, ridge = TRUE)$Theta0
est.prec.copula.thres <- vggm.thr(X=S.copula, p = p, rho=10^-2, ridge = TRUE)$Theta0
est.prec.linear.thres <- vggm.thr(X=S.linear, p = p, rho=10^-2, ridge = TRUE)$Theta0

roc.cmcg.thres <- roc_curve(est.prec.cmcg.thres, Omega0, n_grid = L)
roc.copula.thres <- roc_curve(est.prec.copula.thres, Omega0, n_grid = L)
roc.linear.thres <- roc_curve(est.prec.linear.thres, Omega0, n_grid = L)

###################
## fglasso
#############
# tuning parameter

lambda.cmcg.max <- max(S.cmcg)
lambda.copula.max <- max(S.copula)
lambda.linear.max = max(S.linear)
lambda.cmcg = exp(seq(log(lambda.cmcg.max), log(5*10^-2*lambda.cmcg.max), length = L))
lambda.copula = exp(seq(log(lambda.copula.max), log(5*10^-2*lambda.copula.max), length = L))
lambda.linear = exp(seq(log(lambda.linear.max), log(5*10^-2*lambda.linear.max), length = L))

###############
## first eigenvalue
rho.cmcg = eigen(S.cmcg)$values[1]
rho.copula = eigen(S.copula)$values[1]
rho.linear = eigen(S.linear)$values[1]

library(doParallel)
cl <- makeCluster(L)
registerDoParallel(cl)

tol.Theta <-1e-2
tol.beta <-1e-3
maxit.Theta <-20
maxit.beta <-5
tol.beta_k <-1e-5
maxit.beta_k <-500
method <- 1

roc.cmcg.fglasso <- foreach(l = 1:L, .combine="rbind", .packages=c("matrixcalc","igraph","BB")) %dopar% {
  G.qiao <- ProxAlg_FGM(S.cmcg, p=p, M=M, lambda.cmcg[l])$Support # p*p support
  # G.qiao <- bglasso(S = S.cmcg+10^-3*rho.cmcg*diag(M*p), lambda = lambda.cmcg[l], n = n, p = p,
  #                   tol.Theta=tol.Theta, tol.beta=tol.beta,
  #                   tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
  #                   maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.qiao.l <- roc(est = G.qiao, pop = Omega0, thres = 0)
  TFPR.qiao.l
}

roc.copula.fglasso <- foreach(l = 1:L, .combine="rbind", .packages="matrixcalc") %dopar% {
  G.qiao <- ProxAlg_FGM(S.copula, p=p, M=M, lambda.copula[l])$Support # p*p support
  # G.qiao <- bglasso(S = S.copula+10^-3*rho.copula*diag(M*p), lambda = lambda.copula[l], n = n, p = p,
  #                   tol.Theta=tol.Theta, tol.beta=tol.beta,
  #                   tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
  #                   maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.qiao.l <- roc(est = G.qiao, pop = Omega0, thres = 0)
  TFPR.qiao.l
}

roc.linear.fglasso <- foreach(l = 1:L, .combine="rbind", .packages="matrixcalc") %dopar% {
  G.qiao <- ProxAlg_FGM(S.linear, p=p, M=M, lambda.linear[l])$Support # p*p support
  # G.qiao <- bglasso(S = S.linear+10^-3*rho.linear*diag(M*p), lambda = lambda.linear[l], n = n, p = p,
  #                   tol.Theta=tol.Theta, tol.beta=tol.beta,
  #                   tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
  #                   maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
  TFPR.qiao.l <- roc(est = G.qiao, pop = Omega0, thres = 0)
  TFPR.qiao.l
}

stopCluster(cl)

####################################
#      PART 4: SAVE ROC DATA       #
####################################
roc.cmcg.fglasso = t(roc.cmcg.fglasso)
roc.copula.fglasso = t(roc.copula.fglasso)
roc.linear.fglasso = t(roc.linear.fglasso)

roc.df = data.frame(trans.methods = rep(c("cmcg","copula",'linear',
                                          "cmcg","copula",'linear'),each=dim(roc.cmcg.thres)[2]),
                    methods = rep(c("thresholding","fglasso"),each=3*dim(roc.cmcg.thres)[2]),
                    TP = cbind(roc.cmcg.thres, roc.copula.thres, roc.linear.thres, 
                               roc.cmcg.fglasso, roc.copula.fglasso, roc.linear.fglasso)[1,],
                    FP = cbind(roc.cmcg.thres, roc.copula.thres, roc.linear.thres,
                               roc.cmcg.fglasso, roc.copula.fglasso, roc.linear.fglasso)[2,])
roc.df$TP[roc.df$FP == 0] = 0
# plot<-ggplot(roc.df, aes(x=FP, y=TP)) +
#   geom_line(aes(color=trans.methods,linetype=methods))
# print(plot)
#########################################################
time.end <- proc.time()
time.run <- (time.end - time.start)[3]
print(time.run)

save(roc.df, file=paste(save.path,"/ROC.A.seed",args,"M", M, ".Rdata",sep=""))  

