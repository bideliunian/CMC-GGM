#######################
## Model_b ROC (Model 2)
## p = 100, n = 200
## copula
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_b_nblasso"
save.path <- "~/work/CMC-GGM/Model_b_nblasso/Results"

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
n = 300
p = 100
M = 2
L = 20 ## number of grids in roc curve
model = 'model2'
trans.type = 'copula'
n.expr = 50

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
  
  data.cmcopula = mult.trans(X=data.X, p = p, fun = 'cmcopula')
  data.copula = mult.trans(X=data.X, p = p, fun = 'copula')
  data.linear = mult.trans(X=data.X, p = p, fun = 'scale')
  
  roc.cmc <- get.roc.nblasso(data.cmcopula, p, G, L)
  roc.copula <- get.roc.nblasso(data.copula, p, G, L)
  roc.linear <- get.roc.nblasso(data.linear, p, G, L)
  
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

save(df.list, file=paste(save.path,"/ROC.b.100.copula.nbd.Rdata",sep=""))