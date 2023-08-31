############### PART 0: library, source and pathes ########################

library(dplyr)
library(pdist)
library(pracma)
library(clue)
library(BB)
library(mvtnorm)
library(doParallel)
library(randtoolbox)


reading.path <- "~/work/CMC-GGM/EEG/eeg_processed/"
reading.path.2 <- "~/work/CMC-GGM/EEG/TuneResult/"
working.path <- "~/work/CMC-GGM/EEG/"
function.path <- "~/work/CMC-GGM/Functions/"
save.path <- "~/work/CMC-GGM/EEG/Result/"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

rand.seed <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(rand.seed)

############### PART 1: reading data #######################

n.control <- 45
n.alcohol <- 77
p <- 64
tau <- 256
M <- 6
bt.index.alc <- sample(1:n.alcohol, replace = TRUE)
bt.index.ctr <- sample(1:n.control, replace = TRUE)

fpca.score.alcohol = read.csv(file=paste(reading.path,"alcohol_score.csv",sep=""), header = TRUE, row.names = 1)
fpca.score.control = read.csv(file=paste(reading.path,"control_score.csv",sep=""), header = TRUE, row.names = 1)
fpca.score.alcohol = fpca.score.alcohol[bt.index.alc,]
fpca.score.control = fpca.score.control[bt.index.ctr,]

U.alcohol.01 = read.csv(file=paste(reading.path,"/proj_subspace_alcohol_01.csv",sep=""), header = FALSE)
U.control.01 = read.csv(file=paste(reading.path,"/proj_subspace_control_01.csv",sep=""), header = FALSE)
U.alcohol.02 = read.csv(file=paste(reading.path,"/proj_subspace_alcohol_02.csv",sep=""), header = FALSE)
U.control.02 = read.csv(file=paste(reading.path,"/proj_subspace_control_02.csv",sep=""), header = FALSE)

getUlist <- function(U, p, M, k){
  U.list <- list()
  for (i in 1:p) {
    U.list[[i]] = matrix(as.numeric(U[i,]), nrow = M, ncol = k)
  }
  return(U.list)
}

U.alcohol.01.list <- getUlist(U = U.alcohol.01, p = p, M = M, k = 1)
U.alcohol.02.list <- getUlist(U = U.alcohol.02, p = p, M = M, k = 2)
U.control.01.list <- getUlist(U = U.control.01, p = p, M = M, k = 1)
U.control.02.list <- getUlist(U = U.control.02, p = p, M = M, k = 2)
lambda_list <- seq(1, 2, 0.1)

## transformed score
pcmc.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'pcmcopula', U.list = U.alcohol.01.list)
cmc.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'cmcopula', grid.method = 'halton')
copula.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'copula')

pcmc.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'pcmcopula', U.list = U.control.01.list)
cmc.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'cmcopula', grid.method = 'halton')
copula.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'copula')


# reading tuning parameters
lambda_list_pcmc <- seq(2.3, 2.7, len=20)
lambda_list_cmc <- seq(1.3, 1.7, len=20)
lambda_list_copula <- seq(2, 2.4, len=20)
sparsity.mat <- numeric()
for (i in 1:20) {
  load(paste(reading.path.2,"sparsity", i,".Rdata", sep = ""))
  sparsity.mat = rbind(sparsity.mat, sparsity.list)
}
lambda_index = apply(abs(sparsity.mat-0.05), 2, which.min)
lambda = c(lambda_list_pcmc[lambda_index[1]],  lambda_list_cmc[lambda_index[2]], lambda_list_copula[lambda_index[3]],
           lambda_list_pcmc[lambda_index[4]], lambda_list_cmc[lambda_index[5]], lambda_list_copula[lambda_index[6]])

################## PART 2: Graph Estimation ########################
graph.est <- function(score, p, lambda){
  ## score: fpca score vectors with dimension n by (p*M)
  ## method: {'glasso','thresholding'}
  ## trans: {'cmcopula', 'copula', 'gaussian'}
  ## L: number of lambdas
  M = dim(score)[2] / p
  n = dim(score)[1]
  S = cor(score)
  
  bglasso.list <- bglasso(S = S+10^-3*diag(M*p), lambda = lambda, n = n, p = p)
  est.prec <- bglasso.list$alpha
  est.Theta <- bglasso.list$Theta
  sparsity.graph <- sum(est.prec > 0)/(p*(p-1))
  
  
  return(list(est.prec=est.prec, sparsity=sparsity.graph))
}

time_start <- proc.time()
alc.pcmc <- graph.est(score=pcmc.score.alcohol, p=p, lambda = lambda[1])
alc.cmc <- graph.est(score=cmc.score.alcohol, p=p, lambda = lambda[2])
alc.copula <- graph.est(score=copula.score.alcohol, p=p, lambda = lambda[3])

ctr.pcmc <- graph.est(score=pcmc.score.control, p=p, lambda = lambda[4])
ctr.cmc <- graph.est(score=cmc.score.control, p=p, lambda = lambda[5])
ctr.copula <- graph.est(score=copula.score.control, p=p, lambda = lambda[6])

time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

################### PART 3: save data ############################
alc.pcmc.prec <- alc.pcmc$est.prec
alc.cmc.prec <- alc.cmc$est.prec
alc.copula.prec <- alc.copula$est.prec

ctr.pcmc.prec <- ctr.pcmc$est.prec
ctr.cmc.prec <- ctr.cmc$est.prec
ctr.copula.prec <- ctr.copula$est.prec

prec.list <- list(alc.pcmc.prec, alc.cmc.prec, alc.copula.prec,
                  ctr.pcmc.prec, ctr.cmc.prec, ctr.copula.prec)

save(prec.list, file=paste(save.path,"/prec.bt.seed",rand.seed,".Rdata",sep=""))

alc.cmc.sparsity <- alc.cmc$sparsity
alc.copula.sparsity <- alc.copula$sparsity
alc.pcmc.sparsity <- alc.pcmc$sparsity

ctr.cmc.sparsity <- ctr.cmc$sparsity
ctr.copula.sparsity <- ctr.copula$sparsity
ctr.pcmc.sparsity <- ctr.pcmc$sparsity

cat(paste("pcmc-alc sparsity:", alc.pcmc.sparsity,
          ";cmc-alc sparsity:",alc.cmc.sparsity,
          ";copula-alc sparsity:", alc.copula.sparsity))

cat(paste("pcmc-ctr sparsity:", ctr.pcmc.sparsity,
          ";cmc-ctr sparsity:",ctr.cmc.sparsity,
          ";copula-ctr sparsity:", ctr.copula.sparsity))
