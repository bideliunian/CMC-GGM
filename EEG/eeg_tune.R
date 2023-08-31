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
working.path <- "~/work/CMC-GGM/EEG/"
function.path <- "~/work/CMC-GGM/Functions/"
save.path <- "~/work/CMC-GGM/EEG/TuneResult"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)


arg <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(arg)


################## PART 1: reading data##########################

n.control <- 45
n.alcohol <- 77
p <- 64
tau <- 256
M <- 6

fpca.score.alcohol = read.csv(file=paste(reading.path,"alcohol_score.csv",sep=""), header = TRUE, row.names = 1)
fpca.score.control = read.csv(file=paste(reading.path,"control_score.csv",sep=""), header = TRUE, row.names = 1)

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
lambda_list_pcmc <- seq(2.3, 2.7, len=20)
lambda_list_cmc <- seq(1.3, 1.7, len=20)
lambda_list_copula <- seq(2, 2.4, len=20)
lambda_pcmc <- lambda_list_pcmc[arg]
lambda_cmc <- lambda_list_cmc[arg]
lambda_copula <- lambda_list_copula[arg]


## transformed score
pcmc.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'pcmcopula', U.list = U.alcohol.01.list)
cmc.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'cmcopula', grid.method = 'halton')
copula.score.alcohol = mult.trans(X=fpca.score.alcohol, p = p, fun = 'copula')

pcmc.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'pcmcopula', U.list = U.control.01.list)
cmc.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'cmcopula', grid.method = 'halton')
copula.score.control = mult.trans(X=fpca.score.control, p = p, fun = 'copula')


################## PART 2: Tuning parameter based on sparsity ########################
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
alc.pcmc <- graph.est(score=pcmc.score.alcohol, p=p, lambda = lambda_pcmc)
alc.cmc <- graph.est(score=cmc.score.alcohol, p=p, lambda = lambda_cmc)
alc.copula <- graph.est(score=copula.score.alcohol, p=p, lambda = lambda_copula)

ctr.pcmc <- graph.est(score=pcmc.score.control, p=p, lambda = lambda_pcmc)
ctr.cmc <- graph.est(score=cmc.score.control, p=p, lambda = lambda_cmc)
ctr.copula <- graph.est(score=copula.score.control, p=p, lambda = lambda_copula)

time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)



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

sparsity.list = c(alc.pcmc.sparsity, alc.cmc.sparsity, alc.copula.sparsity, 
                  ctr.pcmc.sparsity, ctr.cmc.sparsity, ctr.copula.sparsity)

save(sparsity.list, file=paste(save.path,"/sparsity",arg,".Rdata",sep=""))
