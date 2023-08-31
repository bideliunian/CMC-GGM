############### PART 0: library, source and pathes ########################

library(dplyr)
library(tidyr)
library(fda)
library(fdapace)
library(pdist)
library(pracma)
library(clue)
library(BB)
library(mvtnorm)
library(doParallel)
library(randtoolbox)
library(fdapace)
library(igraph)
library(ggplot2)


reading.path <- "~/work/CMC-GGM/EEG/eeg_processed/"
working.path <- "~/work/CMC-GGM/EEG/"
function.path <- "~/work/CMC-GGM/Functions/"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)


##################PART 1: reading data##########################

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


################## PART 2: Graph Estimation ########################

fgm.eeg <- function(score, p, trans, U.list=NULL, lambda){
  ## score: fpca score vectors with dimension n by (p*M)
  ## method: {'glasso','thresholding'}
  ## trans: {'cmcopula', 'copula', 'gaussian'}
  ## L: number of lambdas
  M = dim(score)[2] / p
  n = dim(score)[1]
  
  if(trans == 'pcmcopula'){
    score.trans = mult.trans(X=score, p = p, fun = trans, U.list = U.list)
  }
  else{
    score.trans = mult.trans(X=score, p = p, fun = trans, grid.method = 'halton')
  }
  
  S = cor(score.trans)
  
  bglasso.list <- bglasso(S = S+10^-3*diag(M*p), lambda = lambda, n = n, p = p)
  est.prec <- bglasso.list$alpha
  est.Theta <- bglasso.list$Theta
  sparsity.graph <- sum(est.prec > 0)/(p*(p-1))

  
  return(list(est.prec=est.prec, score.trans=score.trans, sparsity=sparsity.graph))
}


## control group glasso
est.pcmc.control = fgm.eeg(score=fpca.score.control, p=p, trans='pcmcopula', lambda = 1.4, 
                                  U.list = U.control.01.list)
est.cmc.control = fgm.eeg(score=fpca.score.control, p=p, trans='cmcopula',lambda = 1)
est.copula.control = fgm.eeg(score=fpca.score.control, p=p, trans='copula', lambda = 1.6)

## alcohol group glasso
est.pcmc.alcohol = fgm.eeg(score=fpca.score.control, p=p, trans='pcmcopula', lambda = 1.4, 
                           U.list = U.alcohol.01.list)
est.cmc.alcohol = fgm.eeg(score=fpca.score.alcohol, p=p, trans='cmcopula', lambda = 1)
est.copula.alcohol = fgm.eeg(score=fpca.score.alcohol, p=p, trans='copula', lambda = 1.6)

## sparsity
print(paste("pcmc alcohol sparsity:", est.pcmc.alcohol$sparsity))
print(paste("pcmc control sparsity:", est.pcmc.control$sparsity))
print(paste("cmc alcohol sparsity:", est.cmc.alcohol$sparsity))
print(paste("cmc control sparsity:", est.cmc.control$sparsity))
print(paste("copula alcohol sparsity:", est.copula.alcohol$sparsity))
print(paste("copula control sparsity:", est.copula.control$sparsity))


################################ PART 3: Draw network #####################
node.names <- channel.list
position.list <- list("FPZ"=c(0, 0.93), "AFZ"=c(0, 0.6), "FZ"=c(0, 0.4), FCZ=c(0, 0.2), 
                      "CZ"=c(0, 0), "CPZ"=c(0, -0.2), "PZ"=c(0, -0.4), "POZ"=c(0, -0.6), 
                      "OZ"=c(0, -0.8), "nd"=c(0,-1), "C2"=c(0.2, 0), "C4"=c(0.4, 0), 
                      "C6"=c(0.6, 0), "T8"=c(0.8, 0), "Y"=c(1, 0), "C1"=c(-0.2, 0), 
                      "C3"=c(-0.4, 0), "C5"=c(-0.6, 0), "T7"=c(-0.8, 0), "X"=c(-1, 0), 
                      "FP2"=c(0.3, 0.85), "FP1"=c(-0.3, 0.85), "AF2"=c(0.26, 0.62), 
                      "AF1"=c(-0.26, 0.62), "AF8"=c(0.55, 0.68), "AF7"=c(-0.55, 0.68), 
                      "F2"=c(0.21, 0.41), "F1"=c(-0.21, 0.41), "F4"=c(0.4, 0.45), 
                      "F3"=c(-0.4, 0.45), "F6"=c(0.55, 0.5), "F5"=c(-0.55, 0.5), 
                      "F8"=c(0.75, 0.55), "F7"=c(-0.75, 0.55), "FC2"=c(0.25, 0.21), 
                      "FC1"=c(-0.25, 0.21), "FC4"=c(0.5, 0.22), "FC3"=c(-0.5, 0.22), 
                      "FC6"=c(0.7, 0.26), "FC5"=c(-0.7, 0.26), "FT8"=c(0.9, 0.31), 
                      "FT7"=c(-0.9, 0.31), "CP2"=c(0.25, -0.21), "CP4"=c(0.5, -0.22), 
                      "CP6"=c(0.7, -0.26), "TP8"=c(0.9, -0.31), "CP1"=c(-0.25, -0.21), 
                      "CP3"=c(-0.5, -0.22), "CP5"=c(-0.7, -0.26), "TP7"=c(-0.9, -0.31), 
                      "P2"=c(0.21, -0.41), "P4"=c(0.4, -0.45), "P6"=c(0.55, -0.5), 
                      "P8"=c(0.65, -0.55), "P1"=c(-0.21, -0.41), "P3"=c(-0.4, -0.45), 
                      "P5"=c(-0.55, -0.5), "P7"=c(-0.65, -0.55), "PO2"=c(0.2, -0.62), 
                      "PO8"=c(0.45, -0.68), "PO1"=c(-0.2, -0.62), "PO7"=c(-0.45, -0.68), 
                      "O2"=c(0.2, -0.85), "O1"=c(-0.2, -0.85))

name.v <- c()
layMat <- matrix(NA, nrow=64, ncol=2)
for (i in 1:length(node.names)){
  x <- node.names[i]
  layMat[i, 1] <- unlist(position.list)[paste(x, 1, sep="")]
  layMat[i, 2] <- unlist(position.list)[paste(x, 2, sep="")]
}
## ploting the connected graph from control group and alcohol group together

## alcohol group and glasso
pdf(paste("~/work/CMC-GGM/EEG/graph_cmc.pdf"))
network.plotting(est.cmc.alcohol$est.prec, est.cmc.control$est.prec, node.names = node.names)
dev.off()
pdf(paste("~/work/CMC-GGM/EEG/graph_copula.pdf"))
network.plotting(est.copula.alcohol$est.prec, est.copula.control$est.prec, node.names = node.names)
dev.off
