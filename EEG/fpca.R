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
func.path <- "~/work/CMC-GGM/Functions/"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

############## PART 1: Loading data #######################

load(file = paste(reading.path, 'eeg_control.Rdata', sep = ''))
load(file = paste(reading.path, 'eeg_alcohol.Rdata', sep = ''))
load(file = paste(reading.path, 'nodes_name.Rdata', sep = ''))

# remove ears and nose
# idx <- c(which(acnames=="X"),which(acnames=="Y"),which(acnames=="nd"))
# eegdata.control.mean <- eegdata.control.mean[,-idx,]
# eegdata.alcohol.mean <- eegdata.alcohol.mean[,-idx,]

n.control <- dim(eegdata.control.array)[1]
n.alcohol <- dim(eegdata.alcohol.array)[1]
p <- dim(eegdata.control.array)[2]
tau <- dim(eegdata.control.array)[3]

##### functional pca ########

fpca.eeg <- function(h){
  # Number of basis functions
  n <- dim(h)[1]
  p <- dim(h)[2]
  tau <- dim(h)[3]
  
  obs.time <- seq(1/tau, 1, 1/tau) # vector of observation time points of delta
  fpca.list <- list()
  
  for(j in 1:p){
    fpca.object <- MakeFPCAInputs(IDs = rep(1:n, each=tau), tVec=rep(obs.time,n), t(h[,j,]))
    fpca.result <- FPCA(fpca.object$Ly, fpca.object$Lt)
    fpca.list[[j]] <- fpca.result
  }
  
  return(fpca.list)
}

################
fpca.control <- fpca.eeg(h=eegdata.control.array)
fpca.alcohol <- fpca.eeg(h=eegdata.alcohol.array)

## select M to explain 90% percentage of variance
M.list.control <- M.list.alcohol <- numeric()
for(j in 1:p){
  M.list.control <- c(M.list.control, SelectK(fpca.control[[j]], criterion = 'FVE', FVEthreshold = 0.9)$K) 
  M.list.alcohol <- c(M.list.alcohol, SelectK(fpca.alcohol[[j]], criterion = 'FVE', FVEthreshold = 0.9)$K) 
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

M.control <- getmode(M.list.control)
M.alcohol <- getmode(M.list.alcohol)
print(M.control)
print(M.alcohol)

## calculating the fpca scores
fpca.score.control <- fpca.score.alcohol <- numeric()
for(j in 1:p){
  fpca.score.control <- cbind(fpca.score.control, fpca.control[[j]]$xiEst[,1:M.control]) 
  fpca.score.alcohol <- cbind(fpca.score.alcohol, fpca.alcohol[[j]]$xiEst[,1:M.alcohol])  
}

write.csv(fpca.score.alcohol, file=paste(reading.path, 'alcohol_score.csv',sep=""))
write.csv(fpca.score.control, file=paste(reading.path, "control_score.csv",sep=""))

