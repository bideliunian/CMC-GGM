#######################
## Gene/Protein Regulatory Network Inference
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Gene_Protein"
save.path <- "~/work/CMC-GGM/Gene_Protein/CVResult"

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

arg <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(arg)
############### PART 1: reading data #######################
load(paste(working.path, "/gpconsensus.RData", sep = ""))
gene <-  scale(t(dat$gene.select), TRUE,TRUE)
prot <-  scale(t(dat$protein.select), TRUE,TRUE)
n <- nrow(gene)
p <- ncol(gene)
M <- 2
lambda_list <- seq(0.2, 0.65, 0.05)
lambda <- lambda_list[arg]
n.fold <- 10
###### bootstrap ###############


data.combined <- matrix(0, n, p*M)
for (i in 1:p) {
  data.combined[, 2*i-1] <- gene[, i]
  data.combined[, 2*i] <- prot[, i]
}

cmc.score = mult.trans(X=data.combined, p = p, fun = 'cmcopula')
copula.score = mult.trans(X=data.combined, p = p, fun = 'copula')
linear.score = mult.trans(X=data.combined, p = p, fun = 'scale')

################## PART 2: CV for parameter ########################

cverror.cmc <- CV4tuning(cmc.score, p, lambda, n.fold)
cverror.copula <- CV4tuning(copula.score, p, lambda, n.fold)
cverror.linear <- CV4tuning(linear.score, p, lambda, n.fold)
cverror.list <- c(cverror.cmc, cverror.copula, cverror.linear)

save(cverror.list, file=paste(save.path,"/cverror",arg,".Rdata",sep=""))
