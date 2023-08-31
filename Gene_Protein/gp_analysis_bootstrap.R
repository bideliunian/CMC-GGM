#######################
## Gene/Protein Regulatory Network Inference
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Gene_Protein"
save.path <- "~/work/CMC-GGM/Gene_Protein/Result"
reading.path <- "~/work/CMC-GGM/Gene_Protein/CVResult"

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

rand.seed <- as.numeric(commandArgs(trailingOnly=TRUE))
set.seed(rand.seed)

############### PART 1: reading data #######################

load(paste(working.path, "/gpconsensus.RData", sep = ""))
gene <-  scale(t(dat$gene.select), TRUE,TRUE)
prot <-  scale(t(dat$protein.select), TRUE,TRUE)
n <- nrow(gene)
p <- ncol(gene)
M <- 2
bt.index <- sample(1:n, replace = TRUE)
###### bootstrap ###############

data.combined <- matrix(0, n, p*M)
for (i in 1:p) {
  data.combined[, 2*i-1] <- gene[bt.index, i]
  data.combined[, 2*i] <- prot[bt.index, i]
}

cmc.score = mult.trans(X=data.combined, p = p, fun = 'cmcopula')
copula.score = mult.trans(X=data.combined, p = p, fun = 'copula')
linear.score = mult.trans(X=data.combined, p = p, fun = 'scale')

S.cmc = cor(cmc.score)
S.copula = cor(copula.score)
S.linear = cor(linear.score)

# reading tuning parameters
lambda_list <- seq(0.2, 0.65, 0.05)
n.lambda <- length(lambda_list)
cverror.mat <- numeric()
for (i in 1:n.lambda) {
  load(paste(reading.path,"/cverror", i,".Rdata", sep = ""))
  cverror.mat = rbind(cverror.mat, cverror.list)
}
lambda = lambda_list[apply(cverror.mat, 2, which.min)]

################## PART 2: Graph Estimation ########################

graph.est <- function(S, n, p, lambda=NULL, sparsity = NULL){
  ## score: fpca score vectors with dimension n by (p*M)
  ## method: {'glasso','thresholding'}
  ## trans: {'cmcopula', 'copula', 'gaussian'}
  ## L: number of lambdas
  M = dim(S)[2] / p
  
  tol.Theta <-1e-2
  tol.beta <-1e-3
  maxit.Theta <-10
  maxit.beta <-10
  tol.beta_k <-1e-7
  maxit.beta_k <-1500
  method <- 1
  L <-  10
  
  if(is.null(lambda)){
    lambda = seq(0.3, 1, length = L)
    
    cl <- makeCluster(L)
    registerDoParallel(cl)
    
    result <- foreach(l = 1:L, .combine="c", .packages=c("matrixcalc","igraph","BB")) %dopar% {
      func.path <- "~/work/CMC-GGM/Functions"
      source(paste(func.path,"ns_lasso.R", sep="/"))
      bglasso.list <- bglasso(S = S+10^-3*diag(M*p), lambda = lambda[l], n = n, p = p,
                              tol.Theta=tol.Theta, tol.beta=tol.beta,
                              tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                              maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)
      list(bglasso.list)
    }
    
    stopCluster(cl)
    # result <- list()
    # for(l in 1:length(lambda)){
    #   result[[l]] <- bglasso(S = S, lambda = lambda[l], n = n, p = p,
    #                           tol.Theta=tol.Theta, tol.beta=tol.beta,
    #                           tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
    #                           maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)
    # }
    
    bic.list <- lapply(result, "[[", 'bic')
    aic.list <- lapply(result, "[[", 'aic')
    est.prec.list <- lapply(result, "[[", 'alpha')
    Theta.list <- lapply(result, "[[", 'Theta')
    iter.list <- lapply(result, "[[", 'iter.Theta')
    print(paste("aic:", as.vector(unlist(aic.list))))
    print(paste("bic:", as.vector(unlist(bic.list))))
    print(paste("iter:", as.vector(unlist(iter.list))))
    
    if(!is.null(sparsity)){
      sparsity.list <- lapply(est.prec.list, function(x) sum(x > 0)/(p*(p-1)))
      sparsity.list <- as.vector(unlist(sparsity.list))
      print(sparsity.list)
      index <- which.min(abs(sparsity.list - sparsity))
      sparsity.graph <- sparsity.list[index]
      est.prec <- est.prec.list[[index]]
    }
    else{
      bic.list <- as.vector(unlist(bic.list))
      index <- which.min(bic.list)
      est.prec <- est.prec.list[[index]]
      sparsity.graph <- sum(est.prec > 0)/(p*(p-1))
    }
    
    est.Theta <- Theta.list[[index]] 
  }
  
  else {
    bglasso.list <- bglasso(S = S+10^-3*diag(M*p), lambda = lambda, n = n, p = p,
                            tol.Theta=tol.Theta, tol.beta=tol.beta,
                            tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                            maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)
    est.prec <- bglasso.list$alpha
    est.Theta <- bglasso.list$Theta
    sparsity.graph <- sum(est.prec > 0)/(p*(p-1))
  }
  return(list(est.prec=est.prec, sparsity=sparsity.graph, Theta = est.Theta))
}

time_start <- proc.time()
est.cmc.glasso <- graph.est(S=S.cmc, n=n, p=p, lambda = lambda[1])
est.copula.glasso <- graph.est(S=S.copula, n=n, p=p, lambda = lambda[2])
est.linear.glasso <- graph.est(S=S.linear, n=n, p=p, lambda = lambda[3])
time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

########
# library(energy)
# 
# p.values.cmc <- p.values.copula <- p.values.linear <- numeric()
# 
# for (i in 1:p) {
#   index.i <- ((i-1)*M+1):(i*M)
#   p.values.cmc <- c(p.values.cmc, mvnorm.etest(cmc.score[,index.i], R=500)$p.value)
#   p.values.copula <- c( p.values.copula, mvnorm.etest(copula.score[,index.i], R=500)$p.value)
#   p.values.linear <- c(p.values.linear, mvnorm.etest(linear.score[,index.i], R=500)$p.value)
# }
# 
# df.pvalues.linear <- data.frame(pvalues = p.values.linear) 
# ## histogram of the p.values after copula transformation
# 
# p.linear <- ggplot(df.pvalues.linear, aes(x=pvalues)) + geom_histogram(breaks=seq(0,0.8,by=0.025)) + 
#   geom_vline(aes(xintercept=0.05, color='red'), linetype='dashed', size=1) + 
#   theme(legend.position="none", text = element_text(size = 20))


################### PART 3: save data ############################
cmc.prec <- est.cmc.glasso$est.prec
copula.prec <- est.copula.glasso$est.prec
linear.prec <- est.linear.glasso$est.prec
prec.list <- list(cmc.prec, copula.prec, linear.prec)

save(prec.list, file=paste(save.path,"/prec.bt.seed",rand.seed,".Rdata",sep=""))

cmc.sparsity <- est.cmc.glasso$sparsity
copula.sparsity <- est.copula.glasso$sparsity
linear.sparsity <- est.linear.glasso$sparsity

cat(paste("cmc-ggm sparsity:", cmc.sparsity,";copula-ggm sparsity:",copula.sparsity,
          ";ggm sparsity:", linear.sparsity))
