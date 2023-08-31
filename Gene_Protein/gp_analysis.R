#######################
## Gene/Protein Regulatory Network Inference
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Gene_Protein"
save.path <- "~/work/CMC-GGM/Gene_Protein"

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
library(ggplot2)
library(RColorBrewer)

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

############### PART 1: reading data #######################

load(paste(working.path, "/gpconsensus.RData", sep = ""))
gene <-  scale(t(dat$gene.select), TRUE,TRUE)
prot <-  scale(t(dat$protein.select), TRUE,TRUE)
n <- nrow(gene)
p <- ncol(gene)
M <- 2

data.combined <- matrix(0, n, p*M)
for (i in 1:p) {
  data.combined[, 2*i-1] <- gene[, i]
  data.combined[, 2*i] <- prot[, i]
}

################## PART 2: Graph Estimation ########################

graph.est <- function(score, p, method, trans, sparsity = NULL, bic = FALSE){
  ## score: fpca score vectors with dimension n by (p*M)
  ## method: {'glasso','thresholding'}
  ## trans: {'cmcopula', 'copula', 'gaussian'}
  ## L: number of lambdas
  M = dim(score)[2] / p
  n = dim(score)[1]
  if(is.null(sparsity)) {bic == TRUE}
  
  score.trans = mult.trans(X=score, p = p, fun = trans)
  S = cor(score.trans)
  
  if(method=='glasso'){
    tol.Theta <-1e-2
    tol.beta <-1e-3
    maxit.Theta <-5
    maxit.beta <-5
    tol.beta_k <-1e-7
    maxit.beta_k <-500
    method <- 1
    L = 10
    
    lambda = seq(0.5, 2, length = L)

    rho = eigen(S)$values[1]
    
    cl <- makeCluster(L)
    registerDoParallel(cl)
    
    est.prec.list <- foreach(l = 1:L, .combine="c", .packages=c("matrixcalc","igraph","BB")) %dopar% {
      func.path <- "~/work/CMC-GGM/Functions"
      source(paste(func.path,"ns_lasso.R", sep="/"))
      G.bglasso <- bglasso(S = S+10^-3*rho*diag(M*p), lambda = lambda[l], n = n, p = p,
                           tol.Theta=tol.Theta, tol.beta=tol.beta,
                           tol.beta_k=tol.beta_k, maxit.Theta=maxit.Theta, maxit.beta=maxit.beta,
                           maxit.beta_k=maxit.beta_k, Theta=NULL, Theta.inv= NULL, method=method)$alpha
      list(G.bglasso)
    }
    
    stopCluster(cl)
    
    if(!bic){
      sparsity.list <- lapply(est.prec.list, function(x) sum(x[upper.tri(x)] > 0)/(p*(p-1)/2))
      sparsity.list <- as.vector(sparsity.list)
      sparsity.list[sparsity.list>sparsity] <- 0
      index <- which.max(sparsity.list)
      sparsity.graph <- sparsity.list[index]
      est.prec <- est.prec.list[[index]] 
    }
    else{
      bic.list <- lapply(est.prec.list, bic, S=S, M=M, n=n)
      print(bic.list)
      bic.list <- as.vector(bic.list)
      index <- which.min(bic.list)
      est.prec <- est.prec.list[[index]]
      sparsity.graph <- sum(est.prec[upper.tri(est.prec)] > 0)/(p*(p-1)/2)
    }
  }
  
  if(method == 'thresholding'){
    est.prec <- vggm.thr(X=S, p=p, rho=10^-3, ridge=TRUE)$Theta0
    thres <- quantile(x=est.prec[upper.tri(est.prec)], probs=1-sparsity)
    est.prec[est.prec<thres] <- 0
    sparsity.graph <- sparsity
  }
  
  return(list(est.prec=est.prec, score.trans=score.trans, sparsity=sparsity.graph))
}

## control group/ thresolding
est.cmc.glasso <- graph.est(score=data.combined, p=p, method='glasso', 
                           trans='cmcopula')
est.copula.glasso <- graph.est(score=data.combined, p=p, method='glasso', 
                              trans='copula')
est.linear.glasso <- graph.est(score=data.combined, p=p, method='glasso', 
                              trans='scale')

cmc.score <- est.cmc.glasso$score.trans
copula.score <- est.copula.glasso$score.trans
linear.score <- est.linear.glasso$score.trans

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


################### PART 3: Visualization ############################
cmc.prec <- est.cmc.glasso$est.prec
copula.prec <- est.copula.glasso$est.prec
linear.prec <- est.linear.glasso$est.prec
prec.list <- list(cmc.prec, copula.prec, linear.prec)

blockwise_Frob <- function(x, M){
  # x: matrix
  # M: block size
  p <- dim(x)[1]
  if (p %% M != 0){
    stop("dimension of matrix cannot be divided by block size")
  }
  n.b <- p / M
  result <- matrix(0, nrow=n.b, ncol=n.b)
  for (i in 1:n.b){
    for (j in 1:n.b){
      result[i, j] <- norm(as.matrix(x[((i-1)*M+1):(i*M), ((j-1)*M+1):(j*M)]), "F")
    }
  }
  return(result)
}

colorpal <- brewer.pal(7,"Blues")
colorpal <- c("#FFFFFF",colorpal)
Methods <- c("cmc", 'copula', 'linear')

for (i in 1:length(Methods)){
  method <- Methods[i]
  prec <- prec.list[[i]]
  png(paste("~/work/CMC-GGM/Gene_Protein/htmap_", method,".png", sep = ""))
  heatmap(prec, Rowv = NA, Colv = NA, col= colorpal, 
          revC=TRUE, labRow = FALSE, labCol = FALSE, symm =T)
  dev.off()
}

cmc.sparsity <- est.cmc.glasso$sparsity
copula.sparsity <- est.copula.glasso$sparsity
linear.sparsity <- est.linear.glasso$sparsity

cat(paste("cmc-ggm sparsity:", cmc.sparsity, ";ggm sparsity:", linear.sparsity))
