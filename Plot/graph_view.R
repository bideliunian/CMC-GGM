#######################
## plot the graph structure
###########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Plot"
save.path <- "~/work/CMC-GGM/Plot"

## packages
library(pracma)
library(mvtnorm)
library(pdist)
library(RColorBrewer)

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

## global parameters
n = 100
p = 50
M = 2
rand.seed = 2022
models = c('model1', 'model2', 'model3')
trans.type = 'no_trans'

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

colorpal<-brewer.pal(7,"Blues")
colorpal<-c("#FFFFFF",colorpal)

for (i in 1:length(models)){
  model = models[i]
  data = gen_data(n=n, p=p, m=M, model=model, run.ind = rand.seed, trans.type = 'no_trans')
  Omega = data$Omega
  Omega_F_norm = blockwise_Frob(Omega, M=M)
  diag(Omega_F_norm) = 0
  png(paste("~/work/CMC-GGM/Plot/htmap_", model,".png"))
  heatmap(Omega_F_norm, Rowv = NA, Colv = NA, col= colorpal, 
          revC=TRUE, labRow = FALSE, labCol = FALSE, symm =T)
  dev.off()
}
