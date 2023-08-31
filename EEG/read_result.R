############### PART 0: library, source and pathes ########################

library(dplyr)
library(pdist)
library(pracma)
library(clue)
library(BB)
library(mvtnorm)
library(doParallel)
library(randtoolbox)

chanel.path <- "~/work/CMC-GGM/EEG/eeg_processed/"
reading.path <- "~/work/CMC-GGM/EEG/Result/"
working.path <- "~/work/CMC-GGM/EEG/"
function.path <- "~/work/CMC-GGM/Functions/"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)


############### PART 1: Reading Data ####################
n.control <- 45
n.alcohol <- 77
p <- 64
tau <- 256
M <- 6
n.expr <- 50
alc.pcmc.prec.freq <- alc.copula.prec.freq <- alc.cmc.prec.freq <- matrix(0, p, p)
ctr.pcmc.prec.freq <- ctr.copula.prec.freq <- ctr.cmc.prec.freq <- matrix(0, p, p)

for (i in 1:n.expr) {
  load(paste(reading.path,"/prec.bt.seed",i,".Rdata",sep=""))
  
  alc.pcmc.prec <- prec.list[[1]]
  alc.cmc.prec <- prec.list[[2]]
  alc.copula.prec <- prec.list[[3]]
  
  ctr.pcmc.prec <- prec.list[[4]]
  ctr.cmc.prec <- prec.list[[5]]
  ctr.copula.prec <- prec.list[[6]]
  
  alc.pcmc.prec.freq <- alc.pcmc.prec.freq + alc.pcmc.prec / n.expr
  alc.cmc.prec.freq <- alc.cmc.prec.freq + alc.cmc.prec / n.expr
  alc.copula.prec.freq <- alc.copula.prec.freq + alc.copula.prec / n.expr
  
  ctr.pcmc.prec.freq <- ctr.pcmc.prec.freq + ctr.pcmc.prec / n.expr
  ctr.cmc.prec.freq <- ctr.cmc.prec.freq + ctr.cmc.prec / n.expr
  ctr.copula.prec.freq <- ctr.copula.prec.freq + ctr.copula.prec / n.expr
}


############# PART 2: Summary statistics #########################
thres.alc <- 0.5
thres.ctr <- 0.5

robust.prec <- function(prec.freq, threshold){
  prec.robust <- prec.freq
  prec.robust[prec.robust > threshold] <- 1
  prec.robust[prec.robust < threshold] <- 0
  
  return(prec.robust)
}

alc.pcmc.prec.robust <- robust.prec(alc.pcmc.prec.freq, 0.9)
alc.cmc.prec.robust <- robust.prec(alc.cmc.prec.freq, 0.2)
alc.copula.prec.robust <- robust.prec(alc.copula.prec.freq, 0.99)

ctr.pcmc.prec.robust <- robust.prec(ctr.pcmc.prec.freq, 0.9)
ctr.cmc.prec.robust <- robust.prec(ctr.cmc.prec.freq, 0.2)
ctr.copula.prec.robust <- robust.prec(ctr.copula.prec.freq, 0.99)

####### number of edges ####################
alc.num.edge.pcmc <- sum(alc.pcmc.prec.robust)/2
alc.num.edge.cmc <- sum(alc.cmc.prec.robust)/2
alc.num.edge.copula <- sum(alc.copula.prec.robust)/2

ctr.num.edge.pcmc <- sum(ctr.pcmc.prec.robust)/2
ctr.num.edge.cmc <- sum(ctr.cmc.prec.robust)/2
ctr.num.edge.copula <- sum(ctr.copula.prec.robust)/2

print(paste("pcmc alcohol sparsity:", 2*alc.num.edge.pcmc/(p^2-p)))
print(paste("pcmc control sparsity:", 2*ctr.num.edge.pcmc/(p^2-p)))
print(paste("cmc alcohol sparsity:", 2*alc.num.edge.cmc/(p^2-p)))
print(paste("cmc control sparsity:", 2*ctr.num.edge.cmc/(p^2-p)))
print(paste("copula alcohol sparsity:", 2*alc.num.edge.copula/(p^2-p)))
print(paste("copula control sparsity:", 2*ctr.num.edge.copula/(p^2-p)))

############# number of shared edges ##################
edge.share.pcmc <- alc.pcmc.prec.robust & ctr.pcmc.prec.robust
edge.share.cmc <- alc.cmc.prec.robust & ctr.cmc.prec.robust
edge.share.copula <- alc.copula.prec.robust & ctr.copula.prec.robust

cat(paste("edge numbers shared by pcmc", sum(edge.share.pcmc)/2))
cat(paste("edge numbers shared by cmc", sum(edge.share.cmc)/2))
cat(paste("edge numbers shared by copula", sum(edge.share.copula)/2))


############ PART 3: plotting ###################################
library(reshape2)
library(ggplot2)

load(file=paste(chanel.path,"nodes_name.Rdata",sep=""))

# "#F8766D" "#00BA38" "#619CFF"
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
network.plotting <- function(AdjMat1, AdjMat2, node.names, layMat){
  ## get the intersection of two adjacant matrices/ two edge adjmat
  p <- nrow(AdjMat1)
  AdjMat.and <- matrix(as.numeric(AdjMat1 & AdjMat2), nrow=p, ncol=p)
  AdjMat1.edge <- AdjMat1 - AdjMat.and
  AdjMat2.edge <- AdjMat2 - AdjMat.and
  
  row.names(AdjMat.and) <- node.names
  colnames(AdjMat.and) <- node.names
  row.names(AdjMat1.edge) <- node.names
  colnames(AdjMat1.edge) <- node.names
  row.names(AdjMat2.edge) <- node.names
  colnames(AdjMat2.edge) <- node.names
  
  net.and <- graph_from_adjacency_matrix(AdjMat.and, mode="undirected")
  net1 <- graph_from_adjacency_matrix(AdjMat1.edge, mode="undirected")
  net2 <- graph_from_adjacency_matrix(AdjMat2.edge, mode="undirected")
  
  #V(net.and)$label.cex <- 1
  plot(net.and, edge.color="#619CFF", edge.width=2, vertex.size=3, margin=0, layout=layMat,  vertex.color='grey90', 
       edge.curved=0, vertex.label.dist=1, vertex.label.cex=1)
  plot(net1, edge.color="#F8766D", edge.width=2, vertex.size=3, margin=0, layout=layMat, vertex.color='grey90',
       edge.curved=0.3, vertex.label.dist=1, vertex.label.cex=1, add=T)
  plot(net2, edge.color="#00BA38", edge.width=2, vertex.size=3, margin=0, layout=layMat, vertex.color='grey90',
       edge.curved=0.3, vertex.label.dist=1, vertex.label.cex=1, add=T)
}
## ploting the connected graph from control group and alcohol group together

pdf(paste("~/work/CMC-GGM/EEG/graph_pcmc.pdf"))
network.plotting(alc.pcmc.prec.robust, ctr.pcmc.prec.robust, node.names = node.names, layMat = layMat)
dev.off()

pdf(paste("~/work/CMC-GGM/EEG/graph_cmc.pdf"))
network.plotting(alc.cmc.prec.robust, ctr.cmc.prec.robust, node.names = node.names, layMat = layMat)
dev.off()

pdf(paste("~/work/CMC-GGM/EEG/graph_copula.pdf"))
network.plotting(alc.copula.prec.robust, ctr.copula.prec.robust, node.names = node.names, layMat = layMat)
dev.off()
