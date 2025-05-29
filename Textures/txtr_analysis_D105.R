################ PART 0: Preparing ###################
## packages
library(pdist)
library(pracma)
library(clue)
library(huge)
library(igraph)
library(BB)
library(mvtnorm)
library(doParallel)
library(randtoolbox)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)

function.path <- "/Users/qiz/CMC-GGM/Functions"
working.path <- "/Users/qiz/CMC-GGM/Textures"
save.path <- "/Users/qiz/CMC-GGM/Textures"
reading.path <- "/Users/qiz/CMC-GGM/Textures/CVResult"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)
source(paste(working.path, "/graph_est.R", sep = ""))
source(paste(working.path, "/network_plotting.R", sep = ""))
# str_name <- paste(reading.path, '/D79_COLORED.tif', sep = "")
# grey <- raster(str_name)
# rgb <- brick(str_name)
# 
# plotRGB(rgb)
# sub_rgb <- crop(rgb, extent(rgb, 481, 640, 481, 640))
# plotRGB(sub_rgb)

################## PART 1: Processing Data ##################
load(paste(save.path, "/D105rgb.RData", sep = ""))

r_data <- rgb_array[,,1]
g_data <- rgb_array[,,2]
b_data <- rgb_array[,,3]

###########dimension of RGB pixels ###############
cat(paste("dimension of RGB pixels:", dim(rgb_array)[1], dim(rgb_array)[2], dim(rgb_array)[3]))

length_block <- 8
M <- 3
p <- length_block^2
n <- (dim(rgb_array)[1] / length_block)^2
cat(paste("p=", p,";n=",n, sep = ""))


split_matrix <- function(M, length_block){
  n_block <- nrow(M) / length_block
  matrix_new <- numeric()
  for (i in c(1:n_block)) {
    index.i <- c(((i-1)*length_block + 1):(i*length_block))
    for (j in 1:n_block) {
      index.j <- c(((j-1)*length_block + 1):(j*length_block))
      matrix_new <- rbind(matrix_new, as.vector(M[index.i, index.j]))
    }
  }
  return(matrix_new)
}

r_data_split <- split_matrix(r_data, length_block)
g_data_split <- split_matrix(g_data, length_block)
b_data_split <- split_matrix(b_data, length_block)

data.combined <- matrix(0, n, p*M)
for (i in 1:p) {
  data.combined[, 3*i-2] <- r_data_split[, i]
  data.combined[, 3*i-1] <- g_data_split[, i]
  data.combined[, 3*i] <- b_data_split[, i]
}

cmc.score = mult.trans(X=data.combined, p = p, fun = 'cmcopula')
copula.score = mult.trans(X=data.combined, p = p, fun = 'copula')
linear.score = mult.trans(X=data.combined, p = p, fun = 'scale')

S.cmc = cor(cmc.score)
S.copula = cor(copula.score)
S.linear = cor(linear.score)

################## PART 2: Graph Estimation ########################

est.cmc.glasso <- graph.est(S=S.cmc, n=n, p=p, sparsity=0.2)
est.copula.glasso <- graph.est(S=S.copula+10^-3*diag(M*p), n=n, p=p, sparsity=0.2)
est.linear.glasso <- graph.est(S=S.linear, n=n, p=p, sparsity=0.2)

Theta.cmc <- est.cmc.glasso$Theta
Theta.copula <- est.copula.glasso$Theta
Theta.linear <- est.linear.glasso$Theta

edge.cmc <- blockwise_Frob(Theta.cmc, M)
edge.copula <- blockwise_Frob(Theta.copula, M)
edge.linear <- blockwise_Frob(Theta.linear, M)


cmc.sparsity <- est.cmc.glasso$sparsity
copula.sparsity <- est.copula.glasso$sparsity
linear.sparsity <- est.linear.glasso$sparsity

cat(paste("cmc-ggm sparsity:", cmc.sparsity, ";ggm sparsity:", linear.sparsity))

###################### PART 3: plotting ###############################
## graph
cols <- brewer.pal(9, 'RdPu')

p.cmc <- network.plotting(edge.cmc, cols)
p.copula <- network.plotting(edge.copula, cols)
p.linear <- network.plotting(edge.linear, cols)

#load(paste(working.path, "/0311D105.RData", sep = ""))
png(paste("/Users/qiz/CMC-GGM/Textures/graphmap_D105_cmc.png"), width = 360, height = 300)
p.cmc
dev.off()

png(paste("/Users/qiz/CMC-GGM/Textures/graphmap_D105_linear.png"), width = 360, height = 300)
p.linear
dev.off()

