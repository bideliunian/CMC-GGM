
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
library(dplyr)


function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Textures"
save.path <- "~/work/CMC-GGM/Textures/CVResult"
reading.path <- "~/work/CMC-GGM/Textures/Colored Brodatz"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)

arg <- as.numeric(commandArgs(trailingOnly=TRUE))

# 
# str_name <- paste(reading.path, '/D79_COLORED.tif', sep = "")
# grey <- raster(str_name)
# rgb <- brick(str_name)
# 
# plotRGB(rgb)
# sub_rgb <- crop(rgb, extent(rgb, 481, 640, 481, 640))
# plotRGB(sub_rgb)

load(paste(working.path, "/D79rgb.RData", sep = ""))
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
lambda_list <- seq(0.2, 0.65, 0.05)
lambda <- lambda_list[arg]
n.fold <- 5

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


cverror.cmc <- CV4tuning(cmc.score, p, lambda, n.fold)
cverror.copula <- CV4tuning(copula.score, p, lambda, n.fold)
cverror.linear <- CV4tuning(linear.score, p, lambda, n.fold)
cverror.list <- c(cverror.cmc, cverror.copula, cverror.linear)

save(cverror.list, file=paste(save.path,"/D79_cverror",arg,".Rdata",sep=""))

