################################
## visualize the subspace estimation error
##################################

function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Model_a_proj"
read.path <- "~/work/CMC-GGM/proj_cmc_k2"
save.path <- "~/work/CMC-GGM/Model_a_proj/Results"

trans.list <- c('power','exp')
d.list <- c(3, 4, 5, 7, 10)
k <- 2
p <- 50

error.mat <- array(0, dim = c(length(trans.list), length(d.list), 5))

for (i in 1:length(trans.list)){
  trans <- trans.list[i]
  for (j in 1:length(d.list)){
    d <- d.list[j]
    
    if(d==10) {dee <- "10"}
    if(d==7) {dee <- "07"}
    if(d==6) {dee <- "06"}
    if(d==5) {dee <- "05"}
    if(d==4) {dee <- "04"}
    if(d==3) {dee <- "03"}
    
    ## true projection matrix
    proj.true <- matrix(0, nrow = d, ncol = d)
    proj.true[1, 1] =  proj.true[2, 2] = 1
    
    error <- rep(0, p)
    
    U.df <- read.csv(file=paste(read.path,"/proj_subspace_model1_", trans, "_d", dee, ".csv",sep=""), 
                    header = FALSE)
    #U.list = list()
    for (l in 1:p) {
      #U.list[[k]] = matrix(as.numeric(U.df[k,]), nrow = d, ncol = k)
      U = matrix(as.numeric(U.df[l,]), nrow = d, ncol = k)
      proj = U %*% t(U)
      error[l] = norm(proj - proj.true, 'F')
    }
    
    error.mat[i, j, ] = quantile(error)
  }
}

df.error.power = as.data.frame(error.mat[1, , ])
df.error.power$d = d.list
df.error.exp = as.data.frame(error.mat[2, , ])
df.error.exp$d = d.list

require(ggplot2)

png("~/work/CMC-GGM/Plot/proj_error_power.png")
ggplot(df.error.power, aes(x = d, y = V3)) +
  geom_point(size = 4) +
  theme_bw() + 
  geom_errorbar(aes(ymax = V4, ymin = V2))
dev.off()

png("~/work/CMC-GGM/Plot/proj_error_exp.png")
ggplot(df.error.exp, aes(x = d, y = V3)) +
  geom_point(size = 4) +
  theme_bw() + 
  geom_errorbar(aes(ymax = V4, ymin = V2))
dev.off()
