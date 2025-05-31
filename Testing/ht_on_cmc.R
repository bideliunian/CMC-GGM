#######################
## Model_a ROC (Model 1)
## p = 50, n = 100
## ggm
###########################
## name of paths
function.path <- "~/Projects/CMC-GGM/Functions"
working.path <- "~/Projects/CMC-GGM/Testing"
save.path <- "~/Projects/CMC-GGM/Testing"

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
library(energy)

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)
source("~/Projects/CMC-GGM/Testing/funcs.r")
source("~/Projects/CMC-GGM/Testing/funcs2.r")

## global parameters
n = 200
p = 50
M = 2
rand.seed = 1
model = 'model1'
trans.type = 'cmcopula'
trans.func = 'exp'

data = gen_data(n=n, p=p, m=M, model=model, run.ind = rand.seed + 1, 
                trans.type = trans.type, trans.func = trans.func)
data.X = data$X

p.values.cmc <- p.values.copula <- p.values.linear <- numeric()

for (i in 1:100) {
  data = gen_data(n=n, p=p, m=M, model=model, run.ind = rand.seed + i, trans.type = trans.type)
  data.X = data$X
  G = data$Omega0
  Omega = data$Omega
  
  data.cmcopula = mult.trans(X=data.X, p = p, fun = 'cmcopula')
  data.copula = mult.trans(X=data.X, p = p, fun = 'copula')
  data.linear = mult.trans(X=data.X, p = p, fun = 'scale')
  
  p.values.cmc <- c(p.values.cmc, mvnorm.etest(data.cmcopula, R=300)$p.value)
  #p.values.cmc.new <- c(p.values.cmc.new, getp(data.cmcopula, L=1, BB=500)$Yp[1])
  p.values.copula <- c(p.values.copula, mvnorm.etest(data.copula, R=300)$p.value)
  #p.values.copula.new <- c(p.values.copula.new, getp(data.copula, L=1, BB=500)$Yp[1])
  p.values.linear <- c(p.values.linear, mvnorm.etest(data.linear, R=300)$p.value)
  #p.values.linear.new <- c(p.values.linear.new, getp(data.linear, L=1, BB=500)$Yp[1])
  cat("number of experiments:", i)
}

cat("cmc transformation:", "\n", mean(p.values.cmc>0.05), "\n")
cat("copula transformation:", "\n", mean(p.values.copula>0.05), "\n")
cat("linear transformation:", "\n", mean(p.values.linear>0.05), "\n")


# df.pvalues.linear <- data.frame(pvalues = p.values.linear)
# ## histogram of the p.values after copula transformation
# 
# p.linear <- ggplot(df.pvalues.linear, aes(x=pvalues)) + geom_histogram(breaks=seq(0,0.8,by=0.025)) +
#   geom_vline(aes(xintercept=0.05, color='red'), linetype='dashed', size=1) +
#   theme(legend.position="none", text = element_text(size = 20))
