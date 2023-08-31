library(xtable)
function.path <- "~/work/CMC-GGM/Functions"
source(paste(function.path,"auc.R", sep="/"))

mod.list <- c("a", "b", "c")
p.list <- c(100)
trans.list <- c("gauss", "copula", "cmcopula")

grid <- expand.grid(p.list, trans.list, mod.list)
names(grid) <- c('p', 'transform', 'model')
auc.summary <- numeric()

for(mod.alp in mod.list){
  for(trans.alp in trans.list){
    for (p in p.list) {
      
      working.path.2 <- paste("~/work/CMC-GGM/Model_",mod.alp,"_rank", sep="")
      save.path <- paste("~/work/CMC-GGM/Model_",mod.alp,"_rank", "/Results/p",p, sep="")
      
      n.par <- 50
      pee <- "100"
      
      auc.mat <- numeric()
      for(i in 1:n.par){
        load(paste(save.path,"/ROC.",mod.alp,".", pee, ".", trans.alp,".seed",i,".Rdata",sep=""))
        
        df = data.frame(trans.methods = rep('rankcmc',60),
                        methods = rep(c("thresholding","group-glasso",'nbd-group-lasso'),each=20), 
                        TP = df$TP[1:60], FP = df$FP[1:60])
        
        roc <- sapply(split(df, list(df$trans.methods, df$methods)), FUN = auc.calculator)
        auc <- unlist(roc[1,])
        auc.mat <- rbind(auc.mat, auc)
      }
      
      save(auc.mat, file = paste(save.path,"/AUC", trans.alp,".Rdata",sep=""))
      auc.mean <- colMeans(auc.mat)
      auc.sd <- apply(auc.mat, 2, sd)
      auc.mean.sd <- c(auc.mean, auc.sd)
      auc.summary <- rbind(auc.summary, auc.mean.sd) 
    }
  }
}

auc.summary <- round(auc.summary,3)
grid.summary.rank <- cbind(grid, auc.summary)
names <- c("p", "transform", "model","rankcmc.group-glasso", "rankcmc.nbd-group-lasso",
                              "rankcmc.thresholding", "rankcmc.group-glasso_sd", "rankcmc.nbd-group-lasso_sd",
                              "rankcmc.thresholding_sd") 
names(grid.summary.rank) <- names
shuffle_names <- c("p", "transform", "model", 
                   "rankcmc.group-glasso",  "rankcmc.group-glasso_sd", 
                  "rankcmc.thresholding", "rankcmc.thresholding_sd",
                  "rankcmc.nbd-group-lasso", "rankcmc.nbd-group-lasso_sd") 
grid.summary.rank <- grid.summary.rank[, shuffle_names]  
row.names(grid.summary.rank) <- c()

#################### Model A ########################
## write out result
model.a.result <- grid.summary.rank[grid.summary.rank$model == 'a',]
model.a.num <- model.a.result[4:9]
result_a <- matrix(0, nrow = 2*nrow(model.a.num), ncol = ncol(model.a.num)/2)
for (i in 1:nrow(model.a.num)) {
  for (j in 1:ncol(result_a)) {
    result_a[2*i-1,j] <- model.a.num[i,2*j-1]
    result_a[2*i, j] <- paste("(", model.a.num[i,2*j],")", sep = "")
  }
}
result_a <- cbind(rep(model.a.result[,1], each = 2), result_a)
print(xtable(result_a))


#################### Model B #####################
## write out model b result
model.b.result <- grid.summary.rank[grid.summary.rank$model == 'b',]
model.b.num <- model.b.result[4:9]
result_b <- matrix(0, nrow = 2*nrow(model.b.num), ncol = ncol(model.b.num)/2)
for (i in 1:nrow(model.b.num)) {
  for (j in 1:ncol(result_b)) {
    result_b[2*i-1,j] <- model.b.num[i,2*j-1]
    result_b[2*i, j] <- paste("(", model.b.num[i,2*j],")", sep = "")
  }
}
result_b <- cbind(rep(model.b.result[,1], each = 2), result_b)
print(xtable(result_b))

####################### Model C #########################
## write out model c result
model.c.result <- grid.summary.rank[grid.summary.rank$model == 'c',]
model.c.num <- model.c.result[4:9]
result_c <- matrix(0, nrow = 2*nrow(model.c.num), ncol = ncol(model.c.num)/2)
for (i in 1:nrow(model.c.num)) {
  for (j in 1:ncol(result_c)) {
    result_c[2*i-1,j] <- model.c.num[i,2*j-1]
    result_c[2*i, j] <- paste("(", model.c.num[i,2*j],")", sep = "")
  }
}
result_c <- cbind(rep(model.c.result[,1], each = 2), result_c)
print(xtable(result_c))

result_a_rank <- result_a[,2:4]
result_b_rank <- result_b[,2:4]
result_c_rank <- result_c[,2:4]
save(result_a_rank, result_b_rank, result_c_rank, 
     file =  "~/work/CMC-GGM/Plot/result_rank.Rdata")
