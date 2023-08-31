library(xtable)
function.path <- "~/work/CMC-GGM/Functions"
source(paste(function.path,"auc.R", sep="/"))

trans.list <- c('exp', 'power')
d.list <- c(3, 5, 7)

grid <- expand.grid(d.list, trans.list)
names(grid) <- c('d', 'transform')
auc.summary <- numeric()

for(trans.alp in trans.list){
  for (d in d.list) {
      
    working.path <- "~/work/CMC-GGM/Model_a_proj"
    save.path <- "~/work/CMC-GGM/Model_a_proj/Results"
    
    files.name <- list.files(path=save.path)
    n.par <- 50
  
    if(d==7) dee <- "07"
    if(d==5) dee <- "05"
    if(d==3) dee <- "03"
    if(d==4) dee <- "04"
    if(d==10) dee <- "10"
    
    auc.mat <- numeric()
    for(i in 1:n.par){
      load(paste(save.path,"/ROC_a_d",dee,"_", trans.alp, ".seed",i,".Rdata",sep=""))
      
      roc <- sapply(split(df, list(df$trans.methods, df$methods)), FUN = auc.calculator)
      auc <- unlist(roc[1,])
      auc.mat <- rbind(auc.mat, auc)
    }
    
    save(auc.mat, file = paste(save.path,"/AUC", trans.alp, d,".Rdata",sep=""))
    auc.mean <- colMeans(auc.mat)
    auc.sd <- apply(auc.mat, 2, sd)
    auc.mean.sd <- c(auc.mean, auc.sd)
    auc.summary <- rbind(auc.summary, auc.mean.sd)
  }
}


auc.summary <- round(auc.summary,3)
grid.summary <- cbind(grid, auc.summary)
row.names(grid.summary) <- c()

summary.num <- grid.summary[3:ncol(grid.summary)]

result <- matrix(0, nrow = 2*nrow(summary.num), ncol = ncol(summary.num)/2)
for (i in 1:nrow(summary.num)) {
  for (j in 1:ncol(result)) {
    result[2*i - 1, j] <- summary.num[i,j]
    result[2*i, j] <- paste("(", summary.num[i, ncol(result) + j],")", sep = "")
  }
}
result <- cbind(rep(grid.summary[,1], each=2), rep(grid.summary[,2], each=2), result)
print(xtable(result))
