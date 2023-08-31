library(xtable)
function.path <- "~/work/CMC-GGM/Functions"
source(paste(function.path,"auc.R", sep="/"))

mod.list <- c("a", "b", "c")
p.list <- c(50, 100)
trans.list <- c("gauss", "copula", "cmcopula")

grid <- expand.grid(p.list, trans.list, mod.list)
names(grid) <- c('p', 'transform', 'model')
auc.summary <- numeric()

for(mod.alp in mod.list){
  for(trans.alp in trans.list){
    for (p in p.list) {
      
      working.path.2 <- paste("~/work/CMC-GGM/Model_",mod.alp,sep="")
      save.path <- paste("~/work/CMC-GGM/Model_",mod.alp,"/Results/p",p, sep="")
      
      n.par <- 50
    
      if(p==50) pee <- "050"
      if(p==100) pee <- "100"
      
      auc.mat <- numeric()
      for(i in 1:n.par){
        load(paste(save.path,"/ROC.",mod.alp,".", pee, ".", trans.alp,".seed",i,".Rdata",sep=""))
        
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
grid.summary <- cbind(grid, auc.summary)
row.names(grid.summary) <- c()

table.names <- c("p","transform","model","cmcopula.g","cmcopula.g.sd","copula.g","copula.g.sd",
                 "linear.g","linear.g'sd","cmcopula.t","cmcopula.t.sd","copula.t","copula.t.sd",
                 "linear.t","linear.t.sd")

#################### Model A ########################
## write out result
model.a.result.pre <- grid.summary[grid.summary$model == 'a',]
model.a.result <- model.a.result.pre
for(i in 1:6){
  model.a.result[,3+2*i-1] <- model.a.result.pre[,3+i]
  model.a.result[,3+2*i] <- model.a.result.pre[,9+i]
}
names(model.a.result) <- table.names

model.a.num <- model.a.result[4:15]
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
model.b.result.pre <- grid.summary[grid.summary$model == 'b',]
model.b.result <- model.b.result.pre
for(i in 1:6){
  model.b.result[,3+2*i-1] <- model.b.result.pre[,3+i]
  model.b.result[,3+2*i] <- model.b.result.pre[,9+i]
}
names(model.b.result) <- table.names

model.b.num <- model.b.result[4:15]
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
model.c.result.pre <- grid.summary[grid.summary$model == 'c',]
model.c.result <- model.c.result.pre
for(i in 1:6){
  model.c.result[,3+2*i-1] <- model.c.result.pre[,3+i]
  model.c.result[,3+2*i] <- model.c.result.pre[,9+i]
}

names(model.c.result) <- table.names

model.c.num <- model.c.result[4:15]
result_c <- matrix(0, nrow = 2*nrow(model.c.num), ncol = ncol(model.c.num)/2)
for (i in 1:nrow(model.c.num)) {
  for (j in 1:ncol(result_c)) {
    result_c[2*i-1,j] <- model.c.num[i,2*j-1]
    result_c[2*i, j] <- paste("(", model.c.num[i,2*j],")", sep = "")
  }
}
result_c <- cbind(rep(model.c.result[,1], each = 2), result_c)
print(xtable(result_c))

load("~/work/CMC-GGM/Plot/result_nb_and.Rdata")
result_a_all <- cbind(result_a, result_a_nb_and)
result_b_all <- cbind(result_b, result_b_nb_and)
result_c_all <- cbind(result_c, result_c_nb_and)

cat("Model A auc table:")
print(xtable(result_a_all))

cat("Model B auc table")
print(xtable(result_b_all))

cat("Model C auc table:")
print(xtable(result_c_all))


################# compare with rank cmc ########################
load("~/work/CMC-GGM/Plot/result_rank.Rdata")

result_a_rankcmc <- cbind(result_a_all[result_a_all[,1]=="100",c(1,2,5,8)], result_a_rank)
result_b_rankcmc <- cbind(result_b_all[result_a_all[,1]=="100",c(1,2,5,8)], result_b_rank)
result_c_rankcmc <- cbind(result_c_all[result_a_all[,1]=="100",c(1,2,5,8)], result_c_rank)

cat("Model A cmc and rankcmc auc table:")
print(xtable(result_a_rankcmc))

cat("Model B cmc and rankcmc auc table:")
print(xtable(result_b_rankcmc))

cat("Model C cmc and rankcmc auc table:")
print(xtable(result_c_rankcmc))
