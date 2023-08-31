library(ggplot2)
library(gridExtra)
library(dplyr)

mod.list <- c("a", "b", "c")
p.list <- c(100)
trans.list <- c("gauss", "copula", "cmcopula")


## aggregated dataframe the beginning and ending
df2 = expand.grid(trans.methods = c("cmc", "copula",'linear'), 
                  methods = c("thres","bglasso"), TP = c(0,1), FP = c(0,1))
df2 = df2[df2$TP == df2$FP,]

######################
## plot the roc curve with median best auc score
##########################

## function to get median index
# get.median.index <- function(x){
#   n <- length(x)
#   idx <- floor(n/2)
#   return(which(rank(x) == idx))
# }
# 
# p = 100
# plot.list = list()
# index = 1
# for(mod.alp in mod.list){
#   for(trans.alp in trans.list){
#     array <- c(1:50)
#     
#     working.path <- paste("~/work/CMC-GGM/Model_",mod.alp,sep="")
#     save.path <- paste("~/work/CMC-GGM/Model_",mod.alp,"/Results/p",p, sep="")
#     
#     n.par <- length(array)
# 
#     if(p==50) pee <- "050"
#     if(p==100) pee <- "100"
#     
#     load(paste(save.path,"/AUC", trans.alp,".Rdata",sep=""))
#     median.index <- apply(auc.mat, 2, get.median.index)
#     run.ind <- c(median.index[4:6],median.index[1:3])
#     df.comb <- df2
#     
#     for (j in 1:6) {
#       load(paste(save.path,"/ROC.",mod.alp,".",p,".", trans.alp,".seed",run.ind[j],".Rdata", sep=""))
#       df.temp <- df[(20*(j-1)+1):(20*j),]
#       df.comb <- rbind(df.comb, df.temp)
#     }
#     
#     df.comb <- df.comb[order(df.comb$trans.methods,df.comb$methods, df.comb$FP, df.comb$TP),]
#       
#     plot <- ggplot(df.comb, aes(x=FP, y=TP)) +
#     geom_line(aes(color=trans.methods,linetype=methods))+ scale_linetype_manual(values=c("dashed","solid"))+
#       xlab("FPR") + ylab("TPR")
#     
#     plot.list[[index]] = plot
#     index = index + 1
#   }
# }


#################################
## plot roc curve by average
#################################

# function: interpolation
interp <- function(TPR, FPR, x.list){
  # input: lists; output: y.list
  TPR <- sort(TPR)
  FPR <- sort(FPR)
  TPR <- c(0,TPR,1)
  FPR <- c(0,FPR,1)
  
  
  len.TPR <- length(TPR)
  len.x <- length(x.list)
  y.list <- rep(1, len.x)
  for(i in 1:(len.x-1)){
    for(j in 1:(len.TPR-1)){
      if((x.list[i] >= FPR[j]) & (x.list[i] < FPR[j+1])){
        int.left <- j
        int.right <- j+1
        break
      }
    }
    y.list[i] <- TPR[int.left] + (x.list[i]-FPR[int.left])*(TPR[int.right] - TPR[int.left])/(FPR[int.right] - FPR[int.left])
  }
  return(y.list)
}

p = 100
plot.list = list()
index = 1
for(mod.alp in mod.list){
  for(trans.alp in trans.list){
    array <- c(1:50)
    
    working.path <- paste("~/work/CMC-GGM/Model_",mod.alp,sep="")
    save.path <- paste("~/work/CMC-GGM/Model_",mod.alp,"/Results/p",p, sep="")
    
    n.par <- length(array)
    
    if(p==50) pee <- "050"
    if(p==100) pee <- "100"
    
    x.list <- c(seq(0, 0.2, by=0.01), seq(0.22, 1, by=0.02))
    y.cmc.g <- 0; y.cmc.t <- 0; y.copula.g <- 0; y.copula.t <- 0; y.na.g <- 0; y.na.t <- 0
    for(i in 1:n.par){
      load(paste(save.path,"/ROC.",mod.alp,".", pee, ".", trans.alp,".seed",i,".Rdata",sep=""))
      
      y.cmc.g <- y.cmc.g + interp(df[df$trans.methods=='cmc'&df$methods=='bglasso',3], df[df$trans.methods=='cmc'&df$methods=='bglasso',4], x.list) / n.par
      y.copula.g <- y.copula.g + interp(df[df$trans.methods=='copula'&df$methods=='bglasso',3], df[df$trans.methods=='copula'&df$methods=='bglasso',4], x.list) / n.par
      y.na.g <- y.na.g + interp(df[df$trans.methods=='linear'&df$methods=='bglasso',3], df[df$trans.methods=='linear'&df$methods=='bglasso',4], x.list) / n.par
      y.cmc.t <- y.cmc.t + interp(df[df$trans.methods=='cmc'&df$methods=='thres',3], df[df$trans.methods=='cmc'&df$methods=='thres',4], x.list) / n.par
      y.copula.t <- y.copula.t + interp(df[df$trans.methods=='copula'&df$methods=='thres',3], df[df$trans.methods=='copula'&df$methods=='thres',4], x.list) / n.par
      y.na.t <- y.na.t + interp(df[df$trans.methods=='linear'&df$methods=='thres',3], df[df$trans.methods=='linear'&df$methods=='thres',4], x.list) / n.par
    }
   
    df.comb = data.frame(trans.methods = rep(c("cmc","copula",'linear',"cmc","copula",'linear'),each=length(x.list)),
                       methods = rep(c("thres","bglasso"),each=3*length(x.list)),
                       TP = c(y.cmc.t, y.copula.t, y.na.t, y.cmc.g, y.copula.g, y.na.g),
                       FP = rep(x.list,6))
    df.comb <- rbind(df.comb,df2)
    plot <- ggplot(df.comb, aes(x=FP, y=TP)) +
    geom_line(aes(color=trans.methods,linetype=methods))+ xlab("FPR") + ylab("TPR") + theme_bw()#scale_linetype_manual(values=c("dashed","solid"))+ 
    
    plot.list[[index]] = plot
    index = index + 1
  }
}


  
pdf("~/work/CMC-GGM/Plot/roc_p100.pdf")
do.call("ggarrange", c(plot.list, ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom"))   
dev.off()

