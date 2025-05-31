library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)

mod.list <- c("a", "b", "c")
p.list <- c(100)
trans.list <- c("gauss", "copula", "cmcopula")


## aggregated dataframe the beginning and ending
df2 = expand.grid(trans.methods = c("cmc", "copula",'linear'), 
                  methods = c("thres","bglasso",'and','or'), TP = c(0,1), FP = c(0,1))
df2 = df2[df2$TP == df2$FP,]
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
    
    #### thresholding and glasso
    working.path <- paste("~/work/CMC-GGM/Model_",mod.alp,sep="")
    save.path <- paste("~/work/CMC-GGM/Model_",mod.alp,"/Results/p",p, sep="")
    
    n.par <- length(array)
    
    if(p==50) pee <- "050"
    if(p==100) pee <- "100"
    
    x.list <- c(seq(0, 0.2, by=0.01), seq(0.22, 1, by=0.02))
    y.cmc.g <-  y.cmc.t <- y.copula.g <- y.copula.t <- y.na.g <- y.na.t <- 0
    for(i in 1:n.par){
      load(paste(save.path,"/ROC.",mod.alp,".", pee, ".", trans.alp,".seed",i,".Rdata",sep=""))
      
      y.cmc.g <- y.cmc.g + interp(df[df$trans.methods=='cmc'&df$methods=='bglasso',3], df[df$trans.methods=='cmc'&df$methods=='bglasso',4], x.list) / n.par
      y.copula.g <- y.copula.g + interp(df[df$trans.methods=='copula'&df$methods=='bglasso',3], df[df$trans.methods=='copula'&df$methods=='bglasso',4], x.list) / n.par
      y.na.g <- y.na.g + interp(df[df$trans.methods=='linear'&df$methods=='bglasso',3], df[df$trans.methods=='linear'&df$methods=='bglasso',4], x.list) / n.par
      y.cmc.t <- y.cmc.t + interp(df[df$trans.methods=='cmc'&df$methods=='thres',3], df[df$trans.methods=='cmc'&df$methods=='thres',4], x.list) / n.par
      y.copula.t <- y.copula.t + interp(df[df$trans.methods=='copula'&df$methods=='thres',3], df[df$trans.methods=='copula'&df$methods=='thres',4], x.list) / n.par
      y.na.t <- y.na.t + interp(df[df$trans.methods=='linear'&df$methods=='thres',3], df[df$trans.methods=='linear'&df$methods=='thres',4], x.list) / n.par
    }
    
    #### neighborhood lasso
    save.path.2 <- paste("~/work/CMC-GGM/Model_",mod.alp,"_nblasso","/Results", sep="")
    load(paste(save.path.2,"/ROC.",mod.alp,".", pee, ".", trans.alp,".nbd.Rdata",sep=""))
    
    x.list <- c(seq(0, 0.2, by=0.01), seq(0.22, 1, by=0.02))
    y.cmc.and <-  y.cmc.or <- y.copula.and <- y.copula.or <- y.na.and <- y.na.or <- 0
    for(i in 1:n.par){
      df <- df.list[[i]]
      y.cmc.and <- y.cmc.and + interp(df[df$trans.methods=='cmc'&df$methods=='and',3], df[df$trans.methods=='cmc'&df$methods=='and',4], x.list) / n.par
      y.copula.and <- y.copula.and + interp(df[df$trans.methods=='copula'&df$methods=='and',3], df[df$trans.methods=='copula'&df$methods=='and',4], x.list) / n.par
      y.na.and <- y.na.and + interp(df[df$trans.methods=='linear'&df$methods=='and',3], df[df$trans.methods=='linear'&df$methods=='and',4], x.list) / n.par
      y.cmc.or <- y.cmc.or + interp(df[df$trans.methods=='cmc'&df$methods=='or',3], df[df$trans.methods=='cmc'&df$methods=='or',4], x.list) / n.par
      y.copula.or <- y.copula.or + interp(df[df$trans.methods=='copula'&df$methods=='or',3], df[df$trans.methods=='copula'&df$methods=='or',4], x.list) / n.par
      y.na.or <- y.na.or + interp(df[df$trans.methods=='linear'&df$methods=='or',3], df[df$trans.methods=='linear'&df$methods=='or',4], x.list) / n.par
    }
    
    
    df.comb = data.frame(trans.methods = rep(rep(c("cmc","copula",'linear'), 4), each=length(x.list)),
                       methods = rep(c("thres","bglasso","and","or"),each=3*length(x.list)),
                       TP = c(y.cmc.t, y.copula.t, y.na.t, y.cmc.g, y.copula.g, y.na.g,
                              y.cmc.and, y.copula.and, y.na.and, y.cmc.or, y.copula.or, y.na.or),
                       FP = rep(x.list,12))
    df.comb <- rbind(df.comb,df2)
    
    ####### choosing plotting methods #############
    df.comb.select <- df.comb[df.comb$methods %in% c("bglasso", "and"),]
    
    plot <- ggplot(df.comb.select, aes(x=FP, y=TP)) +
    geom_line(aes(color=trans.methods, linetype=methods))+ xlab("FPR") + ylab("TPR") + 
      theme_bw() +       scale_linetype_manual(values=c(1,2,3)) +
      scale_linetype_discrete(labels=c('group-glasso', 'nbd-group-lasso')) 
    
    plot.list[[index]] = plot
    index = index + 1
  }
}

  
pdf("~/work/CMC-GGM/Plot/roc_p100_select.pdf")
do.call("ggarrange", c(plot.list, ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom"))   
dev.off()



pdf("~/work/CMC-GGM/Plot/roc_p100_select_gauss.pdf")
do.call("ggarrange", c(list(plot.list[[1]], plot.list[[4]], plot.list[[7]]), ncol = 3, nrow = 3, common.legend = TRUE, legend="top"))   
dev.off()


pdf("~/work/CMC-GGM/Plot/roc_p100_select_copula.pdf")
do.call("ggarrange", c(list(plot.list[[2]], plot.list[[5]], plot.list[[8]]), ncol = 3, nrow = 3, common.legend = TRUE, legend="top"))   
dev.off()

pdf("~/work/CMC-GGM/Plot/roc_p100_select_cmc.pdf")
do.call("ggarrange", c(list(plot.list[[3]], plot.list[[6]], plot.list[[9]]), ncol = 3, nrow = 3, common.legend = TRUE, legend="top"))   
dev.off()
