library(ggplot2)
library(gridExtra)
library(dplyr)

d.list <- c(3, 5, 7)
trans.list <- c("copula","cmc","pcmc")
model.list <- c('power', 'exp')
est.list <- c('group-glasso', 'thresholding', 'nbd-group-lasso')
## aggregated dataframe the beginning and ending
df2 = expand.grid(trans.methods = trans.list, est.methods = est.list,
                  TP = c(0,1), FP = c(0,1))
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

p = 50
plot.list = list()
index = 4

for(model.alp in model.list){
  for(d in d.list){
  
    if(d==7) dee <- "07"
    if(d==5) dee <- "05"
    if(d==3) dee <- "03"
    if(d==4) dee <- "04"
    if(d==10) dee <- "10"
    
    array <- c(1:50)
    
    working.path <- "~/work/CMC-GGM/Model_a_proj"
    save.path <- "~/work/CMC-GGM/Model_a_proj/Results"
    
    n.par <- length(array)
    
    x.list <- c(seq(0, 0.2, by=0.01), seq(0.22, 1, by=0.02))
    y.cmc.g <- 0; y.copula.g <- 0; y.pcmc.g <- 0;
    y.cmc.t <- 0; y.copula.t <- 0; y.pcmc.t <- 0;
    y.cmc.nb <- 0; y.copula.nb <- 0; y.pcmc.nb <- 0;
    for(i in 1:n.par){
      load(paste(save.path,"/ROC_a_d",dee,"_", model.alp, ".seed",i,".Rdata",sep=""))
      
      y.copula.g <- y.copula.g + interp(df[df$trans.methods=='copula'&df$methods=='group-glasso',3], 
                                        df[df$trans.methods=='copula'&df$methods=='group-glasso',4], x.list) / n.par
      y.cmc.g <- y.cmc.g + interp(df[df$trans.methods=='cmc'&df$methods=='group-glasso',3], 
                                  df[df$trans.methods=='cmc'&df$methods=='group-glasso',4], x.list) / n.par
      y.pcmc.g <- y.pcmc.g + interp(df[df$trans.methods=='pcmc'&df$methods=='group-glasso',3], 
                                    df[df$trans.methods=='pcmc'&df$methods=='group-glasso',4], x.list) / n.par
      y.copula.t <- y.copula.t + interp(df[df$trans.methods=='copula'&df$methods=='thresholding',3], 
                                        df[df$trans.methods=='copula'&df$methods=='thresholding',4], x.list) / n.par
      y.cmc.t <- y.cmc.t + interp(df[df$trans.methods=='cmc'&df$methods=='thresholding',3], 
                                  df[df$trans.methods=='cmc'&df$methods=='thresholding',4], x.list) / n.par
      y.pcmc.t <- y.pcmc.t + interp(df[df$trans.methods=='pcmc'&df$methods=='thresholding',3], 
                                    df[df$trans.methods=='pcmc'&df$methods=='thresholding',4], x.list) / n.par
      y.cmc.nb <- y.cmc.nb + interp(df[df$trans.methods=='cmc'&df$methods=='nbd-group-lasso',3], 
                                     df[df$trans.methods=='cmc'&df$methods=='nbd-group-lasso',4], x.list) / n.par
      y.copula.nb <- y.copula.nb + interp(df[df$trans.methods=='copula'&df$methods=='nbd-group-lasso',3], 
                                           df[df$trans.methods=='copula'&df$methods=='nbd-group-lasso',4], x.list) / n.par
      y.pcmc.nb <- y.pcmc.nb + interp(df[df$trans.methods=='pcmc'&df$methods=='nbd-group-lasso',3], 
                                   df[df$trans.methods=='pcmc'&df$methods=='nbd-group-lasso',4], x.list) / n.par
    }
   
    df.comb = data.frame(trans.methods = rep(rep(trans.list, each=length(x.list)), 3),
                         est.methods = rep(c('group-glasso', 'nbd-group-lasso', 'thresholding'), each=3*length(x.list)),
                       TP = c(y.copula.g, y.cmc.g, y.pcmc.g, y.copula.nb, y.cmc.nb, y.pcmc.nb,
                              y.copula.t, y.cmc.t, y.pcmc.t),
                       FP = rep(x.list,9))
    df.comb <- rbind(df.comb,df2)
    
    df.comb.select <- df.comb[df.comb$est.methods %in% c('group-glasso', 'nbd-group-lasso'),]
    plot <- ggplot(df.comb.select, aes(x=FP, y=TP)) +
    geom_line(aes(color=trans.methods, linetype=est.methods))+ xlab("FPR") + ylab("TPR") + 
      theme_bw() + 
      scale_linetype_manual(values=c(2,1)) +
      scale_color_manual(values=c("#F8766D", "#00BA38", "#7570B3"))
      #scale_linetype_discrete(labels=c('group-glasso', 'nbd-group-lasso'))
    plot.list[[index]] = plot
    index = index + 1
  }
}

  
pdf("~/work/CMC-GGM/Plot/roc_proj.pdf")
do.call("ggarrange", c(plot.list, ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom"))   
dev.off()

pdf("~/work/CMC-GGM/Plot/roc_proj_slides.pdf")
do.call("ggarrange", c(list(plot.list[[4]], plot.list[[5]], 
                            plot.list[[7]], plot.list[[8]]),
                       ncol = 2, nrow = 2, font.label = list(size = 16, color = "black"),
                       common.legend = TRUE, legend="bottom"))   
dev.off()
