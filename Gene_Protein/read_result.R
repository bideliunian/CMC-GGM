######## PART 1: reading data #########################
## name of paths
function.path <- "~/work/CMC-GGM/Functions"
working.path <- "~/work/CMC-GGM/Gene_Protein"
save.path <- "~/work/CMC-GGM/Gene_Protein/Result"

load(paste(working.path, "/gpconsensus.RData", sep = ""))
gene <-  scale(t(dat$gene.select), TRUE,TRUE)
prot <-  scale(t(dat$protein.select), TRUE,TRUE)
n <- nrow(gene)
p <- ncol(gene)
M <- 2

n.expr <- 50
cmc.prec.freq <- copula.prec.freq <- linear.prec.freq <- matrix(0, p, p)

for (i in 1:n.expr) {
  load(paste(save.path,"/prec.bt.seed",i,".Rdata",sep=""))
  
  cmc.prec <- prec.list[[1]]
  copula.prec <- prec.list[[2]]
  linear.prec <- prec.list[[3]]
  
  cmc.prec.freq <- cmc.prec.freq + cmc.prec / n.expr
  copula.prec.freq <- copula.prec.freq + copula.prec / n.expr
  linear.prec.freq <- linear.prec.freq + linear.prec / n.expr
}


############# PART 2: Summary statistics #########################
thres <- 0.9

robust.prec <- function(prec.freq, threshold){
  prec.robust <- prec.freq
  prec.robust[prec.robust > threshold] <- 1
  prec.robust[prec.robust < threshold] <- 0
  
  return(prec.robust)
}

cmc.prec.robust <- robust.prec(cmc.prec.freq, thres)
copula.prec.robust <- robust.prec(copula.prec.freq, thres)
linear.prec.robust <- robust.prec(linear.prec.freq, thres)

####### number of edges ####################
num.edge.cmc <- sum(cmc.prec.robust)/2
num.edge.copula <- sum(copula.prec.robust)/2
num.edge.linear <- sum(linear.prec.robust)/2

cat(paste("edge numbers: cmc--", num.edge.cmc, "; copula--", num.edge.copula,
          "; linear--", num.edge.linear))

############# number of shared edges ##################
edge.share.cmc.copula <- cmc.prec.robust & copula.prec.robust
edge.share.cmc.linear <- cmc.prec.robust & linear.prec.robust
edge.share.copula.linear <- copula.prec.robust & linear.prec.robust

cat(paste("edge numbers shared by cmc and linear", sum(edge.share.cmc.linear)/2))
cat(paste("edge numbers shared by cmc and copula", sum(edge.share.cmc.copula)/2))
cat(paste("edge numbers shared by copula and linear", sum(edge.share.copula.linear)/2))

################ avg node degree ########################
avg.node <- function(prec){
  return(round(mean(rowSums(prec)), 3))
}
avg.node.cmc <- avg.node(cmc.prec.robust)
avg.node.copula <- avg.node(copula.prec.robust)
avg.node.linear <- avg.node(linear.prec.robust)

cat(paste("average node degree: cmc--", avg.node.cmc, "; copula--", avg.node.copula,
          "; linear--", avg.node.linear ))


############ PART 3: plotting ###################################
graph.cmc <- as.data.frame(cmc.prec.robust + edge.share.cmc.linear)
graph.linear <- as.data.frame(linear.prec.robust + edge.share.cmc.linear)

rownames(graph.cmc) <- colnames(graph.cmc) <- as.character(c(1:p))
graph.cmc$id <- as.character(c(1:p))
df.cmc <- melt(graph.cmc, id='id')

rownames(graph.linear) <- colnames(graph.linear) <- as.character(c(1:p))
graph.linear$id <- as.character(c(1:p))
df.linear <- melt(graph.linear, id='id')

library(reshape2)
library(ggplot2)

## using default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = c('#FFFFFF', gg_color_hue(3))

# "#F8766D" "#00BA38" "#619CFF"

png(paste("~/work/CMC-GGM/Gene_Protein/htmap_cmc.png"))
ggplot(df.cmc, 
       aes(x = variable, y = id, fill = factor(value))) + 
  geom_tile() + 
  scale_fill_manual(values=c("#FFFFFF","#619CFF","#F8766D")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none") +
  scale_x_discrete(limits = graph.cmc$id) +
  scale_y_discrete(limits = rev(graph.cmc$id)) 

dev.off()

png(paste("~/work/CMC-GGM/Gene_Protein/htmap_linear.png"))
ggplot(df.linear, 
       aes(x = variable, y = id, fill = factor(value))) + 
  geom_tile() + 
  scale_fill_manual(values=c("#FFFFFF","#00BA38","#F8766D")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none") +
  scale_x_discrete(limits = graph.linear$id) +
  scale_y_discrete(limits = rev(graph.linear$id)) 
dev.off()

