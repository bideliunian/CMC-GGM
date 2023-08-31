####################
## plots showing the performance of cm transformation
##########################
function.path <- "~/work/CMC-GGM/Functions"
working_path <- "~/work/CMC-GGM/Plot"
import_path <- "~/work/CMC-GGM/proj_cmc"

# source all function scipts from the function path
function.sources <- list.files(function.path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function.sources, source, .GlobalEnv)


library(ggplot2)
library(gridExtra)

##########################
### CMC vs PCMC d = 2
##########################

####### parameters
n = 400
p = 2
M = 2
data = gen_data(n = n, p = p, m = M, model = 'model1', trans.type = 'cmcopula', 
                run.ind = 2022, trans.func = 'power')
data.pre = data$X.pre
data.X = data$X

data.copula = npn(x = data.X[,1:M])
data.cmcopula = cmcopula(x = data.X[,1:M], fun = 'cm', grid.method = 'random', truncation = TRUE)

df1 = cbind(data.pre[1:n,1:2], data.X[1:n, 1:2], data.copula[1:n, 1:2], data.cmcopula[1:n, 1:2])
df1 = df1[!rowSums(abs(df1) > 2.5),][1:150, ]
colnames(df1) = rep(c('x1','x2'), 4)
df = as.data.frame(rbind(df1[,1:2], df1[,3:4], df1[,5:6], df1[,7:8]))
names(df) = c('x1','x2')
df$trans.type = rep(c('oracle', 'linear', 'copula', 'cmcopula'), each = 150)
df$lable = rep(seq(1:150), times = 4)

df.linear = df[df$trans.type %in% c('oracle','linear'),]
p1 = ggplot(df.linear, aes(x1, x2)) +
  geom_line(aes(group=lable), linetype=2, linewidth=0.15) +
  geom_point(aes(color=trans.type), size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

df.copula = df[df$trans.type %in% c('oracle','copula'),]
p2 = ggplot(df.copula, aes(x1, x2)) +
  geom_line(aes(group=lable), linetype=2, linewidth=0.15) +
  geom_point(aes(color=trans.type), size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

df.cmcopula = df[df$trans.type %in% c('oracle','cmcopula'),]
p3 = ggplot(df.cmcopula, aes(x1, x2)) +
  geom_line(aes(group=lable), linetype=2, linewidth=0.15) +
  geom_point(aes(color=trans.type), size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

pdf("~/work/CMC-GGM/Plot/trans_power.pdf")
ggarrange(p1, ggplot() + theme_void(), p2, ggplot() + theme_void(), p3, 
          ncol = 5, nrow = 3, widths = c(1, -0.02, 1, -0.02, 1),
          legend = "none", heights = 1.5)
dev.off()

