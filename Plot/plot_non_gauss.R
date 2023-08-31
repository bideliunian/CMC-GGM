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
### CMC vs Copula d = 2
##########################

####### parameters
n = 300
p = 2
M = 2

#######################
data.exp = gen_data(n = n, p = p, m = M, model = 'model1', trans.type = 'cmcopula', 
                run.ind = 2022, trans.func = 'exp')
data.X.exp = data.exp$X
data.copula.exp = npn(x = data.X.exp[,1:M])
df.exp = as.data.frame(data.copula.exp[!rowSums(abs(data.copula.exp) > 2.5),][1:150, ])
names(df.exp) = c('x1','x2')


################### 
data.power = gen_data(n = n, p = p, m = M, model = 'model1', trans.type = 'cmcopula', 
                run.ind = 2022, trans.func = 'power')
data.X.power = data.power$X
data.copula.power = npn(x = data.X.power[,1:M])
df.power = as.data.frame(data.copula.power[!rowSums(abs(data.copula.power) > 2.5),][1:150, ])
names(df.power) = c('x1','x2')
###################


data.power1 = gen_data(n = n, p = p, m = M, model = 'model1', trans.type = 'pcmcopula1', 
                      run.ind = 2022, trans.func = 'power')
data.X.power1 = data.power1$X
data.copula.power1 = npn(x = data.X.power1[,1:M])
df.power1 = as.data.frame(data.copula.power1[!rowSums(abs(data.copula.power1) > 2.5),][1:150, ])
names(df.power1) = c('x1','x2')
###################


data.exp1 = gen_data(n = n, p = p, m = M, model = 'model1', trans.type = 'pcmcopula1', 
                       run.ind = 2022, trans.func = 'exp')
data.X.exp1 = data.exp1$X
data.copula.exp1 = npn(x = data.X.exp1[,1:M])
df.exp1 = as.data.frame(data.copula.exp1[!rowSums(abs(data.copula.exp1) > 2.5),][1:150, ])
names(df.exp1) = c('x1','x2')
###################

p1 = ggplot(df.power, aes(x1, x2)) +
  geom_point(size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") + 
  theme(plot.title = element_text(size=10, hjust = 0.5))

p2 = ggplot(df.exp, aes(x2, x1)) +
  geom_point(size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

p3 = ggplot(df.power1, aes(x1, x2)) +
  geom_point(size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

p4 = ggplot(df.exp1, aes(x1, x2)) +
  geom_point(size=0.3) + 
  xlim(-2.5, 2.5) + ylim(-2.5, 2.5) + 
  theme_bw() + 
  labs(x ="", y = "") +
  theme(plot.title = element_text(size=10, hjust = 0.5))


pdf("~/work/CMC-GGM/Plot/non_gaussian.pdf")
ggarrange(p1, ggplot() + theme_void(), p2, ggplot() + theme_void(), p4, 
          ncol = 5, nrow = 3, widths = c(1, -0.02, 1, -0.02, 1),
          legend = "none", heights = 1.5)
dev.off()

