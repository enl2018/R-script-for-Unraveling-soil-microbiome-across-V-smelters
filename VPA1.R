library(vegan)
explain <- function(a,b){a/b}
#读入数据
setwd("E:/R/VPA")
OTU<-read.table("OTU.txt",header=T)
phyla<-read.table("Phyla.txt",header=T)
env<-read.table("Env.txt",header=T)
geo<-read.table("Geo.txt",header=T)
V<-read.table("V_data.txt",header=T)
chem<-read.table("Chem.txt",header=T)
div<-read.table("Diversity.txt",header=T)
chao<-read.table("Chao.txt",header=T)
shann<-read.table("Shannon.txt",header=T)
alpha<-read.table("Alpha.txt",header=T)
beta<-read.table("Beta.txt",header=T)
PC<-read.table("PC.txt",header=T)

# 拟合模型
?varpart # 简单了解函数用法
fit <- varpart (div,env, geo, transfo="hel") 
fit
plot(fit, bg = c("hotpink","skyblue")) # bg表示两个变量的背景颜色

# 三个
mod <- varpart(OTU, ~ SubsDens + WatrCont, ~ Substrate + Shrub + Topo,
               env, , transfo="hel")
mod
plot(mod, bg=2:4)
