library(ape)
library(picante)

workdir="E:/R/bNTI_0115"
heatmap_file="E:/R/bNTI_0115/OTU300.txt"

setwd(workdir)
tree <- read.tree("phylo_tree.tre")
otu=read.table(heatmap_file,check.names =F,comment.char ="",sep = "\t",header = T,row.names = 1)
otu <- t(otu)

#组内过度分散检验，计算MNTD和NTI指数
dis <- prune.sample(otu,tree)
dis <- cophenetic(dis)
otu <- otu[,colnames(dis)]
mntd <- ses.mntd(otu,dis,abundance.weighted = TRUE,null.model = "taxa.labels",runs = 9)
#结果输出
write.csv(ses.mntd,file = "NTI.csv")
#abundance.weighted=TRUE，这将改变对这些指标的解释，从物种间的平均系统发育距离，到个体间的平均系统发育距离�?
#ses.mpd <- ses.mpd(otu,dis, abundance.weighted=TRUE, null.model="taxa.labels",runs=999)
#pd <- pd(otu, tree, include.root=TRUE)

#接下来计算βMNTD和βNTI
mntd.obs <- comdistnt(otu, dis, abundance.weighted = TRUE)
mntd.obs <- t(as.vector(t(mntd.obs)))
mntd.rand <- replicate(99,as.vector(t(comdistnt(otu,taxaShuffle(dis),
                                                abundance.weighted = TRUE))))
mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 1, FUN = mean, 
                        na.rm = TRUE)
mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 1, FUN = sd, 
                      na.rm = TRUE)
mntd.obs.z <- as.vector((mntd.obs - mntd.rand.mean)/mntd.rand.sd)
beta.mntd <- data.frame(as.vector(mntd.obs), mntd.rand.mean, mntd.rand.sd,mntd.obs.z)

a <- rownames(otu)
X <- c()
for (i in 2:length(a)) {
  X <- c(X,a[i:length(a)])    
}

write.csv(ses.mntd,file = "NTI.csv")

write.csv(beta.mntd,file = "betaNTI.csv")
