library(vegan)
#读取上述数据集
setwd("E:/R/Mantel test 7")
df<-read.csv('OTU4.csv',header=TRUE)
head(df,3)
dim(df)

#根据物种丰度数据，计算样方间的 Bray-curtis 距离
abund <- df[ ,21:ncol(df)]
dist.abund <- vegdist(abund, method = 'bray')
#根据环境测量指标，计算样方间的欧几里得距离
#钒
V <- df$V
dist.V <- dist(V, method = 'euclidean')
#有效钒
Va <- df$Va
dist.Va <- dist(Va, method = 'euclidean')
#OM-V
VO <- df$VO
dist.VO <- dist(VO, method = 'euclidean')
#OM-V%
VOF <- df$VOF
dist.VOF <- dist(VOF, method = 'euclidean')
#M-V
VM <- df$VM
dist.VM <- dist(VM, method = 'euclidean')
#M-V%
VMF <- df$VMF
dist.VMF <- dist(VMF, method = 'euclidean')
#重金属

#pH
pH <- df$pH
dist.pH <- dist(pH, method = 'euclidean')
#OM
OM <- df$OM
dist.OM <- dist(OM, method = 'euclidean')
#TN
TN <- df$TN
dist.TN <- dist(TN, method = 'euclidean')
#AP
AP <- df$AP
dist.AP <- dist(AP, method = 'euclidean')
#AS
AS <- df$AS
dist.AS <- dist(AS, method = 'euclidean')
#Shannon
Shannon <- df$Shannon
dist.Shannon <- dist(Shannon, method = 'euclidean')
#Chao1
Chao1 <- df$Chao1
dist.Chao1 <- dist(Chao1, method = 'euclidean')
#GD
GD <- df$GD
dist.GD <- dist(GD, method = 'euclidean')
#Lat
Lat <- df$Lat
dist.Lat <- dist(Lat, method = 'euclidean')
#Long
Long <- df$Long
dist.Long <- dist(Long, method = 'euclidean')
#Temperature
MAT <- df$MAT
dist.MAT <- dist(MAT, method = 'euclidean')
#Precipitation
PL <- df$PL
dist.PL <- dist(PL, method = 'euclidean')


#如果期望关注多种环境的协同作用，就选择一个环境子集，计算样方间的欧几里得距离
env <- df[ ,9:13]
scale.env <- scale(env, center = TRUE, scale = TRUE)
dist.env <- dist(scale.env, method = 'euclidean')

#如果期望关注多种环境的协同作用，就选择一个环境子集，计算样方间的欧几里得距离
climate <- df[ ,19:20]
scale.climate <- scale(climate, center = TRUE, scale = TRUE)
dist.climate <- dist(scale.climate, method = 'euclidean')

install.packages("geosphere")
#根据经纬度，计算样方间实际的地理距离
library(geosphere)
geo <- data.frame(df$Long, df$Lat)
d.geo <- distm(geo, fun = distHaversine)      
dist.geo <- as.dist(d.geo)

##执行 Mantel tests
#物种丰度和V的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_V <- mantel(dist.abund, dist.V, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_V
#物种丰度和VM的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_VM <- mantel(dist.abund, dist.VM, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VM
#物种丰度和VMF的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_VMF <- mantel(dist.abund, dist.VMF, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VMF
#物种丰度和VO的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_VO <- mantel(dist.abund, dist.VO, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VO
#物种丰度和VOF的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_VOF <- mantel(dist.abund, dist.pH, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VOF
#物种丰度和pH的相关性，以 spearman 相关系数为例，9999 次置换检验显著性（Mantel tests 基于随机置换的方法获取 p 值）
abund_pH <- mantel(dist.abund, dist.pH, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_pH
#物种丰度和Available V的相关性
abund_Va <- mantel(dist.abund, dist.Va, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Va
#物种丰度和OM的相关性
abund_OM <- mantel(dist.abund, dist.OM, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_OM
#物种丰度和TN的相关性
abund_TN <- mantel(dist.abund, dist.TN, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_TN
#物种丰度和AP的相关性
abund_AP <- mantel(dist.abund, dist.AP, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_AP
#物种丰度和AS的相关性
abund_AS <- mantel(dist.abund, dist.AS, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_AS
#物种丰度和Shannon的相关性
abund_Shannon <- mantel(dist.abund, dist.Shannon, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Shannon
#物种丰度和Chao1的相关性
abund_Chao1 <- mantel(dist.abund, dist.Chao1, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Chao1
#物种丰度和Lat的相关性
abund_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Lat
#物种丰度和Long的相关性
abund_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Long
#物种丰度和环境组合的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
abund_env <- mantel(dist.abund, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_env
PL_env <- mantel(dist.PL, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
PL_env
MAT_env <- mantel(dist.MAT, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
MAT_env
#物种丰度和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
abund_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_geo
Shannon_geo <- mantel(dist.abund, dist.Shannon, method = 'spearman', permutations = 9999, na.rm = TRUE)
Shannon_geo
abund_PL <- mantel(dist.abund, dist.PL, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_PL
abund_MAT <- mantel(dist.abund, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_MAT
#环境组合和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
env_geo <- mantel(dist.env, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_geo
#TN和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
TN_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_geo
TN_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_Lat
TN_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_Long
#AP和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
AP_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_geo
AP_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_Lat
AP_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_Long
#OM和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
OM_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_geo
OM_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_Lat
OM_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_Long
#AS和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
AS_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_geo
AS_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_Lat
AS_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_Long
#pH和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
pH_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_geo
pH_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_Lat
pH_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_Long
#V和地理距离的相关性，以 spearman 相关系数为例，9999 次置换检验显著性
V_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_geo
V_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_Lat
V_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_Long
#物种丰度和Temperature的相关性
abund_MAT <- mantel(dist.abund, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_MAT
#物种丰度和PL的相关性
abund_PL <- mantel(dist.abund, dist.PL, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_PL
#地理和PL的相关性
PL_geo <- mantel(dist.PL, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
PL_geo
#地理和V的相关性
V_geo <- mantel(dist.V, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_geo
#随便
env_geo <- mantel(dist.enV, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_geo
env_MAT <- mantel(dist.enV, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_MAT
#某物种与温度的相关性，横轴温度，纵轴物种丰度，颜色表示样方的纬度
library(ggplot2)
#Sphingomonas vs V
xx1 = ggplot(df, aes(x = V, y = Sphingomonas)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'Sphingonomas', x = 'V (mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx1

library(ggplot2)

#Actinobacteria vs V
xx2 = ggplot(df, aes(x = V, y = norank_c_Actinobacteria)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'norank_c_Actinobacteria', x = 'V (mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx2

#Microvirga vs TN
xx3 = ggplot(df, aes(x = TN, y = Microvirga)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'Microvirga', x = 'TN (mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'Dark green', low = 'green')
xx3

#Bacillus vs AP
xx4 = ggplot(df, aes(x = AP, y = Bacillus)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'Bacillus', x = 'Ln(AP) (mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'Brown', low = 'red')
xx4

#Bradyrhizobium vs AS
xx5 = ggplot(df, aes(x = AS, y = Bradyrhizobium)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'Bradyrhizobium', x = 'Ln(AS)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'Dark blue', low = 'light blue')
xx5

#Solirubrobacter vs V
xx6 = ggplot(df, aes(x = V, y = Solirubrobacter)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'Solirubrobacter', x = 'V(mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'Dark blue', low = 'light blue')
xx6

#某物种与温度的相关性，横轴温度，纵轴物种丰度，颜色表示样方的纬度
library(ggplot2)
#norank_c_Actinobacteria vs AV
xx7 = ggplot(df, aes(x = AV, y = norank_c_Actinobacteria)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'norank_c_Actinobacteria', x = 'AV Fraction (%)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx7

#norank_c_Actinobacteria vs AV2
xx8 = ggplot(df, aes(x = AV2, y = norank_c_Actinobacteria)) +
  
  geom_smooth(method = 'lm', alpha = 0.2, colour = 'black') +
  
  geom_point(aes(colour = Lat), size = 4) +
  
  labs(y = 'norank_c_Actinobacteria', x = 'Metal associateed V (mg/kg)') +
  
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12),
         
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'),
         
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = 'black'),
         
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  
  scale_colour_continuous(high = 'navy', low = 'salmon')
xx8

#基于物种丰度的距离与基于温度指标的距离之间的相关性散点图，上文已知二者显著相关；同时颜色表示样方间地理距离
library(ggplot2)
mm <- ggplot(mat, aes(y = aa, x = tt)) +
  
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = gg/1000)) +
  
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  
  labs(x = "V (mg/kg)", y = "Bray-Curtis Dissimilarity", fill = "Physical Separation (km)") +
  
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12),
         
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
         
         axis.title= element_text(face = "bold", size = 14, colour = "black"),
         
         panel.background = element_blank(),
         
         panel.border = element_rect(fill = NA, colour = "black"),
         
         legend.position = "top",
         
         legend.text = element_text(size = 10, face = "bold"),
         
         legend.title = element_text(size = 11, face = "bold")) +
  
  scale_fill_continuous(high = "navy", low = "skyblue")


mm
