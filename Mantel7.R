library(vegan)
#��ȡ�������ݼ�
setwd("E:/R/Mantel test 7")
df<-read.csv('OTU4.csv',header=TRUE)
head(df,3)
dim(df)

#�������ַ�����ݣ������������ Bray-curtis ����
abund <- df[ ,21:ncol(df)]
dist.abund <- vegdist(abund, method = 'bray')
#���ݻ�������ָ�꣬�����������ŷ����þ���
#��
V <- df$V
dist.V <- dist(V, method = 'euclidean')
#��Ч��
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
#�ؽ���

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


#���������ע���ֻ�����Эͬ���ã���ѡ��һ�������Ӽ��������������ŷ����þ���
env <- df[ ,9:13]
scale.env <- scale(env, center = TRUE, scale = TRUE)
dist.env <- dist(scale.env, method = 'euclidean')

#���������ע���ֻ�����Эͬ���ã���ѡ��һ�������Ӽ��������������ŷ����þ���
climate <- df[ ,19:20]
scale.climate <- scale(climate, center = TRUE, scale = TRUE)
dist.climate <- dist(scale.climate, method = 'euclidean')

install.packages("geosphere")
#���ݾ�γ�ȣ�����������ʵ�ʵĵ�������
library(geosphere)
geo <- data.frame(df$Long, df$Lat)
d.geo <- distm(geo, fun = distHaversine)      
dist.geo <- as.dist(d.geo)

##ִ�� Mantel tests
#���ַ�Ⱥ�V������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_V <- mantel(dist.abund, dist.V, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_V
#���ַ�Ⱥ�VM������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_VM <- mantel(dist.abund, dist.VM, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VM
#���ַ�Ⱥ�VMF������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_VMF <- mantel(dist.abund, dist.VMF, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VMF
#���ַ�Ⱥ�VO������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_VO <- mantel(dist.abund, dist.VO, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VO
#���ַ�Ⱥ�VOF������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_VOF <- mantel(dist.abund, dist.pH, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_VOF
#���ַ�Ⱥ�pH������ԣ��� spearman ���ϵ��Ϊ����9999 ���û����������ԣ�Mantel tests ��������û��ķ�����ȡ p ֵ��
abund_pH <- mantel(dist.abund, dist.pH, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_pH
#���ַ�Ⱥ�Available V�������
abund_Va <- mantel(dist.abund, dist.Va, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Va
#���ַ�Ⱥ�OM�������
abund_OM <- mantel(dist.abund, dist.OM, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_OM
#���ַ�Ⱥ�TN�������
abund_TN <- mantel(dist.abund, dist.TN, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_TN
#���ַ�Ⱥ�AP�������
abund_AP <- mantel(dist.abund, dist.AP, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_AP
#���ַ�Ⱥ�AS�������
abund_AS <- mantel(dist.abund, dist.AS, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_AS
#���ַ�Ⱥ�Shannon�������
abund_Shannon <- mantel(dist.abund, dist.Shannon, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Shannon
#���ַ�Ⱥ�Chao1�������
abund_Chao1 <- mantel(dist.abund, dist.Chao1, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Chao1
#���ַ�Ⱥ�Lat�������
abund_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Lat
#���ַ�Ⱥ�Long�������
abund_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_Long
#���ַ�Ⱥͻ�����ϵ�����ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
abund_env <- mantel(dist.abund, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_env
PL_env <- mantel(dist.PL, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
PL_env
MAT_env <- mantel(dist.MAT, dist.env, method = 'spearman', permutations = 9999, na.rm = TRUE)
MAT_env
#���ַ�Ⱥ͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
abund_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_geo
Shannon_geo <- mantel(dist.abund, dist.Shannon, method = 'spearman', permutations = 9999, na.rm = TRUE)
Shannon_geo
abund_PL <- mantel(dist.abund, dist.PL, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_PL
abund_MAT <- mantel(dist.abund, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_MAT
#������Ϻ͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
env_geo <- mantel(dist.env, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_geo
#TN�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
TN_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_geo
TN_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_Lat
TN_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
TN_Long
#AP�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
AP_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_geo
AP_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_Lat
AP_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
AP_Long
#OM�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
OM_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_geo
OM_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_Lat
OM_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
OM_Long
#AS�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
AS_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_geo
AS_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_Lat
AS_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
AS_Long
#pH�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
pH_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_geo
pH_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_Lat
pH_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
pH_Long
#V�͵������������ԣ��� spearman ���ϵ��Ϊ����9999 ���û�����������
V_geo <- mantel(dist.abund, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_geo
V_Lat <- mantel(dist.abund, dist.Lat, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_Lat
V_Long <- mantel(dist.abund, dist.Long, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_Long
#���ַ�Ⱥ�Temperature�������
abund_MAT <- mantel(dist.abund, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_MAT
#���ַ�Ⱥ�PL�������
abund_PL <- mantel(dist.abund, dist.PL, method = 'spearman', permutations = 9999, na.rm = TRUE)
abund_PL
#������PL�������
PL_geo <- mantel(dist.PL, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
PL_geo
#������V�������
V_geo <- mantel(dist.V, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
V_geo
#���
env_geo <- mantel(dist.enV, dist.geo, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_geo
env_MAT <- mantel(dist.enV, dist.MAT, method = 'spearman', permutations = 9999, na.rm = TRUE)
env_MAT
#ĳ�������¶ȵ�����ԣ������¶ȣ��������ַ�ȣ���ɫ��ʾ������γ��
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

#ĳ�������¶ȵ�����ԣ������¶ȣ��������ַ�ȣ���ɫ��ʾ������γ��
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

#�������ַ�ȵľ���������¶�ָ��ľ���֮��������ɢ��ͼ��������֪����������أ�ͬʱ��ɫ��ʾ�������������
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