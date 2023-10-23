workdir="C:/Users/yww/Desktop/Vertical/NTI"
setwd(workdir)
#首先加载上述自定义函数 raup_crick()，不多说
source('RaupCrick.txt')

#然后读取示例的分类群丰度表
#上述自定义函数 raup_crick() 要求输入矩阵里面行是样本，列是物种或分类群
spp <- read.delim('Shallow.txt', sep = '\t')

#使用上述自定义函数 raup_crick() 计算所有样本（群落）对之间的 Raup-Crick 相异指数
#详细参数设置，请审阅源代码里面的注释
set.seed(123)
raup_crick.dist <- raup_crick(spXsite = spp, 
    plot_names_in_col1 = TRUE, 
    classic_metric = FALSE, 
    split_ties = TRUE, 
    reps = 999, 
    set_all_species_equal = FALSE, 
    as.distance.matrix = TRUE, 
    report_similarity = FALSE
)

#计算好的 Raup-Crick 相异指数默认以 dist 类型存储，可转换为 matrix 类型输出
raup_crick.matrix <- as.matrix(raup_crick.dist)
write.table(raup_crick.matrix, 'Raup-Crick-Shallow.txt', sep = '\t', col.names = NA)

##基于 Raup-Crick 指数的 NMDS 分析的例子
library(vegan)
library(ggplot2)

#读取计算好的 Raup-Crick 相异指数矩阵
beta_RC <- read.delim('Raup-Crick.txt', row.names = 1, sep = '\t')
beta_RC.dist <- as.dist(beta_RC)	#转为 dist 数据类型

#NMDS 排序，定义 2 个维度，详情 ?metaMDS
nmds_dis <- metaMDS(beta_RC.dist, k = 2)
nmds_dis

#获取 stress 值
stress <- nmds_dis$stress
stress

#提取样本得分（排序坐标）
nmds_dis_site <- data.frame(nmds_dis$points)
nmds_dis_site

#添加分组信息，绘制 NMDS 图
nmds_dis_site$group <- c(rep('ENV1', 10), rep('ENV2', 10))

ggplot(data = nmds_dis_site, aes(MDS1, MDS2)) +
geom_point(aes(color = group)) +
scale_color_manual(values = c('red3', 'orange3')) +
scale_fill_manual(values = c('red', 'orange')) +
stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +	
theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent')) +
labs(x = 'NMDS1', y = 'NMDS2') +
annotate('text', label = paste('Stress =', round(stress, 3)), x = 0.35, y = 0.4, size = 4, colour = 'black') 

##再象征性画个 Raup-Crick 相异指数的组内箱线图
env1 <- as.dist(beta_RC[1:10,1:10])
env2 <- as.dist(beta_RC[11:20,11:20])

plot_data <- data.frame(
    beta_RC = c(env1, env2), 
    group = c(rep('ENV1', length(env1)), rep('ENV2', length(env2)))
)

ggplot(plot_data, aes(group, beta_RC, fill = group)) +
geom_boxplot(show.legend = FALSE) +
scale_fill_manual(values = c('red', 'orange')) +
labs(x = '', y = 'Raup-Crick Dissimilarity')
