#��ȡ OTU ��ȱ�

#����Ԥ��ѡ��õ� 50 ����Ҫ�� OTU ��Է���Լ��� 90������������Ӧ�ķ�����
otu <- read.delim("E:/R/RF3/RF_Alpha_0130.txt", row.names = 1)

################
##randomForest 
library(randomForest)

#���ɭ�ּ��㣨Ĭ������ 500 ������
set.seed(123)
otu_forest <- randomForest(Shannon~., data = otu, importance = TRUE, ntree = 500, nPerm = 1000)
otu_forest

#ʹ�ú��� importance() �鿴��ʾÿ��Ԥ�������ϸ�� OTU����Ҫ�Եĵ÷֣���׼����ĵ÷֣�
importance_otu.scale <- data.frame(importance(otu_forest, scale = TRUE), check.names = FALSE)
importance_otu.scale

#��Ԥ�������ϸ�� OTU������Ҫ�Ե÷��Ÿ���������ݡ�%IncMSE��
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]

#�򵥵���ͼչʾԤ�������ϸ�� OTU���� %IncMSE ֵ
library(ggplot2)

importance_otu.scale$OTU_name <- rownames(importance_otu.scale)
importance_otu.scale$OTU_name <- factor(importance_otu.scale$OTU_name, levels = importance_otu.scale$OTU_name)

p <- ggplot(importance_otu.scale, aes(OTU_name, `%IncMSE`)) +
  geom_col(width = 0.5, fill = '#FFC068', color = NA) +
  labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 16))

p
#���ϽǱ�עģ�͵���֪������

p <- p +
  annotate('text', label = 'PC', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 96.14), x = 9, y = 13, size = 3, parse = TRUE)

p

################
##rfPermute 
library(rfPermute)

#ʹ�ú��� rfPermut() ���¶���������ִ�����ɭ�ַ��������� ?rfPermut
#rfPermut() ��װ�� randomForest() �ķ���������ڸ������ݺ����в���һ�µ�����£������������Ҳ��һ�µ�
#�����������ͨ�� nrep ����ִ�� 1000 �ε�����û����������������Ե� p ֵ
#���������ϴ󣬿�ͨ�� num.cores �������ö��߳�����
set.seed(123)
otu_rfP <- rfPermute(Shannon~., data = otu, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
otu_rfP

#��ȡԤ�������ϸ�� OTU������Ҫ�Ե÷֣���׼����ĵ÷֣�
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale

#��ȡԤ�������ϸ�� OTU������Ҫ�Ե÷ֵ������ԣ��Ա�׼����ĵ÷�Ϊ����
# summary(otu_rfP)
importance_otu.scale.pval <- (otu_rfP$pval)[ , , 2]
importance_otu.scale.pval

#由于基于置换检验的原理，因此我们不妨查看变量随机置换后得分的零分布曲线（图中灰色曲线），以及实际观测值的得分（图中红线），以此比较观测值是否大于绝大多数零分布，以推断显著�?
#观察预测变量（细�? OTU）重要性得分的观测值和零分布（以标准化后的得分为例�?
#plotNull() 中通过使用 perds 参数指定预测变量的名称，例如
plotNull(otu_rfP, scale = TRUE, preds = 'OTU_2') 
plotNull(otu_rfP, scale = TRUE, preds = 'OTU_370') 

#��ͼչʾԤ�������ϸ�� OTU������Ҫ�Ե÷֣���׼����ĵ÷֣������������ĵ÷֣�Ĭ�� p<0.05���Ժ�ɫ��ʾ
plot(rp.importance(otu_rfP, scale = TRUE))

#对预测变量（细菌 OTU）按重要性得分排个序，例如根据�?%IncMSE�?
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]

#简单地作图展示预测变量（细�? OTU）的 %IncMSE �?
library(ggplot2)

importance_otu.scale$OTU_name <- rownames(importance_otu.scale)
importance_otu.scale$OTU_name <- factor(importance_otu.scale$OTU_name, levels = importance_otu.scale$OTU_name)

p <- ggplot() +
  geom_col(data = importance_otu.scale, aes(x = OTU_name, y = `%IncMSE`), width = 0.5, fill = '#FFC068', color = NA) +
  labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 16))

p

#标记预测变量（细�? OTU）的显著性信�?
#默认�? p<0.05 �? *，p<0.01 �? **，p<0.001 �? ***
for (OTU in rownames(importance_otu.scale)) {
  importance_otu.scale[OTU,'%IncMSE.pval'] <- importance_otu.scale.pval[OTU,'%IncMSE']
  if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.05) importance_otu.scale[OTU,'%IncMSE.sig'] <- ''
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.01 & importance_otu.scale[OTU,'%IncMSE.pval'] < 0.05) importance_otu.scale[OTU,'%IncMSE.sig'] <- '*'
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] >= 0.001 & importance_otu.scale[OTU,'%IncMSE.pval'] < 0.01) importance_otu.scale[OTU,'%IncMSE.sig'] <- '**'
  else if (importance_otu.scale[OTU,'%IncMSE.pval'] < 0.001) importance_otu.scale[OTU,'%IncMSE.sig'] <- '***'
}

p <- p +
  geom_text(data = importance_otu.scale, aes(x = OTU_name, y = `%IncMSE`, label = `%IncMSE.sig`), nudge_y = 1)

p

#右上角备注模型的已知解释�?
p <- p +
  annotate('text', label = 'Plant Age', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 96.14), x = 9, y = 13, size = 3, parse = TRUE)

p

################
##A3 包套用随机森林，以评估模�? p �?
library(A3)

#model.fn=randomForest 调用随机森林的方法进行运�?
#p.acc=0.001 表示基于 1000 次随机置换获得对 p 值的估计，p.acc 值越小代表置换次数越多，运算也就越慢，因此如果对全模�? p 值不是很迫切的话还是慎用
#model.args 用于传递参数给 randomForest()，因此里面的参数项根�? randomForest() 的参数项而定，具体可 ?randomForest
#其它详情�? ?a3 查看帮助 
set.seed(123)
otu_forest.pval <- a3(plant_age~., data = otu, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
otu_forest.pval

#继续在右上角备注全模型的显著性信�?
p <- p +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.001), x = 9, y = 12, size = 3, parse = TRUE)

p
