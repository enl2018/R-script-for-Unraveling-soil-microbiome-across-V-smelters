#读取 OTU 丰度表

#包含预先选择好的 50 个重要的 OTU 相对丰度以及这 90个土壤样本对应的钒含量
otu <- read.delim("E:/R/RF3/RF_Alpha_0130.txt", row.names = 1)

################
##randomForest 
library(randomForest)

#随机森林计算（默认生成 500 棵树）
set.seed(123)
otu_forest <- randomForest(Shannon~., data = otu, importance = TRUE, ntree = 500, nPerm = 1000)
otu_forest

#使用函数 importance() 查看表示每个预测变量（细菌 OTU）重要性的得分（标准化后的得分）
importance_otu.scale <- data.frame(importance(otu_forest, scale = TRUE), check.names = FALSE)
importance_otu.scale

#对预测变量（细菌 OTU）按重要性得分排个序，例如根据“%IncMSE”
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]

#简单地作图展示预测变量（细菌 OTU）的 %IncMSE 值
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
#右上角备注模型的已知解释率

p <- p +
  annotate('text', label = 'PC', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 96.14), x = 9, y = 13, size = 3, parse = TRUE)

p

################
##rfPermute 
library(rfPermute)

#使用函数 rfPermut() 重新对上述数据执行随机森林分析，详情 ?rfPermut
#rfPermut() 封装了 randomForest() 的方法，因此在给定数据和运行参数一致的情况下，两个函数结果也是一致的
#并在这里额外通过 nrep 参数执行 1000 次的随机置换以评估变量显著性的 p 值
#若数据量较大，可通过 num.cores 参数设置多线程运算
set.seed(123)
otu_rfP <- rfPermute(Shannon~., data = otu, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
otu_rfP

#提取预测变量（细菌 OTU）的重要性得分（标准化后的得分）
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale

#提取预测变量（细菌 OTU）的重要性得分的显著性（以标准化后的得分为例）
# summary(otu_rfP)
importance_otu.scale.pval <- (otu_rfP$pval)[ , , 2]
importance_otu.scale.pval

#鐢变簬鍩轰簬缃崲妫�楠岀殑鍘熺悊锛屽洜姝ゆ垜浠笉濡ㄦ煡鐪嬪彉閲忛殢鏈虹疆鎹㈠悗寰楀垎鐨勯浂鍒嗗竷鏇茬嚎锛堝浘涓伆鑹叉洸绾匡級锛屼互鍙婂疄闄呰娴嬪�肩殑寰楀垎锛堝浘涓孩绾匡級锛屼互姝ゆ瘮杈冭娴嬪�兼槸鍚﹀ぇ浜庣粷澶у鏁伴浂鍒嗗竷锛屼互鎺ㄦ柇鏄捐憲鎬?
#瑙傚療棰勬祴鍙橀噺锛堢粏鑿? OTU锛夐噸瑕佹�у緱鍒嗙殑瑙傛祴鍊煎拰闆跺垎甯冿紙浠ユ爣鍑嗗寲鍚庣殑寰楀垎涓轰緥锛?
#plotNull() 涓�氳繃浣跨敤 perds 鍙傛暟鎸囧畾棰勬祴鍙橀噺鐨勫悕绉帮紝渚嬪
plotNull(otu_rfP, scale = TRUE, preds = 'OTU_2') 
plotNull(otu_rfP, scale = TRUE, preds = 'OTU_370') 

#作图展示预测变量（细菌 OTU）的重要性得分（标准化后的得分），其中显著的得分（默认 p<0.05）以红色显示
plot(rp.importance(otu_rfP, scale = TRUE))

#瀵归娴嬪彉閲忥紙缁嗚弻 OTU锛夋寜閲嶈鎬у緱鍒嗘帓涓簭锛屼緥濡傛牴鎹�?%IncMSE鈥?
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]

#绠�鍗曞湴浣滃浘灞曠ず棰勬祴鍙橀噺锛堢粏鑿? OTU锛夌殑 %IncMSE 鍊?
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

#鏍囪棰勬祴鍙橀噺锛堢粏鑿? OTU锛夌殑鏄捐憲鎬т俊鎭?
#榛樿浠? p<0.05 涓? *锛宲<0.01 涓? **锛宲<0.001 涓? ***
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

#鍙充笂瑙掑娉ㄦā鍨嬬殑宸茬煡瑙ｉ噴鐜?
p <- p +
  annotate('text', label = 'Plant Age', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 96.14), x = 9, y = 13, size = 3, parse = TRUE)

p

################
##A3 鍖呭鐢ㄩ殢鏈烘．鏋楋紝浠ヨ瘎浼版ā鍨? p 鍊?
library(A3)

#model.fn=randomForest 璋冪敤闅忔満妫灄鐨勬柟娉曡繘琛岃繍绠?
#p.acc=0.001 琛ㄧず鍩轰簬 1000 娆￠殢鏈虹疆鎹㈣幏寰楀 p 鍊肩殑浼拌锛宲.acc 鍊艰秺灏忎唬琛ㄧ疆鎹㈡鏁拌秺澶氾紝杩愮畻涔熷氨瓒婃參锛屽洜姝ゅ鏋滃鍏ㄦā鍨? p 鍊间笉鏄緢杩垏鐨勮瘽杩樻槸鎱庣敤
#model.args 鐢ㄤ簬浼犻�掑弬鏁扮粰 randomForest()锛屽洜姝ら噷闈㈢殑鍙傛暟椤规牴鎹? randomForest() 鐨勫弬鏁伴」鑰屽畾锛屽叿浣撳彲 ?randomForest
#鍏跺畠璇︽儏鍙? ?a3 鏌ョ湅甯姪 
set.seed(123)
otu_forest.pval <- a3(plant_age~., data = otu, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 500))
otu_forest.pval

#缁х画鍦ㄥ彸涓婅澶囨敞鍏ㄦā鍨嬬殑鏄捐憲鎬т俊鎭?
p <- p +
  annotate('text', label = sprintf('italic(P) < %.3f', 0.001), x = 9, y = 12, size = 3, parse = TRUE)

p
