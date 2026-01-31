##Script to Figure 2-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("vegan","ggplot2","tidyverse","ggExtra")
bioconductor_packages=c()

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

if (!require("BiocManager", quietly = T, warn.conflicts = F)) install.packages("BiocManager")

for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}


############################## Data loading ##############################
otu_raw <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2B_Bacts.csv")
group <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2_sample.csv")#读入分组文件

############################## Data processing #####################################
otu_raw$OTU <- gsub(".*s__(\\S+).*", "s__\\1", otu_raw$OTU)
rownames(otu_raw) <- otu_raw$OTU
otu_raw <- otu_raw[,-1]
otu <- t(otu_raw) #PCA分析的目标是探索样本间的差异，因此每个样本应作为一行，使样本成为观测单元。
#计算bray_curtis距离
otu_distance <- vegdist(otu, method = 'bray') #可选euclidean、manhattan、jaccard
bray <- as.matrix(otu_distance) #转成距离矩阵
#pcoa分析
pcoa <- cmdscale(bray,eig=TRUE) #k=2 表示降维到2维，可选k=(nrow(data)-1
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度:eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
#pcl2原来是matrix,转化为data.frame
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$sample <- row.names(pc12)
#将绘图数据和分组合并
df <- merge(pc12,group,by="sample")
#浅看散点图
ggplot(df,aes(x=V1, y=V2, color=treatment))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
#置换多元（因素）方差分析（PERMANOVA） 基于bray-curtis距离进行
otu_adonis <- adonis2(bray ~ treatment, data = df, permutations = 999) #这里用距离矩阵 就不用再加method = 'bray' 能够识别出矩阵，如果用丰度表就要加了【本质上是用距离矩阵在做】
#otu_adonis1 <- adonis2(bray ~ treatment*time, data = df, permutations = 999)
##处理对菌群组成影响显著 该分析独立于PCA PCOA NMDS等降维方法 因此可以采用多种降维方法选择好看的画图 标注显著性即可
# 格式化 PERMANOVA 结果，用于图中显示
adonis_text <- sprintf("PERMANOVA:\ndf = %d\nR² = %.5f\np-value = %s",
                       otu_adonis$Df[1], otu_adonis$R2[1], 
                       format.pval(otu_adonis$`Pr(>F)`[1], digits =5, eps =0.00001)) 
adonis_text <- sprintf("PERMANOVA:\nR² = %.5f\np-value = %s",
                       otu_adonis$R2[1], 
                       format.pval(otu_adonis$`Pr(>F)`[1], digits =4, eps =0.00001)) 

####【使用】NMDS（非度量多维标度分析）
#计算 Bray-Curtis 距离的 NMDS
nmds <- metaMDS(bray, k=2)
#提取NMDS分析的stress值（应力函数值），一般不大于0.2为合理
stress <- round(nmds$stress,4) #越小表示 NMDS 结果越可靠
#NMDS评估，拟合R2越大越合理
stressplot(nmds) #检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
# 提取 NMDS 计算出的坐标，并与样本分组信息合并
plot <- nmds$points %>% as.data.frame() %>% #样方得分 转换为数据框
  rownames_to_column(var="sample") %>% # 添加样本 ID 列
  left_join(., group, by="sample") # 关联样本分组信息


############################## [Fig2-B] Painting #####################################
# 绘制结果【椭圆填充色】
color=c("#8BBCD6","#ACD26A","#DD7F88")
p2 <- ggplot(data=plot,aes(x=MDS1,y=MDS2,
                   color=treatment))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=1.8, shape =16)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", colour="#BEBEBE")+
  geom_hline(yintercept = 0,lty="dashed", colour="#BEBEBE")+#图中虚线
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+#添加数据点的标签
  # guides(color=guide_legend(title=NULL))+#去除图例标题
  stat_ellipse(data=plot,
               geom = "polygon",level=0.95,
               size= NA,
               aes(fill=treatment),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = color) +#点的颜色设置
  scale_fill_manual(values = c("#8BBCD6","#ACD26A","#DD7F88"))+
  xlab("MDS1") + ylab("MDS2") + # 轴标签
  labs(title = paste("bray stress =", stress, sep=" ")) + # 添加 stress 值到标题
  geom_text(aes(x = I(0.5), y = I(0.6), label = adonis_text), # 在图中显示 PERMANOVA 结果
            size =3, color ="black", inherit.aes =F) +
  theme(legend.position=c(0.83,0.18), legend.title=element_blank(),
        axis.title.x=element_text(size=10),#修改X轴标题文本
        axis.title.y=element_text(size=10,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=8),#修改x轴刻度标签文本
        axis.text.x=element_text(size=8),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
ggMarginal(
  p2,
  type=c('density'),
  margins='both',
  size=4.5,
  groupColour=F,
  groupFill=T
)

ggsave(p2, "../fig2-B.pdf") 