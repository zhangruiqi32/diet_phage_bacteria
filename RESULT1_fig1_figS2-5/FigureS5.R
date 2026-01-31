##Script to Supplementary fig3.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","ggtree","ape","treeio","tidytree","dplyr","ggtreeExtra","ggnewscale","RColorBrewer","ggplot2")
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
tree <- read.tree("~/cooperation/202409zhaofengxiang/data/figS5_phylogenetic_tree.nwk")
metadata <- read.csv("~/cooperation/202409zhaofengxiang/data/figS5_metadata.csv")


############################## Data processing ##############################
# 获取叶节点元数据
num_tips <- length(tree$tip.label)  # 叶节点总数
tip_indices <- 1:num_tips          # 叶节点编号范围（1到24）

# 筛选需要删除的分支
threshold <- 1.5  # 设置长度阈值

# 创建筛选矩阵
is_tip_edge <- tree$edge[,2] %in% tip_indices  # 标记叶节点分支（TRUE/FALSE）
long_edges <- which(tree$edge.length > threshold & is_tip_edge)  # 双重条件筛选

# 获取待删除的叶节点标签
tips_to_remove <- tree$edge[long_edges, 2] %>% 
  unique() %>% 
  sapply(function(x) tree$tip.label[x])  # 将节点编号转换为实际标签

# 执行删除操作
if(length(tips_to_remove) > 0){
  clean_tree <- drop.tip(tree, tip = tips_to_remove)
  message("成功删除 ", length(tips_to_remove), " 个过长分支:")
  print(tips_to_remove)
} else {
  clean_tree <- tree
  message("未发现长度超过", threshold, "的分支")
}
ggtree(clean_tree, layout = "circular")

# 注释表
metadata <- metadata[is.na(metadata$gene)==F,]
valid_ids <- metadata$gene
metadata$ProteinID<- metadata$gene
tree <- keep.tip(tree, valid_ids)


############################## [FigS4] Painting #####################################
# 基础树图
# p <- ggtree(tree, layout = "fan", size= 0.05,branch.length = "none") 
p<- ggtree(clean_tree, layout = "circular", size= 0.25)
p<- ggtree(clean_tree, layout = "fan", size= 0.25)

# 病毒
metadata$NR_annotations_taxon_classfied <-sapply(strsplit(metadata$species_annotation,split = " ",fixed = T),"[",1)
a <- table(metadata$NR_annotations_taxon_classfied )
a <- a[a>=9]
metadata$NR_annotations_taxon_classfied[metadata$NR_annotations_taxon_classfied%in%names(a)==F] <- "Others"
library(ggstar)
p1 <- p + geom_fruit(
  data=metadata,
  geom=geom_star,
  mapping=aes(y = gene,fill=NR_annotations_taxon_classfied),
  starshape= 15,
  size=1,color=NA,
  position="identity",
  offset = 0.02 # 增加offset避免重叠
)+ 
  scale_fill_manual(name = "Virus Family", 
                    values = list("Myoviridae"="#c4a5cc","Podoviridae"="#8A9CC4",
                                  "Siphoviridae"="#abd3a3","Microviridae"="#FDC086","Others"="gray"),na.value = "gray")
# 重置fill标度以兼容第二个连续型变量
p1 <- p1 + ggnewscale::new_scale_fill()

# 宿主
metadata$family2 <-sapply(strsplit(metadata$host_annotation,split = " ",fixed = T),"[",1)
a <- table(metadata$family2 )
a <- a[a>=10]
metadata$family2[metadata$family2%in%names(a)==F] <- "Others"
metadata$family2[metadata$family2%in%c("unclassified_Bacteria_phylum")] <- "Others"
p2 <- p1 +
  geom_fruit(
    data = metadata,
    geom = geom_col, # 使用geom_col代替geom_bar以直接映射数值
    mapping = aes(
      y = gene,
      x = 0.02,         # 固定 x 值以匹配 width = 0.3
      fill = family2 ),
    width = 1,pwidth = 0.04,
    na.rm = TRUE,
    offset = 0.01 # 增加offset避免重叠
  )+
  scale_fill_manual(name = "Host annotation", 
                    values = list("Actinomycetota"="#7CA5C2","Bacillota"="#DDEAC4",
                                  "Bacteroidota"="#85C3B9","Pseudomonadota"="#F8D7E7","Others"="#EFEFEF"),na.value = "gray")
# 重置fill标度以兼容第二个连续型变量
p2 <- p2 + ggnewscale::new_scale_fill()

#生活方式
metadata$Lytic[metadata$Lytic<0.2]=0.2
metadata$Lytic[metadata$Lytic>0.8]=0.8
# 添加第二个fruit图层（连续型变量）
p3 <- p2 +
  geom_fruit(
    data = metadata,
    geom = geom_col, # 使用geom_col代替geom_bar以直接映射数值
    mapping = aes(
      y = gene,
      x = 0.02,         # 固定 x 值以匹配 width = 0.3
      fill = log10(as.numeric(Lytic) / (1 - as.numeric(Lytic))) ),
    width = 1,pwidth = 0.04,
    na.rm = TRUE,
    offset = 0.01 # 增加offset避免重叠
  ) +
  scale_fill_gradient2(
    name = "Lytic Score",          # 图例标题
    low = "#BA90C0",                  # 负值颜色（最小）
    mid = "#F7F8F8",                 # 中间值颜色（0点）
    high = "#EF7E1F",                  # 正值颜色（最大）
    midpoint = 0  ,              # 明确指定中间点为0
    na.value = "#F7F8F8"
  )
# 重置fill标度以兼容第二个连续型变量
p3 <- p3 + ggnewscale::new_scale_fill()


#饮食
metadata$Diet_Type[metadata$Diet_Type%in%c("Fermented","Gluten","Sugar","Protein")] <- "Others"
#control:灰 #D3D3D3，Fat:红 #F3A59A，Fiber:绿 #80D0C3，Special population diet:蓝 #A6DDEA，Other:紫 #9EAAC4（Fermented，Gluten，Protein，Sugar）
p4 <- p3 +
  geom_fruit(
    data = metadata,
    geom = geom_col,
    mapping = aes(y = gene,x =0.02, fill = Diet_Type),
    width = 1,pwidth = 0.02,
    offset = 0.01
    # orientation = "y"
  ) +
  scale_fill_manual(values = list("Control"= NA,
                                  "Fiber"=NA,
                                  "Special population diet"=NA,
                                  "Fat"="#FE9C9D",
                                  "Others"=NA),
                    na.value = NA) # 为离散变量命名标度
p4 <- p4 + ggnewscale::new_scale_fill()

p5 <- p4 +
  geom_fruit(
    data = metadata,
    geom = geom_col,
    mapping = aes(y = gene,x =0.02, fill = Diet_Type),
    width = 1,pwidth = 0.02,
    offset = 0
    # orientation = "y"
  ) +
  scale_fill_manual(values = list("Control"= NA,
                                  "Fiber"="#98C897",
                                  "Special population diet"=NA,
                                  "Fat"=NA,
                                  "Others"=NA),
                    na.value = NA) # 为离散变量命名标度
p5 <- p5 + ggnewscale::new_scale_fill()

p6 <- p5 +
  geom_fruit(
    data = metadata,
    geom = geom_col,
    mapping = aes(y = gene,x =0.02, fill = Diet_Type),
    width = 1,pwidth = 0.02,
    offset = 0
    # orientation = "y"
  ) +
  scale_fill_manual(values = list("Control"= NA,
                                  "Fiber"=NA,
                                  "Special population diet"="#A6DDEA",
                                  "Fat"=NA,
                                  "Others"=NA),
                    na.value = NA) # 为离散变量命名标度
p6 <- p6 + ggnewscale::new_scale_fill()

p7 <- p6 +
  geom_fruit(
    data = metadata,
    geom = geom_col,
    mapping = aes(y = gene,x =0.02, fill = Diet_Type),
    width = 1,pwidth = 0.02,
    offset = 0
    # orientation = "y"
  ) +
  scale_fill_manual(values = list("Control"= NA,
                                  "Fiber"=NA,
                                  "Special population diet"=NA,
                                  "Fat"=NA,
                                  "Others"="#9EAAC4"),
                    na.value = NA) # 为离散变量命名标度

p7

ggsave(p7, "../figS5.pdf") 


