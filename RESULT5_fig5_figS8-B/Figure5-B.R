##Script to Figure 5-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","tidyverse","ggtreeExtra","ggtree","ggstar","ggnewscale","ape","RColorBrewer","ggplot2")
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
tree<- read.tree("../data/fig5B_BV_AMG_trimal_fasta.treefile")
dat1 <- read.table("../data/fig5B_BV_anno_2",sep = "\t",header = T)

########################数据导入+处理##################################
tree<- drop.tip(tree,c("NODE_284108_length_84_cov_4.586207_1"))
tree<-drop.tip(tree,c("NODE_195514_length_101_cov_4.217391_1"))
tree<-drop.tip(tree,c("NODE_288280_length_84_cov_3.206897_1"))
tree <- root(tree,1)
# drop.tip(tree,tree$tip.label[tree$edge.length])
data = fortify(tree)

a <- unique(dat1[c("Phacts_result","Contig")])
df_id <- data.frame(id=tree$tip.label)
df_id$Contig2 <- substring(df_id$id,1,30)
dat1 <- merge(dat1,df_id,by="Contig2",all.x = T)
dat1$Con_change <- as.numeric(dat1$Con_D10)-as.numeric(dat1$Con_D0)
dat1$HFD_change <- as.numeric(dat1$HFD_D10)-as.numeric(dat1$HFD_D0)
dat1$HFD.Fuco_change <- as.numeric(dat1$HFD.Fuco_D10)-as.numeric(dat1$HFD.Fuco_D0)
dat1[is.na(dat1)] <- 0
dat1$Con_change[dat1$Con_change>0] <- log10(abs(dat1$Con_change[dat1$Con_change>0])+1)
dat1$Con_change[dat1$Con_change<0] <- -1*log10(abs(dat1$Con_change[dat1$Con_change<0])+1)
dat1$HFD_change[dat1$HFD_change>0] <- log10(abs(dat1$HFD_change[dat1$HFD_change>0])+1)
dat1$HFD_change[dat1$HFD_change<0] <- -1*log10(abs(dat1$HFD_change[dat1$HFD_change<0])+1)
dat1$HFD.Fuco_change[dat1$HFD.Fuco_change>0] <- log10(abs(dat1$HFD.Fuco_change[dat1$HFD.Fuco_change>0])+1)
dat1$HFD.Fuco_change[dat1$HFD.Fuco_change<0] <- -1*log10(abs(dat1$HFD.Fuco_change[dat1$HFD.Fuco_change<0])+1)
dat1 <- dat1[order(dat1$Contig),]
dat1$contig_type <- gsub("k__","",dat1$contig_type)
dat1$contig_type[dat1$contig_type=="AMG"] <- "Viruses"
color_point <- c(brewer.pal(9,"Set1")[2],"gray",brewer.pal(9,"Set2")[2])
color_point_name <- unique(dat1$contig_type)
color_point_name <- color_point_name[order(color_point_name)]
names(color_point) <-   color_point_name
# dat1 <- dat1[!duplicated(dat1$Contig),]
dat2 <-dat1
dat2$pheatmap <- 1
dat2$pheatmap_content <- dat1$metabolism
dat2$pheatmap_color <- dat1$metabolism_color
dat2$pheatmap_color <- as.character(factor(dat2$pheatmap_color,levels = unique(dat2$pheatmap_color),labels= brewer.pal(9,"Set3")  )   )
# dat2$pheatmap_content <- as.numeric(factor(dat2$pheatmap_color,levels = unique(dat2$pheatmap_color),labels= brewer.pal(9,"Set1")  )   )
dat2$pheatmap_color[dat2$pheatmap_color=="gray"] <- "white"
dat4 <-dat2
dat4 <-unique(dat4[c("id","pheatmap","pheatmap_content","pheatmap_color")])
dat4 <-dat4[order(dat4$pheatmap_content),]
dat4 <-dat4[!duplicated(dat4$id,dat4$pheatmap_content),]
color2 <- unique(dat4$pheatmap_color)
names(color2) <- unique(dat4$pheatmap_content)
dat3 <-dat1
dat3$pheatmap <- 2
dat3$pheatmap_content <- dat1$Phacts_result
dat3 <-unique(dat3[c("Contig","id","pheatmap","pheatmap_content")])
dat3 <-dat3[order(dat3$pheatmap_content),]
dat3 <- dat3[dat3$pheatmap_content!="None",]
dat3 <- aggregate(1:nrow(dat3),by=list(id=dat3$id,pheatmap_content=dat3$pheatmap_content),function(x) length(x) )

dat3 <-spread(dat3,key="pheatmap_content",value="x")
dat3[is.na(dat3)] <- 0
dat3$Ratio_of_Temperate_and_Lytic  <-log10( (as.numeric(dat3$Temperate)+1)/(as.numeric(dat3$Lytic)+1))
dat3$pheatmap <- 2
dat3$pheatmap_content <- dat3$Ratio_of_Temperate_and_Lytic
dat <-dat1
dat$pheatmap <-3
dat$pheatmap_content <- dat1$Con_change
dat$pheatmap_color <- dat1$Con_change
dat2 <-dat
dat <-dat1
dat$pheatmap <-4
dat$pheatmap_content <- dat1$HFD_change
dat$pheatmap_color <- dat1$HFD_change
dat2 <-rbind(dat2,dat)
dat <-dat1
dat$pheatmap <-5
dat$pheatmap_content <- dat1$HFD.Fuco_change
dat$pheatmap_color <- dat1$HFD.Fuco_change
dat2 <-rbind(dat2,dat)


############################## Data analysis #####################################
#溶源态下辅因子和维生素代谢相关的AMG占比
result <- dat1 %>%
  filter(Phacts_result == "Temperate") %>%  # 筛选Temperate行
  summarise(
    total_temperate = n(),  # 总Temperate计数
    target_count = sum(metabolism == "Metabolism_of_cofactors_and_vitamins"),  # 目标代谢计数
    percentage = target_count / total_temperate * 100  # 计算百分比
  )

# 计算每种metabolism的占比
metabolism_stats <- dat1 %>%
  filter(Phacts_result == "Temperate") %>%  # 筛选溶原性病毒
  count(metabolism, name = "count") %>%    # 按代谢类型计数
  mutate(
    total = sum(count),                    # 计算总数
    percentage = count / total * 100,      # 计算百分比
    percentage_label = sprintf("%.2f%%", percentage)  # 格式化百分比
  ) %>%
  arrange(desc(count))                     # 按计数降序排列

#氨基酸代谢中烈性/溶源占比
aa_phacts_stats <- dat1 %>%
  filter(metabolism == "Amino_acid_metabolism") %>%  # 筛选氨基酸代谢的行
  count(Phacts_result, name = "count") %>%           # 按Phacts_result计数
  mutate(
    total = sum(count),                             # 计算总数
    percentage = count / total * 100,                # 计算百分比
    percentage_label = sprintf("%.2f%%", percentage) # 格式化百分比
  ) %>%
  arrange(desc(count))                              # 按计数降序排列

#氨基酸代谢中HFD_change不同值的占比
HFD_change_stats <- dat1 %>%
  filter(metabolism == "Amino_acid_metabolism") %>%  # 筛选氨基酸代谢的行
  mutate(
    change_group = case_when(
      HFD_change > 0 ~ ">0 (丰度增加)",
      HFD_change < 0 ~ "<0 (丰度减少)",
      HFD_change == 0 ~ "=0 (丰度不变)",
      TRUE ~ "NA"  # 处理缺失值
    )
  ) %>%
  count(change_group, name = "count") %>%           # 按变化组计数
  filter(change_group != "NA") %>%                  # 排除缺失值
  mutate(
    total = sum(count),                             # 计算总数
    percentage = count / total * 100,               # 计算百分比
    percentage_label = sprintf("%.2f%%", percentage) # 格式化百分比
  ) %>%
  arrange(match(change_group, c("<0 (减少)", "=0 (不变)", ">0 (增加)")))  # 按逻辑顺序排列

############################## [Fig5-B] Painting #####################################
phacts_color <- c(brewer.pal(11,"Spectral")[2],brewer.pal(11,"Spectral")[10] )
phacts_color <- c(brewer.pal(11,"Set2")[5],brewer.pal(11,"Set2")[6] )
phacts_color <- c(brewer.pal(9,"Set1")[4],brewer.pal(9,"Set1")[5] )
Abudance_color <- c(brewer.pal(11,"RdYlBu")[2],brewer.pal(11,"RdYlBu")[10] )

p<- ggtree(tree, layout="fan", open.angle=90, size=1)+
  geom_fruit(data=dat1,geom=geom_point,
             mapping=aes(y=id, color=contig_type),position="identity", size=1.5)+
              scale_color_manual(values = color_point)+
              # scale_fill_manual(values = color_point)+
              theme(legend.key = element_rect(fill = NA,color = 'transparent'))+
              guides(color=guide_legend(#title = "Gene source",
                                        override.aes = list(size=2)   )     )+
  new_scale_fill() +
  geom_fruit(
    data=dat4,
    geom=geom_bar,
    mapping=aes(y=id, x=pheatmap, fill=pheatmap_content),
    orientation="y",
    stat="identity",
    pwidth=0.1,
    offset=0,
    width=1.2,alpha=0.8
  ) +
  scale_fill_manual(values = color2)+
  guides(fill=guide_legend(title = "Metabolic type"))+
  new_scale_fill() +
  geom_fruit(
    data=dat3,
    geom=geom_bar,
    mapping=aes(y=id, x=pheatmap, fill=pheatmap_content),
    orientation="y",
    stat="identity",
    pwidth=0.1,
    offset=0,
    width=1.2
  )+
  # scale_fill_manual(values = color)+
  scale_fill_gradient2(low = phacts_color[2], mid = "white", high = phacts_color[1])+
  guides(fill=guide_legend(title = "Metabolic type"))+
  new_scale_fill() +
  geom_fruit(
      data=dat2,
      geom=geom_tile,
      mapping=aes(y=id, x=pheatmap, fill=pheatmap_content),
      offset=0,
      pwidth=0.1
    )+
  guides(fill=guide_legend(title = "Gene abundance changes"))+
  scale_fill_gradient2(low = Abudance_color[2], mid = "white", high = Abudance_color[1])
p  
ggsave("../fig5-B.pdf")









