##Script to Figure 2-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("dplyr","pheatmap","RColorBrewer","ggplot2")
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
data <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2C_DifferentBacts_viruse.csv", header=TRUE)
ann_col <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2_sample.csv", header=TRUE)  


############################## [Fig2C] Processing and Painting #####################################
data$viruses <- gsub(".*s__(\\S+).*", "s__\\1", data$viruses)
data$host <- sapply(strsplit(data$viruses,split = "_",fixed = T),"[",3)
data$virus <- sapply(strsplit(data$viruses,split = "_virus_",fixed = T),"[",2)
data$phage <- sapply(strsplit(data$viruses,split = "_phage_",fixed = T),"[",2)
data$type <- ifelse(is.na(data$virus),"phage","virus")
data$name <- ifelse(is.na(data$virus),data$phage,data$virus)

order_data <- arrange(data, host,type)
paintdata <- order_data[,2:58]
rownames(paintdata) <- order_data[,63]
paintdata <- log10(paintdata+0.000000001) #丰度太低了 log一下

rownames(ann_col) <- colnames(paintdata) #必须有这行
ann_col$time <- factor(ann_col$time)
ann_col$treatment <- factor(ann_col$treatment)

ann_row <- order_data[,c(59,62)]
rownames(ann_row) <- rownames(paintdata)
ann_row$host <- factor(ann_row$host)
ann_row$type <- factor(ann_row$type)

ann_color <- list(time = c(D0="#F8F3EB",D10="#EAE1EF",
                           D14="#C2DCE9",D18="#C2DACF"),
                  treatment = c(CON="#8cbcd4",
                                HFD="#dc7c8c",
                                FUC="#abd16d"),
                  type= c(phage="#C47070",
                          virus="#BEBAB9"),
                  host=c(Acinetobacter="#F08961",
                         Bacteroides="#DA87B6",
                         Bacillus="#5E4FA2",
                         Flavobacterium="#E6F598",
                         Fusobacterium= "#ABDDA4",
                         Stenotrophomonas="#66C2A5",
                         Serratia="#237B7A",
                         Clostridioides="#A6CEE3",
                         Campylobacter="#3288BD",
                         Gordonia="#9FAAD1"
                  ))
mycol <- colorRampPalette(c("#DCDDDD","white", "orange"))(10) #更改色差

pheatmap(paintdata,scale="row",#均一化方向
         fontsize = 5,
         fontsize_row = 5,
         fontsize_col = 4,
         cellwidth = 4,
         cellheight = 6,
         color = mycol,
         cluster_rows = F, cluster_cols = F, #对行/列聚类
         #cutree_rows = NA, cutree_cols = T, #根据行/列聚类数量分隔热图行/列
         gaps_col=c(21,38), #未进行列聚类使用，在第21列和第38列添加gap
         treeheight_row = 20, treeheight_col=30, #对聚类树高度调整
         border=NA, #热图每个小单元格边框颜色
         annotation_col=ann_col[,2:3], annotation_row =ann_row,
         annotation_color=ann_color,
         annotation_legend = T,
         annotation_names_col = T, annotation_names_row = T)

ggsave("../fig2-C.pdf") 