##Script to Supplementary fig10-A.

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

library(tidyr)
############################## Data loading ##############################
data_g_B <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_B")
data_g_B <- data_g_B[grep("s__",rownames(data_g_B),fixed = T),]
data_g_B<- apply(data_g_B, 2, function(col) { col / sum(col, na.rm = TRUE)}) %>% as.data.frame()

data_g_V <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_V")
data_g_V <- data_g_V[grep("s__",rownames(data_g_V),fixed = T),]
data_g_V<- apply(data_g_V, 2, function(col) { col / sum(col, na.rm = TRUE)}) %>% as.data.frame()

prj2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv", header = TRUE)

############################## Data processing ##############################
lefse_results2 <- read.csv("~/cooperation/202409zhaofengxiang/data/figS9B_merged_result.txt",
                           sep="\t", header = F)
colnames(lefse_results2) <- c("Taxonomy","LDA_score","Group","LDA_score_adj","P_value","Filename")
lefse_results2$Filename <- sapply(strsplit(lefse_results2$Filename,split = "/out_",fixed = T),"[",2)
lefse_results2$Filename <- gsub(".res","",lefse_results2$Filename)
a <- unique(lefse_results2$Filename[grep("High_fat",lefse_results2$Group)])
lefse_results2 <- lefse_results2[lefse_results2$Filename%in%a,]
lefse_results2 <- lefse_results2[is.na(lefse_results2$LDA_score_adj)==F,]
lefse_results2 <- pivot_wider(lefse_results2[c("Taxonomy","Group","Filename")],names_from = "Filename",values_from = "Group",values_fill = "")
lefse_results2<-lefse_results2[order(rowSums(lefse_results2=="High_fat_diet"),decreasing = T),]
lefse_results2<-lefse_results2[rowSums(lefse_results2=="High_fat_diet")>=3,]

prjs_duotang <- c("PRJNA615253","PRJNA731974","mgp6153")

result <- data.frame(matrix(nrow = 0, ncol = ncol(data_g_V)))
colnames(result) <- colnames(data_g_V)
for(tax in lefse_results2$Taxonomy) {
  tax_name <- tail(strsplit(tax, "\\.")[[1]], 1)
  level_prefix <- substr(tax_name, 1, 3)
  if(level_prefix == "s__") {
    sp_name <- gsub("s__", "", tax_name) %>% gsub("\\.", "_", .)
    if(sp_name %in% rownames(data_g_V)) {
      result[tax_name,] <- data_g_V[sp_name, ]
    }
  } else {
    pattern <- paste0(level_prefix, gsub(paste0("^", level_prefix), "", tax_name))
    matching_rows <- grep(pattern, rownames(data_g_V))
    if(length(matching_rows) > 0) {
      result[tax_name,] <- colSums(data_g_V[matching_rows, ], na.rm = TRUE)
    }
  }
}

prj2 <- prj2[prj2$BioProject %in% prjs_duotang, ]
prj2 <- prj2[prj2$Diet1=="High_fat_diet",]
prj2 <- prj2[prj2$Nutrition1%in%c("probiotics")==F,]
annotation <- data.frame(row.names = colnames(result), Diet1 = prj2$Nutrition1[match(colnames(result), prj2$Run)],
                         BioProject=prj2$BioProject[match(colnames(result), prj2$Run)] )
result2 <- result[, order(annotation$BioProject,annotation$Diet1)]
result2 <- result2[,colnames(result2)%in%prj2$Run]
result2 <- result2[rowSums(result2)>0,colSums(result2)>0]
result2 <- unique(result2)
result2[result2<0] <- 0
annotation <- annotation[order(annotation$BioProject,annotation$Diet1), , drop = FALSE]

################################ [FigS10-A] Painting ################################
library(pheatmap)
ann_col <- drop_na(annotation)
ann_col$BioProject <- factor(ann_col$BioProject,
                             levels = c("PRJNA615253","PRJNA731974","mgp6153"))
ann_col$Diet1 <- factor(ann_col$Diet1,
                                 levels = c("Control","APS", "Inulin", "FOS"))
ann_col <- arrange(ann_col, BioProject,Diet1)
new_order <- rownames(ann_col)
order_data <- result2[,c(new_order)]

ann_color <- list(Diet1 = c(Control ="#F8F3EB", APS = "#2F5B60",
                           Inulin = "#AAD3C7", FOS = "#9CB4AD"),
                  BioProject = c(PRJNA615253="#2F5B60",
                                 PRJNA731974="#AAD3C7",
                                 mgp6153="#9CB4AD"))
mycol <- colorRampPalette(c("#DCDDDD","#fcf8f8", "#C47070"))(10) #更改色差

pheatmap(log10(order_data + 1),scale="row",#均一化方向
         fontsize = 6,
         fontsize_row = 6,
         #fontsize_col = 4,
         cellwidth = 8,
         cellheight = 10,
         color = mycol,
         cluster_rows = T, cluster_cols = F, #对行/列聚类
         #cutree_rows = NA, cutree_cols = T, #根据行/列聚类数量分隔热图行/列
         gaps_col=c(6,14), #未进行列聚类使用，在第6列和第14列添加gap
         treeheight_row = 6, treeheight_col=20, #对聚类树高度调整
         border=NA, #热图每个小单元格边框颜色
         annotation_col=ann_col, 
         #annotation_row =ann_row,
         annotation_color=ann_color,
         annotation_legend = T,
         annotation_names_col = T, annotation_names_row = T)

Diet1 = c(Control ="#F8F3EB", APS = "#829da0",
          Inulin = "#cce5dd", FOS = "#c4d2ce")


