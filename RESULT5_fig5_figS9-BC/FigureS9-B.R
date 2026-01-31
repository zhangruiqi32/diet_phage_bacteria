##Script to Supplementary fig9-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("dplyr","tidyr","stringr","ggplot2","pheatmap","reshape2")
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

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(reshape2)
############################## Data loading ##############################
meta <- read.table("~/cooperation/202409zhaofengxiang/data/sample_meta")
AMG_contig_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_contig_contents.csv")
AMG_contents <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_contents.csv")
AMG_contig_length <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_contig_length.csv")
AMG_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_reads.csv")

############################## Data processing ##############################
AMG_contig_contents$scaffold2 <- substring(AMG_contig_contents$scaffold, 1, last = 30)
AMG_contig_contents<-merge(AMG_contig_contents,meta,by="Sample_id")
AMG_contig_contents$Phacts_result[grep("Temperate|Lytic",AMG_contig_contents$Probability)]<- AMG_contig_contents$Probability[grep("Temperate|Lytic",AMG_contig_contents$Probability)]

AMG_contents2 <- aggregate(AMG_contents$Present.AMG.KOs,by=list(AMG_contents$Metabolism),function(x){paste(unique(x),collapse = ",")})
AMG_contents3 <- data.frame(matrix(ncol = 2,nrow = 0))
for (i in 1:nrow(AMG_contents2)){
  kos <- unlist(strsplit(AMG_contents2[i,2],split = ","))
  AMG_contents4 <-data.frame(metabolism=AMG_contents2[i,1],AMG.KO=kos)
  AMG_contents3 <- rbind(AMG_contents3,AMG_contents4)
}
AMG_contents3 <- unique(AMG_contents3)
AMG_contents_metabolism <-merge(AMG_contents3,AMG_contig_contents[,colnames(AMG_contig_contents)%in%grep("Sample_id",colnames(AMG_contig_contents),value = T)[-1]==F],by="AMG.KO")

######### 差异KO
AMG_contig_contents2 <- AMG_contig_contents
AMG_contig_contents2$protein <- sapply(strsplit(AMG_contig_contents2$protein,split = " ",fixed = T),"[",1)
AMG_reads2 <- merge(AMG_reads[c("mapped_Contig","mapped_Sample_id","mapped_Normalized_Mapped_reads","mapped_GD","mapped_Diet")],
                    AMG_contig_contents2,by.x = "mapped_Contig",by.y = "protein")

ko_vector <-  c("K00789","K01582","K01953","K01953","K05396","K14157","K14157","K17103","K00558","K05396","K12960",
                "K17398")

ko_pathway <- AMG_contents[,c(2,3,5)]  %>%
  separate_rows(Present.AMG.KOs, sep = ",") %>%
  unique() %>%
  filter(Present.AMG.KOs %in% ko_vector) %>%
  filter(Metabolism %in% "Amino acid metabolism") 
  

ko_name <- AMG_contig_contents[, 4:5] %>%
  filter(AMG.KO %in% ko_vector) %>%
  unique() %>%
  mutate(KO_name = str_extract(AMG.KO.name,"(?<=;)[^\\[]+(?=\\[)"))


# 准备绘图数据
plot_data <- AMG_reads2 %>%
  # 汇总每个KO在每个样本中的表达量
  group_by(AMG.KO, mapped_Sample_id, mapped_Diet, mapped_GD) %>%
  summarise(expr = sum(mapped_Normalized_Mapped_reads), .groups = "drop") %>%
  # 添加代谢通路信息
  left_join(
    unique(AMG_contents_metabolism[, 1:2]),
    by = "AMG.KO"
  ) %>%
  # 计算每个时间点的均值和误差
  group_by(AMG.KO, mapped_Diet, mapped_GD, metabolism) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    sd_expr = sd(expr, na.rm = TRUE),
    se_expr = sd_expr / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  # 确保时间点是因子并正确排序
  mutate(
    mapped_GD = factor(mapped_GD, levels = paste0("D", sort(unique(as.numeric(gsub("D", "", mapped_GD)))))),
    AMG.KO = factor(AMG.KO)
  ) 
  
plot_data_log2 <- plot_data %>%
  mutate( mean_expr_log2= log2(mean_expr + 1),
          se_expr_log2= sd_expr / (sqrt(n_samples) * (mean_expr + 1) * log(2)))
plot_data_filtered <- plot_data_log2 %>%
  filter(AMG.KO %in% ko_vector) %>%
  filter(metabolism %in% "Amino acid metabolism") %>%
  # add KO name & pathway
  left_join(ko_name[,c(1,3)], by = "AMG.KO") %>%
  left_join(ko_pathway[,2:3], by = c("AMG.KO"="Present.AMG.KOs"))

plot_data_filtered$mapped_Diet <- factor(plot_data_filtered$mapped_Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC") )
# 创建分面折线图
ggplot(plot_data_filtered, aes(x = mapped_GD, y = mean_expr_log2, color = mapped_Diet, group = mapped_Diet)) +
  # 折线
  geom_line(size = 1) +
  # 点
  geom_point(size = 2) +
  # 误差线
  geom_errorbar(aes(ymin = mean_expr_log2 - se_expr_log2, ymax = mean_expr_log2 + se_expr_log2), 
                width = 0.1, size = 0.5) +
  # 分面：KO为行，组别为列
  facet_grid(rows = vars(AMG.KO), cols = vars(mapped_Diet), scales = "free_y") +
  # 美化主题
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  # 坐标轴标签
  labs(
    x = "Time (Day)",
    y = "log2(Normalized Expression)",
    color = "Diet",
    title = "KO Expression Time Course by Diet"
  ) +
  # 颜色设置
  scale_color_manual(values = c(CON="#8BBCD6","HFD" = "#EE7C70", "FUC" = "#ACD26A"))

############################## [FigS9-B] Processing and Painting #####################################
df <- plot_data_filtered
##### 相同Pathway放一起，且升高+降低都热图显示（D10）排序
# 1. 创建函数计算相对于D0的变化
calculate_d0_changes <- function(data) {
  # 添加一列表示每个KO在不同饮食组中D0的表达量
  data_with_d0 <- data %>%
    group_by(AMG.KO, mapped_Diet, KO_name, Pathway) %>%
    mutate(
      D0_expr_log2 = mean_expr_log2[mapped_GD == "D0"],
      D0_expr = mean_expr[mapped_GD == "D0"],
      change_vs_D0_log2 = mean_expr_log2 - D0_expr_log2,
      change_vs_D0_percent = round((mean_expr - D0_expr) / D0_expr * 100, 1)
    ) %>%
    ungroup()
  
  return(data_with_d0)
}

# 计算相对于D0的变化
df_changes <- calculate_d0_changes(df)

# 2. 筛选条件：
# a) HFD组中表达升高（相比D0）
# b) 在相同时间点，FUC组表达低于HFD组

# 首先，筛选HFD组表达升高的数据
hfd_increased <- df_changes %>%
  filter(
    mapped_Diet == "HFD",
    mapped_GD != "D0",  # 排除D0自身比较
    change_vs_D0_log2 > 0  # HFD组相比D0表达升高
  ) %>%
  select(
    AMG.KO, KO_name, Pathway, mapped_GD,
    HFD_expr_log2 = mean_expr_log2,
    HFD_change_vs_D0_log2 = change_vs_D0_log2,
    HFD_change_vs_D0_percent = change_vs_D0_percent
  )

# 然后，获取相同KO和时间点的FUC组数据
hfd_fuc_comparison <- hfd_increased %>%
  left_join(
    df_changes %>%
      filter(mapped_Diet == "FUC", mapped_GD != "D0") %>%
      select(
        AMG.KO, KO_name, Pathway, mapped_GD,
        FUC_expr_log2 = mean_expr_log2,
        FUC_change_vs_D0_log2 = change_vs_D0_log2,
        FUC_change_vs_D0_percent = change_vs_D0_percent
      ),
    by = c("AMG.KO", "KO_name", "Pathway", "mapped_GD")
  ) %>%
  # 筛选FUC表达低于HFD的条目
  filter(!is.na(FUC_expr_log2) & FUC_expr_log2 < HFD_expr_log2) %>%
  # 计算逆转幅度
  mutate(
    # FUC相比HFD的降低幅度（log2）
    reduction_vs_HFD_log2 = HFD_expr_log2 - FUC_expr_log2,
    # FUC相比HFD的降低百分比
    reduction_vs_HFD_percent = round((2^HFD_expr_log2 - 2^FUC_expr_log2) / 2^HFD_expr_log2 * 100, 1)
  ) %>%
  # 添加HFD相比FUC的倍数变化
  mutate(
    HFD_vs_FUC_ratio = round(2^(HFD_expr_log2 - FUC_expr_log2), 2)
  )

# 3. 创建HFD升高百分比热图数据
# 首先，选择出现在至少2个时间点的KO用于热图
ko_timepoints <- hfd_fuc_comparison %>%
  group_by(AMG.KO, KO_name, Pathway) %>%
  summarise(
    n_timepoints = n(),
    .groups = "drop"
  ) %>%
  filter(n_timepoints >= 2)  # 至少出现在2个时间点

# 构建热图数据（HFD升高百分比）
heatmap_hfd_data <- hfd_fuc_comparison %>%
  filter(AMG.KO %in% ko_timepoints$AMG.KO) %>%
  # 创建行标签
  mutate(
    row_label = paste(KO_name, " (", AMG.KO, ")", sep = "")
  ) %>%
  select(row_label, Pathway, mapped_GD, HFD_change_vs_D0_percent)

# 转换为宽格式
heatmap_hfd_wide <- heatmap_hfd_data %>%
  dcast(row_label + Pathway ~ mapped_GD, value.var = "HFD_change_vs_D0_percent")

# 处理缺失值（填充为0）
heatmap_hfd_wide[is.na(heatmap_hfd_wide)] <- 0

# 按Pathway分组，并在每个Pathway内按D10升高百分比降序排序
heatmap_hfd_sorted <- heatmap_hfd_wide %>%
  group_by(Pathway) %>%
  arrange(desc(D10), .by_group = TRUE) %>%
  ungroup()

# 提取排序后的行标签顺序
row_order_hfd <- heatmap_hfd_sorted$row_label

# 准备热图矩阵
heatmap_hfd_matrix <- heatmap_hfd_sorted %>%
  select(-row_label, -Pathway) %>%
  as.matrix()

rownames(heatmap_hfd_matrix) <- heatmap_hfd_sorted$row_label

# 准备行注释（Pathway信息）
row_annotation_hfd <- heatmap_hfd_sorted %>%
  select(Pathway) %>%
  as.data.frame()
rownames(row_annotation_hfd) <- heatmap_hfd_sorted$row_label

# 确保列顺序一致（D10, D14, D18）
all_timepoints <- c("D10", "D14", "D18")
missing_cols_hfd <- setdiff(all_timepoints, colnames(heatmap_hfd_matrix))

if(length(missing_cols_hfd) > 0) {
  for(col in missing_cols_hfd) {
    heatmap_hfd_matrix <- cbind(heatmap_hfd_matrix, rep(0, nrow(heatmap_hfd_matrix)))
    colnames(heatmap_hfd_matrix)[ncol(heatmap_hfd_matrix)] <- col
  }
}

# 重新排列列顺序
heatmap_hfd_matrix <- heatmap_hfd_matrix[, all_timepoints, drop = FALSE]

# 4. 绘制HFD升高百分比热图（按Pathway聚类）
if(nrow(heatmap_hfd_matrix) > 0) {
  # 创建颜色映射（使用蓝色渐变色系表示升高百分比）
  color_palette_hfd <- colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDBB84", "#EF6548"))(50)
  
  # 确定颜色范围，使0值显示为白色
  max_value <- max(heatmap_hfd_matrix, na.rm = TRUE)
  
  # 绘制热图
  pheatmap(
    heatmap_hfd_matrix,
    cluster_rows = FALSE,  # 不进行行聚类，使用我们的排序
    cluster_cols = FALSE,  # 不进行列聚类
    color = color_palette_hfd,
    main = "HFD相比D0的表达升高百分比(%) - 按Pathway分组",
    display_numbers = TRUE,
    number_format = "%.0f",
    fontsize_row = 9,
    fontsize_col = 10,
    annotation_row = row_annotation_hfd,
    gaps_row = cumsum(table(row_annotation_hfd$Pathway))[-length(table(row_annotation_hfd$Pathway))],
    show_rownames = TRUE,
    show_colnames = TRUE,
    breaks = seq(0, max_value, length.out = 51)  # 确保颜色从0开始
  )
}

# 5. 创建组合热图：左侧HFD升高百分比，右侧FUC降低百分比
# 首先获取FUC降低百分比的数据（与之前相同，但使用相同的行顺序）
heatmap_fuc_data <- hfd_fuc_comparison %>%
  filter(AMG.KO %in% ko_timepoints$AMG.KO) %>%
  mutate(
    row_label = paste(KO_name, " (", AMG.KO, ")", sep = "")
  ) %>%
  select(row_label, mapped_GD, reduction_vs_HFD_percent = reduction_vs_HFD_percent)

# 转换为宽格式
heatmap_fuc_wide <- heatmap_fuc_data %>%
  dcast(row_label ~ mapped_GD, value.var = "reduction_vs_HFD_percent")

# 处理缺失值（填充为0）
heatmap_fuc_wide[is.na(heatmap_fuc_wide)] <- 0

# 确保行顺序与HFD热图相同
heatmap_fuc_wide <- heatmap_fuc_wide %>%
  mutate(row_label = factor(row_label, levels = row_order_hfd)) %>%
  arrange(row_label) %>%
  mutate(row_label = as.character(row_label))

# 准备FUC热图矩阵
heatmap_fuc_matrix <- heatmap_fuc_wide %>%
  select(-row_label) %>%
  as.matrix()

rownames(heatmap_fuc_matrix) <- heatmap_fuc_wide$row_label

# 确保列顺序一致
missing_cols_fuc <- setdiff(all_timepoints, colnames(heatmap_fuc_matrix))

if(length(missing_cols_fuc) > 0) {
  for(col in missing_cols_fuc) {
    heatmap_fuc_matrix <- cbind(heatmap_fuc_matrix, rep(0, nrow(heatmap_fuc_matrix)))
    colnames(heatmap_fuc_matrix)[ncol(heatmap_fuc_matrix)] <- col
  }
}

heatmap_fuc_matrix <- heatmap_fuc_matrix[, all_timepoints, drop = FALSE]

# 6. 创建并排热图展示HFD升高和FUC降低的对比
library(gridExtra)
library(grid)

# 创建HFD热图
p1 <- pheatmap(
  heatmap_hfd_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#FFFFFF", "#FC8D59", "#EF6548", "#D7301F"))(50),
  border=NA, #热图每个小单元格边框颜色
  main = "HFD升高百分比(%)",
  display_numbers = TRUE,
  number_format = "%.0f",
  fontsize_row = 8,
  fontsize_col = 9,
  annotation_row = row_annotation_hfd,
  show_rownames = F,
  show_colnames = TRUE,
  silent = TRUE  # 不直接显示，用于组合
)

# 创建FUC热图
p2 <- pheatmap(
  heatmap_fuc_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#FAFDFF", "#C6DBEF", "#6BAED6", "#2171B5"))(50),
  border=NA, #热图每个小单元格边框颜色
  main = "FUC降低百分比(%)",
  display_numbers = TRUE,
  number_format = "%.0f",
  fontsize_row = 8,
  fontsize_col = 9,
  show_rownames = TRUE,
  show_colnames = TRUE,
  silent = TRUE  # 不直接显示，用于组合
)

# 使用grid.arrange组合两个热图
grid.arrange(
  p1[[4]], p2[[4]],
  ncol = 2,
  top = textGrob("HFD引起的表达变化及FUC的逆转效应", gp = gpar(fontsize = 16, fontface = "bold")),
  bottom = textGrob("注：左侧为HFD相比D0的升高百分比，右侧为FUC相比HFD的降低百分比", 
                    gp = gpar(fontsize = 10))
)

# 7. 按通路展示关键发现（提供每个时间点对应的百分比）
cat("\n=== 按通路展示关键发现 ===\n\n")

# 获取所有有逆转效应的通路
pathways <- unique(hfd_fuc_comparison$Pathway)

for(pathway in pathways) {
  pathway_data <- hfd_fuc_comparison %>%
    filter(Pathway == pathway) %>%
    arrange(AMG.KO, mapped_GD)
  
  if(nrow(pathway_data) > 0) {
    cat(paste0("=== ", pathway, "通路 ===\n"))
    
    # 获取该通路下的所有基因
    genes_in_pathway <- unique(pathway_data$AMG.KO)
    
    for(gene_ko in genes_in_pathway) {
      gene_data <- pathway_data %>%
        filter(AMG.KO == gene_ko) %>%
        arrange(mapped_GD)
      
      if(nrow(gene_data) > 0) {
        # 获取基因名称
        gene_name <- unique(gene_data$KO_name)[1]
        
        cat(paste0("\n基因: ", gene_name, " (", gene_ko, ")\n"))
        
        # 为每个时间点显示百分比
        for(i in 1:nrow(gene_data)) {
          time_point <- gene_data$mapped_GD[i]
          hfd_increase <- gene_data$HFD_change_vs_D0_percent[i]
          fuc_reduction <- gene_data$reduction_vs_HFD_percent[i]
          
          cat(paste0("  ", time_point, "时间点: ",
                     "HFD相比D0升高", hfd_increase, "%",
                     "，FUC相比HFD降低", fuc_reduction, "%\n"))
        }
      }
    }
    cat("\n" , paste(rep("-", 60), collapse = ""), "\n\n")
  }
}

