##Script to Supplementary fig7.

############################## Packages loading ############################## 
cran_packages=c("vegan","ggplot2","tidyverse","ggExtra")
bioconductor_packages=c()

for(package in cran_packages){
  library(package, character.only=T,quietly = T, warn.conflicts = F)
}

library(ggExtra)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyverse)
############################## Data loading ##############################
group <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2_sample.csv")

#VIRUSE
#otu_raw <- read.csv("/wang/students/zhaofengxiang/Fucoidan_5.4/foradonis/Viruses",sep = "\t")
#otu_raw <- data.frame(OTU=rownames(otu_raw),otu_raw)
#otu_raw <- otu_raw[-grep("s__uncultured_crAssphage",rownames(otu_raw)),]

otu_raw <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2B_Bacts.csv")
otu_raw <- otu_raw[grep("s__", otu_raw$OTU), ]

otu_raw$OTU <- gsub(".*s__(\\S+).*", "s__\\1", otu_raw$OTU)
rownames(otu_raw) <- otu_raw$OTU
otu_raw <- otu_raw[,-1]
otu <- t(otu_raw)
############################## [FigS7-A] Processing and Painting #####################################
# 定义颜色
color=c("#8BBCD6","#ACD26A","#DD7F88")

### 使用所有数据计算MDS坐标
# 使用所有样本计算Bray-Curtis距离
all_otu_distance <- vegdist(otu, method = 'bray')
all_bray <- as.matrix(all_otu_distance)

# 计算所有样本的NMDS坐标
all_nmds <- metaMDS(all_bray, k=2)
all_stress <- round(all_nmds$stress, 4)

# 提取所有样本的NMDS坐标
all_plot_data <- all_nmds$points %>% 
  as.data.frame() %>%
  rownames_to_column(var="sample") %>%
  left_join(., group, by="sample")

# 重命名MDS坐标列
colnames(all_plot_data)[2:3] <- c("MDS1", "MDS2")
cat("Global NMDS stress (所有样本一起计算):", all_stress, "\n")

### 为每个时间点分别计算PERMANOVA
# 获取唯一的时间点
unique_times <- unique(group$time)

# 为每个时间点计算PERMANOVA结果
permanova_results <- data.frame()

for(time_point in unique_times) {
  # 筛选当前时间点的样本
  time_samples <- group$sample[group$time == time_point]
  time_otu <- otu[rownames(otu) %in% time_samples, ]
  time_group <- group[group$sample %in% time_samples, ]
  
  # 确保样本顺序一致
  time_otu <- time_otu[match(time_group$sample, rownames(time_otu)), ]
  
  # 计算当前时间点的距离矩阵
  otu_distance <- vegdist(time_otu, method = 'bray')
  bray <- as.matrix(otu_distance)
  
  # 为当前时间点进行PERMANOVA分析
  time_df <- data.frame(sample = rownames(time_otu))
  time_df <- merge(time_df, time_group, by = "sample")
  otu_adonis <- adonis2(bray ~ treatment, data = time_df, permutations = 999)
  
  # 格式化 PERMANOVA 结果
  p_value_text <- ifelse(otu_adonis$`Pr(>F)`[1] < 0.001, "< 0.001", 
                         sprintf("%.3f", otu_adonis$`Pr(>F)`[1]))
  
  # 存储结果
  result <- data.frame(
    time_point = time_point,
    R2 = otu_adonis$R2[1],
    p_value = otu_adonis$`Pr(>F)`[1],
    p_value_text = p_value_text,
    df = otu_adonis$Df[1],
    F_value = otu_adonis$F[1],
    adonis_text = sprintf("PERMANOVA:\nR² = %.3f\np = %s", 
                          otu_adonis$R2[1], p_value_text)
  )
  
  permanova_results <- rbind(permanova_results, result)
  
  cat("Time:", time_point, 
      "| R² =", round(otu_adonis$R2[1], 3), 
      "| p =", p_value_text, "\n")
}

# 将PERMANOVA结果合并到绘图数据中
all_plot_data <- all_plot_data %>%
  left_join(permanova_results %>% select(time_point, adonis_text), by = c("time" = "time_point"))

# 为每个分面计算文本位置
text_positions <- all_plot_data %>%
  group_by(time) %>%
  summarise(
    x_min = min(MDS1, na.rm = TRUE),
    x_max = max(MDS1, na.rm = TRUE),
    y_min = min(MDS2, na.rm = TRUE),
    y_max = max(MDS2, na.rm = TRUE)
  ) %>%
  mutate(
    # 将文本放在每个分面的右上角（相对位置）
    text_x = x_min + (x_max - x_min) * 0.7,
    text_y = y_min + (y_max - y_min) * 0.8
  )

# 合并文本位置
all_plot_data <- all_plot_data %>%
  left_join(text_positions %>% select(time, text_x, text_y), by = "time")

### 绘制分面图
# 创建分面图
combined_plot <- ggplot(data=all_plot_data, aes(x=MDS1, y=MDS2, color=treatment)) +
  theme_bw() +
  geom_point(size=2.5, shape=16, alpha=0.8) +
  facet_wrap(~ time, ncol = 2) +  # 按时间分面
  geom_vline(xintercept = 0, linetype="dashed", colour="#BEBEBE", linewidth=0.5) +
  geom_hline(yintercept = 0, linetype="dashed", colour="#BEBEBE", linewidth=0.5) +
  stat_ellipse(geom = "polygon", level=0.90, linewidth=NA,
               aes(fill=treatment), alpha=0.2, show.legend = TRUE) +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  xlab("NMDS1") + ylab("NMDS2") +
  # 添加全局stress值到标题
  labs(title = paste("Global NMDS Stress =", all_stress)) +
  # 为每个分面添加对应的PERMANOVA结果
  geom_text(aes(x = text_x, y = text_y, label = adonis_text),
            size = 3, color = "black", inherit.aes = T, hjust = 0.6, vjust = -0.4) +
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.box = "horizontal",
    legend.margin = margin(t = -5, b = 5),
    axis.title.x = element_text(size=10),
    axis.title.y = element_text(size=10, angle=90),
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(size=8),
    panel.grid = element_blank(),
    strip.text = element_text(size=10, face="plain"),
    strip.background = element_rect(fill="lightgray", color="black"),
    panel.spacing = unit(1.2, "lines")
  )

# 显示图形
print(combined_plot)

############################## [FigS7-B] Processing and Painting #####################################
# 计算全局NMDS
all_bray <- as.matrix(vegdist(otu, method = 'bray'))
all_nmds <- metaMDS(all_bray, k=2)
all_stress <- round(all_nmds$stress, 4)

# 准备绘图数据
all_plot_data <- all_nmds$points %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(group, by = "sample") %>%
  rename(MDS1 = MDS1, MDS2 = MDS2)

# 为每个diet计算PERMANOVA
diets <- unique(group$treatment)
permanova_results <- data.frame()

for(d in diets) {
  sub_data <- group[group$treatment == d, ]
  sub_otu <- otu[rownames(otu) %in% sub_data$sample, ]
  
  if(nrow(sub_otu) > 1) {
    sub_bray <- as.matrix(vegdist(sub_otu, method = 'bray'))
    time_df <- data.frame(sample = rownames(sub_otu))
    time_df <- merge(time_df, sub_data, by = "sample")
    otu_adonis <- adonis2(sub_bray ~ time, data = time_df, permutations = 999)
    
    p_value_text <- ifelse(otu_adonis$`Pr(>F)`[1] < 0.001, "< 0.001", 
                           sprintf("%.3f", otu_adonis$`Pr(>F)`[1]))
    
    result <- data.frame(
      diet = d,
      R2 = otu_adonis$R2[1],
      p_value = otu_adonis$`Pr(>F)`[1],
      p_value_text = p_value_text,
      adonis_text = sprintf("PERMANOVA:\nR² = %.3f\np = %s", 
                            otu_adonis$R2[1], p_value_text)
    )
    permanova_results <- rbind(permanova_results, result)
  }
}

# 合并PERMANOVA结果
all_plot_data <- all_plot_data %>%
  left_join(permanova_results %>% dplyr::select(diet, adonis_text), 
            by = c("treatment" = "diet"))

# 计算文本位置
text_positions <- all_plot_data %>%
  group_by(treatment) %>%
  summarise(
    x_min = min(MDS1, na.rm = TRUE),
    x_max = max(MDS1, na.rm = TRUE),
    y_min = min(MDS2, na.rm = TRUE),
    y_max = max(MDS2, na.rm = TRUE)
  ) %>%
  mutate(
    text_x = x_min + (x_max - x_min) * 0.7,
    text_y = y_min + (y_max - y_min) * 0.8
  )

all_plot_data <- all_plot_data %>%
  left_join(text_positions %>% dplyr::select(treatment, text_x, text_y), by = "treatment")

# 颜色设置
color_map <- c(CON_D0 = "grey", CON_D10 = "#D4E5F0", CON_D14 = "#AAC8DD", CON_D18 = "#8BBCD6",
               HFD_D0 = "grey", HFD_D10 = "#F9D6D2", HFD_D14 = "#F4A098", HFD_D18 = "#EE7C70",
               FUC_D0 = "grey", FUC_D10 = "#E5F0C1", FUC_D14 = "#C6E088", FUC_D18 = "#ACD26A")

# 绘图
all_plot_data$treatment <- factor(all_plot_data$treatment,levels = c("CON","HFD","FUC"))
ggplot(all_plot_data, aes(x = MDS1, y = MDS2, color = paste(treatment, time, sep = "_"))) +
  geom_point(size = 2.5, alpha = 0.8) +
  facet_wrap(~ treatment) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#BEBEBE", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#BEBEBE", linewidth = 0.5) +
  stat_ellipse(geom = "polygon", level = 0.90, linewidth = NA,
               aes(fill = paste(treatment, time, sep = "_")), alpha = 0.2) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  xlab("NMDS1") + ylab("NMDS2") +
  labs(title = paste("Global NMDS Stress =", all_stress)) +
  geom_text(aes(x = text_x, y = text_y, label = adonis_text),
            size = 3, color = "black", inherit.aes = FALSE, hjust = 0.6, vjust = -0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = "black"))
