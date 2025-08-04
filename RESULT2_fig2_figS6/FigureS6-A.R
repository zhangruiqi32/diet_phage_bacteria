##Script to Supplementary fig6-A.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("reshape2","ggplot2","tidyr","dplyr","FSA")
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
weight_data <- read.csv("../data/figS6A_weight.csv", header = T)  

############################## Data processing #####################################
weight_data$week <- as.factor(weight_data$week)
# 计算每周各组体重的平均值（忽略缺失值）
weight_mean <- aggregate(list(CON = weight_data$CON, FUC = weight_data$FUC, HFD = weight_data$HFD),
                by = list(week = weight_data$week),
                function(x) mean(x[is.na(x) == F]))  # 使用is.na(x)==F过滤缺失值

# 将宽格式数据转换为长格式（便于ggplot绘图）
weight_mean <- melt(weight_mean, 
           variable.name = "Group",  # 新列名用于存储组别
           value.name = "mean")   # 新列名用于存储平均值

# 计算每周各组体重的标准差（忽略缺失值）
weight_sd <- aggregate(list(CON = weight_data$CON, FUC = weight_data$FUC, HFD = weight_data$HFD),
                by = list(week = weight_data$week),
                function(x) sd(x[is.na(x) == F]))  # 计算标准

# 将标准差数据转换为长格式
weight_sd <- melt(weight_sd, 
           variable.name = "Group", 
           value.name = "sd")  # 新列名用于存储标准差

# 将标准差合并到平均值数据框中
weight_mean <- data.frame(weight_mean, sd = weight_sd[, 3])  # 添加sd列

# 计算误差范围（平均值±标准差）
weight_mean$max <- weight_mean$mean + weight_mean$sd  # 误差上限
weight_mean$min <- weight_mean$mean - weight_mean$sd  # 误差下限

# 将Group列转换为因子类型（ggplot绘图需要）
weight_mean$Group <- as.factor(weight_mean$Group)


############################## [FigS6-A] Painting #####################################
ggplot(data = weight_mean, aes(x = week, y = mean,  
                      group = Group, color = Group, shape = Group)) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),  # 添加误差棒
                width = 0.2) +  # 误差棒宽度
  geom_point(size = 3) +         # 添加数据点，大小=3
  geom_line(size = 1) +          # 添加连线，线宽=1
  labs(x = "week", y = "Body weight") +  # 设置坐标轴标签
  theme_classic() +              # 使用经典主题（简洁无背景网格）
  scale_color_manual(values = c("#8BBCD6", "#ACD26A", "#DD7F88"))  # 自定义颜色

ggsave("../figS6-A.pdf") 

############################## Difference analysis #####################################
# 提取第6周数据
week6 <- subset(weight_data, week == 6)

# 创建长格式数据
long_data <- pivot_longer(week6, 
                          cols = c(CON, FUC, HFD),
                          names_to = "group",
                          values_to = "weight",
                          values_drop_na = TRUE)

# 1. 正态性检验 (Shapiro-Wilk)
shapiro_CON <- shapiro.test(long_data$weight[long_data$group == "CON"])
shapiro_FUC <- shapiro.test(long_data$weight[long_data$group == "FUC"])
shapiro_HFD <- shapiro.test(long_data$weight[long_data$group == "HFD"])

# 2. 方差齐性检验 (Bartlett检验)
bartlett_result <- bartlett.test(weight ~ group, data = long_data)

# 3. 描述性统计
desc_stats <- long_data %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(weight),
    sd = sd(weight)
  )

# 4. 参数检验选择
if(all(c(shapiro_CON$p.value, shapiro_FUC$p.value, shapiro_HFD$p.value) > 0.05) &&
   bartlett_result$p.value > 0.05) {
  # 满足参数检验条件
  anova_result <- aov(weight ~ group, data = long_data)
  anova_summary <- summary(anova_result)
  p_value_anova <- anova_summary[[1]]$`Pr(>F)`[1]
  
  # 事后检验 (Tukey HSD)
  if(p_value_anova < 0.05) {
    tukey_result <- TukeyHSD(anova_result)
  }
} else {
  # 不满足参数检验条件，使用非参数检验
  kruskal_result <- kruskal.test(weight ~ group, data = long_data)
  
  # 事后检验 (Dunn检验)
  if(kruskal_result$p.value < 0.05) {
    library(FSA)
    dunn_result <- dunnTest(weight ~ group, data = long_data, method = "bonferroni")
  }
}

# 打印所有结果
print(desc_stats)
if(exists("anova_result")) {
  cat("\nANOVA results:\n")
  print(anova_summary)
  if(exists("tukey_result")) {
    cat("\nTukey HSD post-hoc tests:\n")
    print(tukey_result)
  }
} else if(exists("kruskal_result")) {
  cat("\nKruskal-Wallis test results:\n")
  print(kruskal_result)
  if(exists("dunn_result")) {
    cat("\nDunn's post-hoc tests:\n")
    print(dunn_result)
  }
}
