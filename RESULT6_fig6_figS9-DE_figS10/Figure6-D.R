##Script to Figure 6-D.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","dplyr","ggridges","ggsignif","viridis","ggpubr","RColorBrewer","ggplot2")
bioconductor_packages=c()

for(package in cran_packages){
    library(package, character.only=T,quietly = T, warn.conflicts = F)
}

if (!require("BiocManager", quietly = T, warn.conflicts = F)) install.packages("BiocManager")

for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}


library(tidyr)
library(dplyr)
library(ggplot2)
############################## Data loading ##############################
prj_data <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv")

AMGs_bidui2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_AMGs_bidui2.csv")
AMGs_out_id <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_AMGs_id.csv")
VIBRANT_genbank_table <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_genbank.csv")
VIBRANT_AMG_pathways <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_AMG_pathways.csv")
VIBRANT_AMG_counts <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_AMG_counts.csv")

AMG_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5_AMG_reads.csv")
AMG_contents_metabolism2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig5C_metabolism.csv")


############################## Data processing ##############################
AMG_reads <- merge(AMG_reads,AMG_contents_metabolism2[c("protein","metabolism")],by.x="mapped_Contig",by.y="protein")

length(intersect(AMGs_out_id$V1,unique(AMGs_bidui2$AMG_name)))
AMGs_bidui3 <- AMGs_bidui2[AMGs_bidui2$AMG_name%in%AMGs_out_id$V1,]
prjs <- unique(prj_data$BioProject[prj_data$Nutrition%in%c("Polysaccharide","Oligosaccharide")] )
AMGs_bidui3 <- AMGs_bidui3[AMGs_bidui3$BioProject%in%prjs,]

VIBRANT_genbank_table <- VIBRANT_genbank_table[VIBRANT_genbank_table$protein!="protein",]
VIBRANT_AMG_pathways <- VIBRANT_AMG_pathways[VIBRANT_AMG_pathways$KEGG.Entry!="KEGG Entry",]
VIBRANT_AMG_counts<- VIBRANT_AMG_counts[VIBRANT_AMG_counts$AMG.count!="AMG count",]

VIBRANT_AMG_pathways2<- aggregate(VIBRANT_AMG_pathways$Present.AMG.KOs,by=list(Metabolism=VIBRANT_AMG_pathways$Metabolism,Pathway=VIBRANT_AMG_pathways$Pathway),
                                  function(x) paste(unique(unlist(strsplit(x,split = ","))),collapse = ","))
colnames(VIBRANT_AMG_pathways2) <- c("Metabolism","Pathway","KO")
VIBRANT_AMG_pathways2 <- data.frame(Metabolism= rep(VIBRANT_AMG_pathways2$Metabolism, lengths(strsplit(VIBRANT_AMG_pathways2$KO, ","))),
                                    Pathway= rep(VIBRANT_AMG_pathways2$Pathway, lengths(strsplit(VIBRANT_AMG_pathways2$KO, ","))),
                                    KO=unlist(strsplit(VIBRANT_AMG_pathways2$KO, ",")))
VIBRANT_genbank_table <- merge(VIBRANT_genbank_table,VIBRANT_AMG_pathways2,by.x="accession",by.y = "KO")
AMGs_bidui3 <- merge(AMGs_bidui3 ,VIBRANT_genbank_table,by.x = "AMG_name",by.y="protein")
AMGs_bidui3$Mapped_number_nomalized <- AMGs_bidui3$Mapped_number/(as.numeric(AMGs_bidui3$num_seqs_before)*AMGs_bidui3$AMG_length*2)*1000000000 #FPKM计算（双端要/2）因为比对上该基因的reads/（所有能比对上的reads*该基因长度）→因为不同基因长度不同，比对上的概率也不同

AMGs_bidui3 <- AMGs_bidui3[AMGs_bidui3$Nutrition1!="probiotics",]
#AMG丰度改变对比的是 多糖+高脂 与 高脂 ！！！！需要Diet==High_fat_diet 
AMGs_bidui3 <- AMGs_bidui3[AMGs_bidui3$Diet =="High_fat_diet",]
AMGs_bidui3 <- AMGs_bidui3[AMGs_bidui3$BioProject!="PRJNA704567",]
AMGs_bidui3 <- AMGs_bidui3[AMGs_bidui3$Mapped_number>=10,]




# 处理AMG_reads数据
AMG_reads$Group <- as.character(AMG_reads$mapped_Group)
AMG_reads$Type <- ifelse(
  grepl("^HFD_D(10|14|18)$", AMG_reads$Group), "control",
  ifelse(grepl("^FUC_D(10|14|18)$", AMG_reads$Group), "case", NA)
)

# 筛选有效数据
AMG_reads_filtered <- AMG_reads[!is.na(AMG_reads$Type) & 
                                  AMG_reads$Group %in% c("HFD_D10","HFD_D14","HFD_D18","FUC_D10","FUC_D14","FUC_D18"), ]

# 添加必要列
AMG_reads_filtered$BioProject <- "Fucoidan"
AMG_reads_filtered$SRA <- AMG_reads_filtered$mapped_Sample_id
AMG_reads_filtered$Metabolism <- AMG_reads_filtered$metabolism
AMG_reads_filtered$Nutrition1 <- AMG_reads_filtered$Type
AMG_reads_filtered$Mapped_number_nomalized <- AMG_reads_filtered$mapped_Normalized_Mapped_reads

# 准备AMGs_bidui3数据
AMGs_bidui3$Mapped_number_nomalized <- AMGs_bidui3$Mapped_number_nomalized
AMGs_bidui3$Type <- as.character(factor(AMGs_bidui3$Nutrition1,
                                        levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR","case","control"),
                                        labels = c("case","control","control","case","case","case","case","control")))

# 合并数据
common_cols <- c("BioProject", "SRA", "Metabolism", "Nutrition1", "Mapped_number_nomalized", "Type")
AMG_combined <- rbind(AMGs_bidui3[, common_cols], 
                      AMG_reads_filtered[, c("BioProject", "SRA", "Metabolism", "Nutrition1", "Mapped_number_nomalized", "Type")])

# 聚合数据
AMG_sampleAb <- aggregate(list(Mapped_number_nomalized = AMG_combined$Mapped_number_nomalized),
                          by = list(BioProject = AMG_combined$BioProject,
                                    Sample = AMG_combined$SRA,
                                    Metabolism = AMG_combined$Metabolism,
                                    Nutrition1 = AMG_combined$Nutrition1),
                          sum)

AMG_sampleAb$Type <- as.character(factor(AMG_sampleAb$Nutrition1,
                                         levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR", "case", "control"),
                                         labels = c("case","control","control","case","case","case","case","control")))

# 筛选特定代谢通路
R_paintdata <- filter(AMG_sampleAb, Metabolism %in% c("Amino acid metabolism", 
                                                      "Biosynthesis of other secondary metabolites",
                                                      "Carbohydrate metabolism",
                                                      "Energy metabolism",
                                                      "Folding, sorting and degradation",
                                                      "Glycan biosynthesis and metabolism",
                                                      "Lipid metabolism",
                                                      "Metabolism of cofactors and vitamins",
                                                      "Metabolism of other amino acids"))

R_paintdata$Metabolism <- factor(R_paintdata$Metabolism, 
                                 levels = c("Carbohydrate metabolism",
                                            "Amino acid metabolism", 
                                            "Metabolism of cofactors and vitamins",
                                            "Energy metabolism",
                                            "Lipid metabolism",
                                            "Glycan biosynthesis and metabolism",
                                            "Biosynthesis of other secondary metabolites",
                                            "Metabolism of other amino acids",
                                            "Folding, sorting and degradation"))

R_paintdata$Type <- factor(R_paintdata$Type, levels = c("control","case"))
R_paintdata$Mapped_number_nomalized_abundance<- R_paintdata$Mapped_number_nomalized
R_paintdata$Mapped_number_nomalized <- log10(R_paintdata$Mapped_number_nomalized+1)
R_paintdata$BioProject <- factor(R_paintdata$BioProject, 
                                 levels = c("PRJNA704567","PRJNA615253","PRJNA731974","mgp6153","Fucoidan"))


############################## [Fig 6-D] Painting ##############################
R_paintdata$BioProject <- factor(R_paintdata$BioProject,
                                 levels = c("PRJNA704567","PRJNA615253","PRJNA731974","mgp6153","Fucoidan"),
                                 labels = c("Control","APS","Inulin","FOS","FUC"))

# 安全计算p值
sig_data <- R_paintdata %>%
  group_by(BioProject, Metabolism) %>%
  summarise(
    control_n = sum(Type == "control", na.rm = TRUE),
    case_n = sum(Type == "case", na.rm = TRUE),
    p = if(control_n >= 2 & case_n >= 2) {
      t.test(Mapped_number_nomalized[Type == "control"], 
             Mapped_number_nomalized[Type == "case"])$p.value
    } else NA
  ) %>%
  mutate(
    signif = ifelse(!is.na(p) & p < 0.05, "*", ""),
    y_pos = NA
  )

# 为有显著性的数据添加y位置
for(i in 1:nrow(sig_data)) {
  if(sig_data$signif[i] == "*") {
    sub_data <- R_paintdata[R_paintdata$BioProject == sig_data$BioProject[i] & 
                              R_paintdata$Metabolism == sig_data$Metabolism[i], ]
    sig_data$y_pos[i] <- max(sub_data$Mapped_number_nomalized, na.rm = TRUE) * 1.1
  }
}

# 绘图数据
plot_data <- R_paintdata %>%
  group_by(BioProject, Metabolism, Type) %>%
  summarise(
    mean = mean(Mapped_number_nomalized, na.rm = TRUE),
    sd = sd(Mapped_number_nomalized, na.rm = TRUE),
    ymin = mean - sd,
    ymax = mean + sd
  )

# 绘图
ggplot(plot_data, aes(x = BioProject, y = mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = Type), 
                  position = position_dodge(0.45)) +
  geom_text(data = sig_data[sig_data$signif == "*", ], 
            aes(x = BioProject, y = y_pos*0.91, label = signif), 
            size = 6, vjust = 0) +
  geom_segment(data = sig_data[sig_data$signif == "*", ], 
               aes(x = as.numeric(BioProject) - 1.2, 
                   xend = as.numeric(BioProject) -0.8,
                   y = y_pos * 0.92, 
                   yend = y_pos * 0.92), 
               size = 0.5) +
  scale_color_manual(values = c("#91CAE8", "#F48892")) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  facet_wrap(~Metabolism, nrow = 2, scales = "free_y") +
  labs(x = "Treatment", y = "log10(Normalized Mapped Reads + 1)", color = "Group")
