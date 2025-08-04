##Script to Figure 6-D.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","dplyr","ggridges","ggsignif","viridis","ggpubr","RColorBrewer","ggplot2")
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
prj_data <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv")

AMGs_bidui2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_AMGs_bidui2.csv")
AMGs_out_id <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_AMGs_id.csv")
VIBRANT_genbank_table <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_genbank.csv")
VIBRANT_AMG_pathways <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_AMG_pathways.csv")
VIBRANT_AMG_counts <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6D_VIBRANT_AMG_counts.csv")

############################## Data processing ##############################
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

a<- aggregate(list(Mapped_number_nomalized= AMGs_bidui3$Mapped_number_nomalized ),
                         by=list(BioProject=AMGs_bidui3$BioProject,
                                 Metabolism=AMGs_bidui3$Metabolism,
                                 Nutrition1=AMGs_bidui3$Nutrition1),
                         median )
a$Type <- as.character(factor(a$Nutrition1,
                 levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR"),
                 labels = c("case","control","control","case","case","case")))
AMGs_bidui3$Type <- as.character(factor(AMGs_bidui3$Nutrition1,
                                        levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR"),
                                        labels = c("case","control","control","case","case","case")))

AMG_sampleAb <- aggregate(list(Mapped_number_nomalized= AMGs_bidui3$Mapped_number_nomalized ),
                          by=list(BioProject=AMGs_bidui3$BioProject,
                                  Sample=AMGs_bidui3$SRA,
                                  Metabolism=AMGs_bidui3$Metabolism,
                                  Nutrition1=AMGs_bidui3$Nutrition1),
                          sum ) #等价于使用dplyr的group_by和summarise函数进行聚合求和
AMG_sampleAb$Type <- as.character(factor(AMG_sampleAb$Nutrition1,
                                         levels = c("APS", "Baseline" ,"Control",  "FOS", "Inulin","NR"),
                                         labels = c("case","control","control","case","case","case")))
##用丰度简单看一下差异性和效果
R_paintdata <- filter(AMG_sampleAb, Metabolism %in% c("Amino acid metabolism", 
                                                      "Biosynthesis of other secondary metabolites",
                                                      "Carbohydrate metabolism",
                                                      "Energy metabolism",
                                                      "Folding, sorting and degradation",
                                                      "Glycan biosynthesis and metabolism",
                                                      "Lipid metabolism",
                                                      "Metabolism of cofactors and vitamins",
                                                      "Metabolism of other amino acids")) #%in% 运算符用于只过滤出包含向量中提供的数据的列
R_paintdata$Metabolism <- factor(R_paintdata$Metabolism, levels = c("Carbohydrate metabolism",
                                                                    "Amino acid metabolism", 
                                                                    "Metabolism of cofactors and vitamins",
                                                                    "Energy metabolism",
                                                                    "Lipid metabolism",
                                                                    "Glycan biosynthesis and metabolism",
                                                                    "Biosynthesis of other secondary metabolites",
                                                                    "Metabolism of other amino acids",
                                                                    "Folding, sorting and degradation"))#按大致丰度排序
R_paintdata$Type <- factor(R_paintdata$Type, levels = c("control","case"))
R_paintdata$BioProject <- factor(R_paintdata$BioProject, levels = c("PRJNA704567","PRJNA615253","PRJNA731974","mgp6153")) #相当于按4种多糖分组

############################## [Fig 6-D] Painting ##############################
paintdata_summary <- R_paintdata %>%
  group_by(BioProject,Metabolism,Type) %>%
  summarise(
    mean = mean(Mapped_number_nomalized, na.rm = TRUE),
    sd = sd(Mapped_number_nomalized, na.rm = TRUE),
    n = n(),
    sem = sd / sqrt(n),
    ymin = mean - sd,
    ymax = mean + sd
  )
ggplot(paintdata_summary, aes(x=BioProject, y=mean))+
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color=Type), position=position_dodge(width = 0.45))+
  scale_color_manual(values = c("#91CAE8","#F48892"))+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.text.x = element_text(angle = 45))+
  facet_wrap(.~Metabolism,nrow=2, scales = "free") #只有mgp这个项目具有显著性

ggsave("../fig6-D.pdf")

