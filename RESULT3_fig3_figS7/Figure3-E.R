##Script to Figure 3-E.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","dplyr","FSA","rstatix","RColorBrewer","ggplot2")
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
viruse_contig <- read.table("../data/viruse_contig")
rongyuan_reads <- read.csv("../data/fig3_S7_reads.csv")
rongyuan_contig_length <- read.csv("../data/fig3_S7_length.csv")
viruse_IMG_contig_anno_chayijun <- read.csv("../data/fig3_S7_anno.csv")


############################## Data processing ##############################
rongyuan_reads$Contig2 <- substring(rongyuan_reads$Contig,1,30)
colnames(rongyuan_contig_length) <- c("Contig","rongyuan_contig_length")
rongyuan_contig_length$Contig2 <- sapply(strsplit(rongyuan_contig_length$Contig,split ="-"),"[",1)
rongyuan_reads<- merge(rongyuan_reads,rongyuan_contig_length,by.x="Contig",by.y="Contig2")
rongyuan_reads$Normalized_Mapped_reads <- as.numeric(rongyuan_reads$reads_number)/(as.numeric(rongyuan_reads$rongyuan_contig_length)*as.numeric(rongyuan_reads$num_seqs))*10000000000
colnames(rongyuan_reads) <- paste("mapped",colnames(rongyuan_reads),sep="_")

rongyuan_reads_chayi <- merge(rongyuan_reads,viruse_IMG_contig_anno_chayijun,by.x="mapped_Contig2",by.y="Contig_number2")
rongyuan_reads_chayi$fill <- paste(rongyuan_reads_chayi$Bacteria_Host_G,rongyuan_reads_chayi$GD_mapped,sep="_")

a <- rongyuan_reads_chayi[c("mapped_Normalized_Mapped_reads","mapped_GD","Phacts_result","mapped_Contig2",
                            "mapped_Diet","Bacteria_Host_G","mapped_Sample_id")]
a <- aggregate(list(Contig_Abundance=rongyuan_reads_chayi$mapped_Normalized_Mapped_reads),
               by=list(Time=rongyuan_reads_chayi$mapped_GD,
                       Phacts_result=rongyuan_reads_chayi$Phacts_result,
                       # Phacts_result=rongyuan_reads_chayi6$bac_result,
                       mapped_Contig2=rongyuan_reads_chayi$mapped_Contig2,
                       Diet=rongyuan_reads_chayi$mapped_Diet,
                       Bacteria_Host_G=rongyuan_reads_chayi$Bacteria_Host_G),mean)
a1 <- unique(a[,c("Phacts_result","Bacteria_Host_G")])
a1 <- a1[duplicated(a1$Bacteria_Host_G),]
a <- a[a$Bacteria_Host_G%in%a1$Bacteria_Host_G,]
a <- spread(a, key= "Phacts_result",value = "Contig_Abundance")
a[is.na(a)] <- 0
a$bili <- as.numeric(a$Lytic+1)/as.numeric(a$Temperate+1)
analysisdata <- a
analysisdata$log10_bili <- log10(analysisdata$bili)
#write.csv(analysisdata,"~/cooperation/202409zhaofengxiang/code/fig3-E_analysisdata.csv")
a <- aggregate(list(bili=a$bili),
               by=list(Time=a$Time,
                       Diet=a$Diet,
                       Bacteria_Host_G=a$Bacteria_Host_G),median)
a <- a[a$Time%in%c("D0","D10"),]
rongyuan_reads_chayi4<- a
rongyuan_reads_chayi4 <- spread(rongyuan_reads_chayi4[,c("Bacteria_Host_G","Diet",'Time','bili')], key= 'Time', value =  'bili')
rongyuan_reads_chayi4 <- rongyuan_reads_chayi4[is.na(rongyuan_reads_chayi4$D0)==F&is.na(rongyuan_reads_chayi4$D10)==F,]
rongyuan_reads_chayi4$Change <- as.numeric(rongyuan_reads_chayi4$D10)/as.numeric(rongyuan_reads_chayi4$D0)
rongyuan_reads_chayi4$Change2[as.numeric(rongyuan_reads_chayi4$Change)<0] <- -1*(log10(abs(rongyuan_reads_chayi4$Change[as.numeric(rongyuan_reads_chayi4$Change)<0])))
rongyuan_reads_chayi4$Change2[as.numeric(rongyuan_reads_chayi4$Change)>0] <- log10(abs(rongyuan_reads_chayi4$Change[as.numeric(rongyuan_reads_chayi4$Change)>0]))
rongyuan_reads_chayi4 <-rongyuan_reads_chayi4[order(rongyuan_reads_chayi4$Diet,rongyuan_reads_chayi4$Change),]
rongyuan_reads_chayi4$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_chayi4$Bacteria_Host_G,split = ";"),"[",1)
rongyuan_reads_chayi4$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_chayi4$Bacteria_Host_G,split = "_",fixed = T),"[",1)
rongyuan_reads_chayi4$Bacteria_Host_G <- paste(rongyuan_reads_chayi4$Bacteria_Host_G,"Phage",sep="_")
rongyuan_reads_chayi4$Bacteria_Host_G <- factor(rongyuan_reads_chayi4$Bacteria_Host_G,levels = unique(rongyuan_reads_chayi4$Bacteria_Host_G),labels = unique(rongyuan_reads_chayi4$Bacteria_Host_G))

#浮点数精度问题
rongyuan_reads_chayi_paint <- rongyuan_reads_chayi4 %>% 
  filter(!(Bacteria_Host_G == "Helicobacter_Phage" & 
             near(Change2, -0.396526240, tol = 1e-8))) %>%
  filter(!(Bacteria_Host_G == "Helicobacter_Phage" & 
             near(Change2, 0.041348005, tol = 1e-8)))  %>%
  filter(!(Bacteria_Host_G == "Helicobacter_Phage" & 
             near(Change2, 0.282861638, tol = 1e-8))) 
#更改diet顺序
rongyuan_reads_chayi_paint$Diet <- as.factor(rongyuan_reads_chayi_paint$Diet)
rongyuan_reads_chayi_paint$Diet <- factor(rongyuan_reads_chayi_paint$Diet, levels = c("CON","HFD","FUC"))
#更改Y轴显示顺序
con_data <- rongyuan_reads_chayi_paint %>%
  filter(Diet == "CON") %>%
  arrange(Change2)   # 默认从小到大排序
host_order <- con_data$Bacteria_Host_G # 获取排序后的Bacteria_Host_G向量
rongyuan_reads_chayi_paint$Bacteria_Host_G <- factor(rongyuan_reads_chayi_paint$Bacteria_Host_G, 
                                                     levels = host_order)


############################## [Fig3-E] Painting #####################################
ggplot(rongyuan_reads_chayi_paint,aes(x=Change2,
                                      y = Bacteria_Host_G,
                                      fill=Diet))+
  geom_bar(stat = "identity")+facet_wrap(.~Diet,nrow = 1)+theme_bw()+
  scale_fill_manual(values =  c("#8BBCD6","#DD7F88","#A0CC58") )+
  xlab("Phage lysis and lysogenic ratio")+ylab("")+
  geom_vline(aes(xintercept=0,color=Diet),size=0.8)+
  scale_color_manual(values =  c("#8BBCD6","#DD7F88","#A0CC58") )+
  theme(panel.grid=element_blank())
ggsave("../fig3-E.pdf")


############################## Difference analysis #####################################
analysisdata$Diet <- as.factor(analysisdata$Diet)
# 按细菌宿主分组进行检验
results <- analysisdata %>%
  group_by(Bacteria_Host_G) %>%
  filter(n() >= 5) %>%  # 仅分析样本量≥5的宿主（确保检验效力）
  summarise(
    # Kruskal-Wallis检验（非参数版单因素ANOVA）：比较三组（CON/FUC/HFD）的bili中位数差异;不要求正态分布和方差齐性，适合小样本
    #kruskal_p = kruskal.test(bili ~ Diet)$p.value,
    kruskal_p = kruskal.test(log10_bili ~ Diet)$p.value,
    
    # Dunn事后检验（若整体显著）:两两比较具体哪组存在差异,使用Benjamini-Hochberg方法校正p值（控制假阳性）
    dunn_test = if (kruskal_p < 0.05) {
      list(dunnTest(bili ~ Diet, method = "bh")$res)
    } else {
      NA
    }
  ) %>%
  arrange(kruskal_p)  # 按p值排序

# 查看显著结果
significant_results <- results %>% 
  filter(kruskal_p < 0.05)

# 输出结果
print(significant_results, n = Inf)

# 提取具体宿主的Dunn检验结果（示例：Parabacteroides）
parabacteroides_dunn <- results %>%
  filter(Bacteria_Host_G == "Parabacteroides") %>%
  pull(dunn_test)
print(parabacteroides_dunn)

