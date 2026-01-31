##Script to Supplementary fig10-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("dplyr","ggpubr","ggplot2")
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

library(ggplot2)
library(dplyr)
library(ggpubr)
############################## Data loading ##############################
bacphlip <- read.csv("~/cooperation/202409zhaofengxiang/data/figS10C_bacphlip.results",sep="\t")
clust <- read.table("~/cooperation/202409zhaofengxiang/data/figS10C_clust.tsv"
                    ,col.names = c("Cluster","Contigs"))
prj<- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv",header=T,fileEncoding = 'GBK')

############################ Data processing ##############################
bacphlip <- bacphlip[bacphlip$X!="X",]
bacphlip <- bacphlip[!duplicated(bacphlip$X),]
rownames(bacphlip) <- bacphlip$X
bacphlip <- data.frame(Contigs=rownames(bacphlip),bacphlip)
bacphlip$Samples <- sapply(strsplit(bacphlip$Contigs,split = ".NODE",fixed = T),"[",1)
bacphlip$Cov <- sapply(strsplit(bacphlip$Contigs,split = "cov_",fixed = T),"[",2)
bacphlip$Cov <- sapply(strsplit(bacphlip$Cov,split = "||",fixed = T),"[",1)
bacphlip$Type <- "Lytic"
bacphlip$Type[bacphlip$Temperate>=0.5] <- "Temperate"


clust <- merge(clust,bacphlip,by.x = "Contigs",by.y="Contigs")
clust <- merge(clust,prj,by.x = "Samples",by.y = "Run")
rongyuan_TT <- aggregate(list(Contig_Abundance=as.numeric(clust$Cov)),
                         by=list(Nutrition1=clust$Nutrition1,
                                 BioProject=clust$BioProject,
                                 Samples=clust$Samples,
                                 PhageLifestyle=clust$Type,
                                 Viruse_Family=clust$Cluster
                                 
                         ), median)
# rongyuan_TT2 <- pivot_wider(rongyuan_TT,)
rongyuan_TT2 <- aggregate(list(Contig_Abundance=clust$Contigs),
                          by=list(Viruse_Family=clust$Cluster,
                                  Nutrition1=clust$Nutrition1,
                                  PhageLifestyle=clust$Type,
                                  BioProject=clust$BioProject
                          ), length)
rongyuan_TT2 <- rongyuan_TT2[rongyuan_TT2$Contig_Abundance>=2,]

rongyuan_TT2$a <- paste(rongyuan_TT2$Viruse_Family,rongyuan_TT2$Nutrition1,rongyuan_TT2$PhageLifestyle,sep="\t")
rongyuan_TT$a <- paste(rongyuan_TT$Viruse_Family,rongyuan_TT$Nutrition1,rongyuan_TT$PhageLifestyle,sep="\t")
# rongyuan_TT2 <- rongyuan_TT2[rongyuan_TT2$a>=2,]
rongyuan_TT <- rongyuan_TT[rongyuan_TT$a%in%rongyuan_TT2$a,]
rongyuan_TT <- rongyuan_TT[is.na(rongyuan_TT$Contig_Abundance)==F,]


rongyuan_TT <- pivot_wider(rongyuan_TT,names_from = "PhageLifestyle",values_from = "Contig_Abundance")
value=0.1
rongyuan_TT$Lytic[is.na(rongyuan_TT$Lytic)] <- 0
rongyuan_TT$Temperate[is.na(rongyuan_TT$Temperate)] <- 0
rongyuan_TT$Temperate_chu_Lytic <- (as.numeric(rongyuan_TT$Lytic)+value)/(as.numeric(rongyuan_TT$Temperate)+value)
rongyuan_TT$Temperate_chu_Lytic <-abs(log10(rongyuan_TT$Temperate_chu_Lytic))
# rongyuan_TT$Temperate_chu_Lytic <-log10(rongyuan_TT$Temperate_chu_Lytic)
rongyuan_TT <- rongyuan_TT[is.na(rongyuan_TT$Temperate_chu_Lytic)==F,]

############################## [FigS10-C] Painting #####################################
paintdata <- rongyuan_TT %>%
  filter(!grepl("PRJNA704567",BioProject)) %>%
  filter(!grepl("probiotics",Nutrition1)) %>%
  mutate(Group = paste0(BioProject,"_",Nutrition1)) %>%
  mutate(BioProject = factor(BioProject,levels = c("PRJNA615253","PRJNA731974","mgp6153"))) %>%
  mutate(Nutrition1 = factor(Nutrition1,levels = c("Control","APS","Inulin","FOS"))) %>%
  arrange(BioProject, Nutrition1) %>%
  mutate(Group = factor(Group)) 

color_Group <- c("#DDDDDD","#2F5B60",
                 "#DDDDDD","#AAD3C7",
                 "#DDDDDD","#9CB4AD")
names(color_Group)<- unique(paintdata$Group)

stat.test <- paintdata %>%
  group_by(BioProject) %>%
  wilcox_test(Temperate_chu_Lytic ~ Nutrition1) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Nutrition1", dodge = 0.01, group = "BioProject") %>%
  mutate(xmin=1,xmax=1+xmin)

ggplot(paintdata,aes(x=Nutrition1,y=Temperate_chu_Lytic,fill=Group))+
  geom_boxplot(position = position_dodge(width=0.01,preserve = "single"), width=0.45,
               outlier.colour = NA,alpha=0.8)+
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values =color_Group)+
  xlab("")+
  ylab("LLI")+
  facet_wrap(~ BioProject, scales = "free")+
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  # 显示显著性符号：*，**，***
    tip.length = 0.03,
    #step.increase = 0.01,  # 调整标记位置，避免重叠
    size = 5
  ) 

