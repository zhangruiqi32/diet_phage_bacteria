##Script to Supplementary fig6-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("reshape2","ggpubr","RColorBrewer","ggplot2")
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
x_forp2_gg <- read.csv("~/cooperation/202409zhaofengxiang/data/figS6B_1.csv")
lefse1 <- read.csv("~/cooperation/202409zhaofengxiang/data/DABs_lefse.csv",sep=",",header = T )

############################## [FigS6-B] Processing and Painting #####################################
lefse1$Spe <- gsub("\\.","\\|",lefse1$Spe,fixed = F)
x_forp2_gg$Spe <- factor(x_forp2_gg$Spe,levels = lefse1$Spe,labels = lefse1$Spe)
value <- c("#88AFCB","#D2878D")
p1<- ggplot(x_forp2_gg,aes(x=Spe,y=Group,color=Change,fill=Change,size=Bacteria))+
  geom_point(alpha=1,shape = 21)+theme_bw()+
  scale_color_manual(values = value)+
  scale_fill_manual(values =value)+
  scale_size_continuous(range = c(3,8))+
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+
  xlab("Differentiated bacteria in Three Groups")+ylab("")

annotation_colors <- c("#7CA5C2","#85C3B9","#DDEAC4","#DA87B6","#E9C4DA","#F5D23B","#E1C192")
x_forp2_gg2 <- x_forp2_gg[!duplicated(x_forp2_gg$Spe),]

p2<- ggplot(x_forp2_gg2,aes(x=Spe,y=1,fill=Phylum,color=Phylum))+geom_bar(stat="identity",position="dodge",width=1)+
  scale_fill_manual(values = annotation_colors)+theme_bw()+
  scale_color_manual(values = annotation_colors)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    legend.position="bottom"
  )+xlab("")+ylab("")

ggarrange(p2,p1,ncol = 1)
ggsave("../figS6-B.pdf") 


