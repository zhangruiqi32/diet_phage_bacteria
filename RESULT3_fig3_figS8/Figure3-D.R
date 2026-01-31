##Script to Figure 3-D.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyr","RColorBrewer","ggpubr","ggplot2")
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
viruse_contig <- read.table("~/cooperation/202409zhaofengxiang/data/viruse_contig")

############################## Data processing ##############################
viruse_contig2<- viruse_contig[grep("Temperate|Lytic",viruse_contig$Phacts_result),][grep("Temperate|Lytic",viruse_contig$Phacts_result),]
viruse_contig2<- viruse_contig2[is.na(viruse_contig2$GD)==F,]
viruse_contig2$Diet <- factor(viruse_contig2$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC")   )
a <- aggregate(list(Contig_Abundance=viruse_contig2$Contig_Abundance),
               by=list(GD=viruse_contig2$GD,
                       Phacts_result=viruse_contig2$Phacts_result,
                       Diet=viruse_contig2$Diet),sum )

############################## [Fig3-D] Painting #####################################
color <- c("#CD776D","#DCDDDD")
names(color) <- c("Lytic","Temperate")

ggplot(a,aes(x=GD,y=Contig_Abundance,fill=Phacts_result,group=Phacts_result))+geom_area( position = "fill",alpha=0.8)+facet_wrap(.~Diet,)+
  geom_line(position = "fill",color="black")+
  scale_fill_manual(values = color)+theme_bw()+
  theme(panel.grid=element_blank())+
  xlab("")+ylab("")
ggsave("../fig3-D.pdf")

