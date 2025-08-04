##Script to Supplementary fig7-D.

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
viruse_contig <- read.table("../data/viruse_contig")
rongyuan_reads_forzhexian <- read.csv("../data/fig3_S7_reads.csv")
rongyuan_contig_length <- read.csv("../data/fig3_S7_length.csv")
viruse_IMG_contig_anno_chayijun <- read.csv("../data/fig3_S7_anno.csv")
x <- read.csv("../data/DABs_lefse.csv")

############################## Data processing ##############################
chayijun <- sapply(strsplit(x$Spe,split = "g__"),"[",2)
chayijun <- paste(chayijun,collapse = "|")

rongyuan_reads_forzhexian$Contig2 <- substring(rongyuan_reads_forzhexian$Contig,1,30)
colnames(rongyuan_contig_length) <- c("Contig","rongyuan_contig_length")
rongyuan_contig_length$Contig2 <- sapply(strsplit(rongyuan_contig_length$Contig,split ="-"),"[",1)
rongyuan_reads_forzhexian<- merge(rongyuan_reads_forzhexian,rongyuan_contig_length,by.x="Contig",by.y="Contig2")
rongyuan_reads_forzhexian$Normalized_Mapped_reads <- as.numeric(rongyuan_reads_forzhexian$reads_number)/((rongyuan_reads_forzhexian$rongyuan_contig_length)*as.numeric(rongyuan_reads_forzhexian$num_seqs))*10000000000
colnames(rongyuan_reads_forzhexian) <- paste("mapped",colnames(rongyuan_reads_forzhexian),sep="_")
viruse_contig$Contig_number2 <- substring(viruse_contig$Contig_number,1,30)
rongyuan_reads_forzhexian_chayi_forzhexian <- merge(rongyuan_reads_forzhexian,viruse_contig,by.x="mapped_Contig2",by.y="Contig_number2")
rongyuan_reads_forzhexian_chayi_forzhexian$chayi <- "Normal bacteria"
rongyuan_reads_forzhexian_chayi_forzhexian$chayi[grep(chayijun,rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host)] <- "Differential bacteria"
rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host,split = ";g__"),"[",2)
rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G,split = ";s__"),"[",1)
rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G,split = ";"),"[",1)
rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G <- sapply(strsplit(rongyuan_reads_forzhexian_chayi_forzhexian$Bacteria_Host_G,split = "_"),"[",1)

rongyuan_reads_chayi6<- rongyuan_reads_forzhexian_chayi_forzhexian[rongyuan_reads_forzhexian_chayi_forzhexian$chayi=="Differential bacteria",]
a <- aggregate(list(Contig_Abundance=rongyuan_reads_chayi6$mapped_Normalized_Mapped_reads),
               by=list(GD=rongyuan_reads_chayi6$mapped_GD,
                       Phacts_result=rongyuan_reads_chayi6$Phacts_result,
                       Diet=rongyuan_reads_chayi6$mapped_Diet),sum )


############################## [FigS7-D] Painting #####################################
color <- c("#CD776D","#DCDDDD")
names(color) <- c("Lytic","Temperate")
a$Diet <- factor(a$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC")   )

ggplot(a,aes(x=GD,y=Contig_Abundance,fill=Phacts_result,group=Phacts_result))+geom_area( position = "fill",alpha=0.8)+facet_wrap(.~Diet,)+
  geom_line(position = "fill",color="black")+
  scale_fill_manual(values = color  )+theme_bw()+xlab("")+ylab("")
ggsave("../figS7-D.pdf")
