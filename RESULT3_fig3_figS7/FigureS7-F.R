##Script to Supplementary fig7-F.

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
                       # Phacts_result=rongyuan_reads_chayi6$bac_result,
                       mapped_Contig2=rongyuan_reads_chayi6$mapped_Contig2,
                       Diet=rongyuan_reads_chayi6$mapped_Diet,
                       Bacteria_Host_G=rongyuan_reads_chayi6$Bacteria_Host_G),mean)
a1 <- unique(a[,c("Phacts_result","Bacteria_Host_G")])
a1 <- a1[duplicated(a1$Bacteria_Host_G),]
a <- a[a$Bacteria_Host_G%in%a1$Bacteria_Host_G,]
a <- spread(a, key= "Phacts_result",value = "Contig_Abundance")
a[is.na(a)] <- 0
a$bili <- as.numeric(a$Lytic+1)/as.numeric(a$Temperate+1)
a <- aggregate(list(bili=a$bili),
               by=list(GD=a$GD,
                       Diet=a$Diet,
                       Bacteria_Host_G=a$Bacteria_Host_G),median)
a$group <- paste(a$Diet,a$Bacteria_Host_G,sep="_")
# hist(a$bili)
a1 <- a[a$GD=="D0",]
a$bili_D0 <- factor(a$group,levels = a1$group,labels = a1$bili)
# a <- a[is.na(a$bili_D0)==F,]
a$bili_D0 <- as.numeric(as.character(a$bili_D0))
a[is.na(a)] <- 1
a$bili_change <- a$bili/as.numeric(as.character(a$bili_D0))
a1 <- a[order(a$Diet,a$GD),]
a1$Group <- paste(a1$Diet,a1$GD,sep="_")
a1 <- spread(a1[c("Group","bili_change","Bacteria_Host_G")], key="Group",value = "bili_change")
rownames(a1) <- a1[,1]
a1 <- a1[,-1]
result <- dist(a1, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")
groups<- cutree(result_hc, k=nrow(a1)/2 )
result_hc2 <- names(sort(groups))
a$Bacteria_Host_G <- factor(a$Bacteria_Host_G,levels = result_hc2,labels =result_hc2)
a$Diet <- factor(a$Diet,levels =c("CON","HFD","FUC") ,labels = c("CON","HFD","FUC"))


############################## [FigS7-F] Painting #####################################
ggplot(a,aes(x=GD,y=Bacteria_Host_G,fill=log10(bili_change)))+
  geom_point(shape=21,size=4)+
  theme_test()+
  facet_wrap(.~Diet)+
  scale_fill_gradient2(low = "#2C7FB8",mid = "#F7F8F8",high = "#D7301F")
ggsave("../figS7-F.pdf")
