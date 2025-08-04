##Script to Figure 3-B and Supplementary fig7-AB.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggtern","RColorBrewer","ggpubr","ggplot2")
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
reads_number2s <- read.csv("../data/fig3B_S7AB_readnum.csv")
contigs_phacts <- read.csv("../data/fig3B_S7AB_lifestyle.csv")
cluster <- read.csv("../data/fig3B_S7AB_cluster.csv")
clstr_name_change <- read.csv("../data/fig3B_S7AB_name.csv")
Cluster_Linage <- read.csv("../data/fig3B_S7AB_linage.csv")

viruse_contig <- read.table("../data/viruse_contig")
meta <- read.table("../data/sample_meta")
############################## Data processing ##############################
Rcluster <- data.frame(matrix(nrow = 0,ncol=3))
for (i in 1:nrow(cluster)){
  if (length(grep("Clu",cluster[i,]))>0){
    cluster_name <- cluster[i,]
  }else{
    Rcluster <- rbind(Rcluster,data.frame(Cluster_name=cluster_name,Contig_number=cluster[i,]))
  }
}
Rcluster_beifen <- Rcluster
Rcluster <- Rcluster_beifen
Rcluster$Cluster_name <- sapply(strsplit(Rcluster$Cluster_name,split=">"),"[",2)
Rcluster$Cluster_number <- sapply(strsplit(Rcluster$Contig_number,split="\t"),"[",1)
colnames(Rcluster)[which(colnames(Rcluster)=="Contig_number")] <- "Rcluster_Content"
Rcluster$Contig_number <- sapply(strsplit(Rcluster$Rcluster_Content,split=">"),"[",2)
Rcluster$Similarity <- sapply(strsplit(Rcluster$Contig_number,split=" "),"[",3)
Rcluster$Contig_number_part <- sapply(strsplit(Rcluster$Contig_number,split=" "),"[",1)
Rcluster$Contig_number_part <- sapply(strsplit(Rcluster$Contig_number_part,split="\\."),"[",1)
Rcluster$Similarity[is.na(Rcluster$Similarity)] <- "*"


clstr_name_change$Contig_number_part <- substring(clstr_name_change$V2, 1, last = 19)
clstr_name_change$full_name <-  substring(clstr_name_change$V1, 1)
Rcluster <- merge(Rcluster,unique(clstr_name_change[,c(2:3)]),by="Contig_number_part")
colnames(Rcluster)[which(colnames(Rcluster)=="V2")] <- "full_name"
Rcluster$full_name <- sapply(strsplit(Rcluster$full_name,split=" "),"[",1)
Rcluster$full_name <- paste0("NODE_",Rcluster$full_name)
Rcluster <-Rcluster[is.na(Rcluster$full_name)==F,]
Rcluster3 <- merge(Rcluster,contigs_phacts[,c("Contig_number","Sample_id")],by.x="full_name",by.y="Contig_number")
Rcluster3 <-merge(Rcluster3,meta,by="Sample_id")


viruse_contig <- data.frame(viruse_contig)
Rcluster2 <- merge(Rcluster,viruse_contig,by.x="full_name",by.y="Contig_number")
Rcluster2$Bacteria_Host_G <- sapply(strsplit(Rcluster2$Bacteria_Host,split = ";g__"),"[",2)
Rcluster2$Bacteria_Host_G <- sapply(strsplit(Rcluster2$Bacteria_Host_G,split = ";s__"),"[",1)
Cluster_Linage<- aggregate(list(Bacteria_Host=Rcluster2$Bacteria_Host),by=list(Cluster_name=Rcluster2$Cluster_name),function(x) ifelse(length(table(x)[is.na(table(x))==F])>0,
                                                                                                                                       names(table(x)[order(table(x),decreasing = T)])[1],NA)  )   

Cluster_Linage$Lineage <- sapply(strsplit(Cluster_Linage$Bacteria_Host,split = ";g__"),"[",2)
Cluster_Linage$Lineage <- sapply(strsplit(Cluster_Linage$Lineage,split = ";s__"),"[",1)
Cluster_Linage$Lineage[is.na(Cluster_Linage$Lineage)] <- paste("unclassfied_genus_",1:length(Cluster_Linage$Lineage[is.na(Cluster_Linage$Lineage)] ),sep = "_")
Cluster_Linage$Lineage <- sapply(strsplit(Cluster_Linage$Lineage,split = ";"),"[",1)


Rcluster3<-merge(Rcluster3,reads_number2s,by.x="full_name",by.y="Contig_number")
Rcluster3$Contig_Abundance <- Rcluster3$`Mapped_Reads_Number`/Rcluster3$num_seqs*1000000
Rcluster3$length <- as.numeric(sapply(strsplit(Rcluster3$full_name,split = "_",fixed = T),"[",4))
Rcluster3$Contig_Abundance2 <- Rcluster3$Contig_Abundance/Rcluster3$length*300#对长度进行归一化
Rcluster3 <- merge(Rcluster3,Cluster_Linage,by="Cluster_name",all.x=T)
contigs_phacts$Phacts_result[-grep("Lytic|Temperate",contigs_phacts$Phacts_result)] <- contigs_phacts$Probability[-grep("Lytic|Temperate",contigs_phacts$Phacts_result)] 
Rcluster3 <- merge(Rcluster3,contigs_phacts,by.x="full_name",by.y="Contig_number")
a<- Rcluster3[grep("Crass",Rcluster3$Lineage,ignore.case = T),]


rongyuan_TT <- aggregate(list(Contig_Abundance=Rcluster3$Contig_Abundance2),
                         by=list(Group=Rcluster3$Group,
                                 PhageLifestyle=Rcluster3$Phacts_result,
                                 Viruse_Family=Rcluster3$Cluster_name,
                                 Sample_id=Rcluster3$Sample_id.x), sum)
rongyuan_TT <-merge(rongyuan_TT,meta[c("Sample_id","num_seqs")],by="Sample_id")
rongyuan_TT$Time <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",2)
rongyuan_TT$Diet <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",1)
rongyuan_TT <- rongyuan_TT[grep("D0|D10",rongyuan_TT$Group),]
rongyuan_TT <- rongyuan_TT[order(rongyuan_TT $Diet,rongyuan_TT$Time),]
rongyuan_TT$Group <- paste(rongyuan_TT$Diet,rongyuan_TT$Time,sep = "_")

rongyuan_TT <- spread(rongyuan_TT[colnames(rongyuan_TT)%in%c("num_seqs")==F],key =  'PhageLifestyle', value = 'Contig_Abundance')
rongyuan_TT[is.na(rongyuan_TT)==T] <- "0"
rongyuan_TT <- aggregate(list(Temperate=as.numeric(rongyuan_TT$Temperate),
                              Lytic=as.numeric(rongyuan_TT$Lytic)),
                         by=list(Group=rongyuan_TT$Group,Viruse_Family=rongyuan_TT$Viruse_Family), sum)
rongyuan_TT$Time <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",2)
rongyuan_TT$Diet <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",1)
rongyuan_TT$Diet2 <- factor(rongyuan_TT$Diet,levels = c("CON","HFD","FUC") ,labels =c("Con","HFD","HFD+Fuco") )
rongyuan_TT$Group2 <- paste(as.character(rongyuan_TT$Diet),rongyuan_TT$Time,sep="_")

value=0.1
rongyuan_TT$Group_clu<- paste(rongyuan_TT$Group2,rongyuan_TT$Viruse_Family,sep="_")
rongyuan_TT$Temperate_chu_Lytic <- (as.numeric(rongyuan_TT$Lytic)+value)/(as.numeric(rongyuan_TT$Temperate)+value)
rongyuan_TT$Temperate_chu_Lytic <-abs(log(rongyuan_TT$Temperate_chu_Lytic,base = 2))
rongyuan_TT2 <- rongyuan_TT
rongyuan_TT <- spread(rongyuan_TT[,c("Time","Temperate_chu_Lytic","Viruse_Family","Diet")], key= 'Time', value =  'Temperate_chu_Lytic')
rongyuan_TT$D10_D0 <-(as.numeric(rongyuan_TT$D10)+0.1)/(as.numeric(rongyuan_TT$D0)+0.1)
rongyuan_TT$D10_D0 <-log10(rongyuan_TT$D10_D0+1)
rongyuan_TT <- spread(rongyuan_TT[,colnames(rongyuan_TT)%in%c("D0","D10")==F], key= 'Diet', value =  'D10_D0')
i=1

colnames(rongyuan_TT)[which(colnames(rongyuan_TT)=="Viruse_Family")] <- "Cluster_name"
rongyuan_TT <- merge(rongyuan_TT,Cluster_Linage,by="Cluster_name")
rongyuan_TT$Viruse_Family <- sapply(strsplit(rongyuan_TT$Bacteria_Host,split = ";"),"[",2)
rongyuan_TT$Viruse_Family <- sapply(strsplit(rongyuan_TT$Viruse_Family,split = ";"),"[",1)


rongyuan_TT$Viruse_Family[rongyuan_TT$Viruse_Family%in%names(table(rongyuan_TT$Viruse_Family)[table(rongyuan_TT$Viruse_Family)<15])] <- "Others"
rongyuan_TT <- rongyuan_TT[is.na(rongyuan_TT$Viruse_Family)==F,]
Cluster_Abundance <- aggregate(list(Contig_Abundance=Rcluster3$Contig_Abundance2),by=list(Viruse_Family=Rcluster3$Cluster_name,Sample_id=Rcluster3$Sample_id.x), sum)
Cluster_Abundance <- aggregate(list(Contig_Abundance=Cluster_Abundance$Contig_Abundance),by=list(Viruse_Family=Cluster_Abundance$Viruse_Family), median)
rongyuan_TT <-merge(rongyuan_TT,Cluster_Abundance,by.x="Cluster_name",by.y="Viruse_Family")

wilcox.test(rongyuan_TT$HFD,rongyuan_TT$FUC)
wilcox.test(rongyuan_TT$HFD,rongyuan_TT$CON)
wilcox.test(rongyuan_TT$CON,rongyuan_TT$FUC)
median(rongyuan_TT$HFD)
median(rongyuan_TT$FUC)
median(rongyuan_TT$CON)

rongyuan_TT5 <- data.frame(matrix(ncol=2,nrow=0))
for (i in 2:4){
  Group=colnames(rongyuan_TT)[i]
  rongyuan_TT4 <- data.frame(Group=Group,D10_D0=rongyuan_TT[,i],
                             Viruse_Family=rongyuan_TT[,c("Viruse_Family")],
                             Bacteria_Host=rongyuan_TT[,c("Bacteria_Host")] )
  rongyuan_TT5 <- rbind(rongyuan_TT5,rongyuan_TT4)
}


rongyuan_TT$Viruse_Family2 <- rongyuan_TT$Viruse_Family
rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"] <- substring(rongyuan_TT$Viruse_Family[rongyuan_TT$Viruse_Family!="Others"],4,50)
rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"] <- sapply(strsplit(rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"],split = "_A",fixed = T),"[",1)
rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"] <- sapply(strsplit(rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"],split = "_I",fixed = T),"[",1)
rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"] <- sapply(strsplit(rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"],split = "_C",fixed = T),"[",1)
# rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"] <- paste0("p__",rongyuan_TT$Viruse_Family2[rongyuan_TT$Viruse_Family!="Others"])
rongyuan_TT$Viruse_Family2 <- paste0(rongyuan_TT$Viruse_Family2,"_phage")
color <- c('gray',brewer.pal(8,'Set3')[1],brewer.pal(8,'Set1')[1],brewer.pal(11,'Set3')[6:8])
names(color) <- sort(unique(rongyuan_TT$Viruse_Family2))
max_for_color<- names(table(rongyuan_TT$Viruse_Family2)[order(table(rongyuan_TT$Viruse_Family2),decreasing = T)])[1]
names(color) <- c(max_for_color,unique(rongyuan_TT$Viruse_Family2)[unique(rongyuan_TT$Viruse_Family2)!=max_for_color])
color2 <- color


############################## [Fig3-B] Painting #####################################
rongyuan_TT <- rongyuan_TT[rowSums(is.na(rongyuan_TT[,2:4])==F)==3,]
ggtern(data = rongyuan_TT,aes(x=CON, y=HFD, z=FUC  ) )+
  theme_rgbw()+scale_color_manual(values =color )+
  stat_density_tern(geom='polygon',bdl = 0.04,alpha = 0.2,aes(fill=after_stat(level)))+
  scale_fill_gradient(low = brewer.pal(8,'YlOrRd')[1],high = brewer.pal(8,'YlOrRd')[5])+
  geom_point(aes(color=Viruse_Family2,size=Contig_Abundance),shape=21)
ggsave("../fig3-B.pdf")


############################## [FigS7-A] Painting #####################################
rongyuan_TT5$Group <- factor(rongyuan_TT5$Group,levels =c("CON","HFD","FUC") ,labels = c("CON","HFD","FUC"))
ggplot(rongyuan_TT5,aes(x=Group,y=log10(D10_D0),fill=Group))+
  geom_boxplot()+
  theme_light()+ylab("Changes of Ratio of lytic phages to lysogenic phages")+
  # scale_fill_manual(values = c(brewer.pal(3,'Greens')[3],brewer.pal(3,'Reds')[3],brewer.pal(3,'Blues')[3]))
  scale_fill_manual(values =   c("#8BBCD6","#DD7F88","#A0CC58"))+
  geom_signif(comparisons = list(c("CON","HFD"),c("HFD","FUC"),c("CON","FUC")),color="black",
              # y_position = c(2.5,2.8,3.2))
              y_position = c(1,1.4,1.8))
wilcox.test(rongyuan_TT5$D10_D0[rongyuan_TT5$Group=="FUC"],rongyuan_TT5$D10_D0[rongyuan_TT5$Group=="CON"])
ggsave("../figS7-A.pdf")


############################## [FigS7-B] Painting #####################################
rongyuan_TT6<- rongyuan_TT5[grep("Bacteroidota",rongyuan_TT5$Bacteria_Host),]
rongyuan_TT6 <- rongyuan_TT6[is.na(rongyuan_TT6$D10_D0)==F,]
rongyuan_TT6$D10_D0 <- log10(rongyuan_TT6$D10_D0)
rongyuan_TT6$Group <- factor(rongyuan_TT6$Group,levels =c("CON","HFD","FUC") ,labels = c("CON","HFD","FUC"))
ggplot(rongyuan_TT6,aes(x=Group,y=D10_D0,fill=Group))+
  geom_boxplot()+
  theme_light()+ylab("Changes of Ratio of lytic phages to lysogenic phages")+
  # scale_fill_manual(values = c(brewer.pal(3,'Greens')[3],brewer.pal(3,'Reds')[3],brewer.pal(3,'Blues')[3]))
  scale_fill_manual(values =   c("#8BBCD6","#DD7F88","#A0CC58"))+
  geom_signif(comparisons = list(c("CON","HFD"),c("HFD","FUC"),c("CON","FUC")),color="black",y_position = c(2.5,2.8,3.2))+
  facet_wrap(.~Viruse_Family,)
ggsave("../figS7-B.pdf")

