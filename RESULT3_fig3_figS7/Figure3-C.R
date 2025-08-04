##Script to Figure 3-C.

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
viruse_contig2 <- read.table("../data/viruse_contig")
meta <- read.table("../data/sample_meta")


############################## Data processing ##############################
rongyuan_TT <- aggregate(list(Contig_Abundance=viruse_contig2$Contig_Abundance),
                         by=list(Group=viruse_contig2$Group,
                                 PhageLifestyle=viruse_contig2$bac_result,
                                 Sample_id=viruse_contig2$Sample_id),
                         sum)

rongyuan_TT <- merge(rongyuan_TT,meta[c("Sample_id","num_seqs")],by="Sample_id")
rongyuan_TT <- rongyuan_TT[grep("D0|D10",rongyuan_TT$Group),]
rongyuan_TT <- spread(rongyuan_TT, key= 'PhageLifestyle', value =  'Contig_Abundance')
rongyuan_TT[is.na(rongyuan_TT)==T] <- "0"
rongyuan_TT$Temperate_chu_Lytic <- as.numeric(rongyuan_TT$Lytic)/(as.numeric(rongyuan_TT$Temperate)+as.numeric(rongyuan_TT$Lytic))
rongyuan_TT$Temperate_chu_Lytic <- as.numeric(rongyuan_TT$Lytic)/(as.numeric(rongyuan_TT$Temperate))
rongyuan_TT$Time <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",2)
rongyuan_TT$Diet <- sapply(strsplit(as.character(rongyuan_TT$Group),split = "_"),"[",1)
rongyuan_TT$Group <- paste(rongyuan_TT$Diet,rongyuan_TT$Time,sep="_")
rongyuan_TT$Temperate_chu_Lytic <- log(rongyuan_TT$Temperate_chu_Lytic,base =2 )

outlines_test <- function(x){
  Q1 <- summary(x)[2]
  Q3 <- summary(x)[5]
  IQR <- Q3-Q1
  min <- Q1-1.2*IQR
  max <- Q3+1.2*IQR
  x<max&x>min
}
i=1
rongyuan_TT$Group2 <- rongyuan_TT$Group
rongyuan_TT2 <- data.frame(matrix(nrow = 0,ncol=10))
for (i in 1:length(unique(rongyuan_TT$Group2))){
  Group <- unique(rongyuan_TT$Group2)[i]
  rongyuan_TT2_group <-rongyuan_TT[rongyuan_TT$Group2==Group,]
  rongyuan_TT2_group <- rongyuan_TT2_group[outlines_test(rongyuan_TT2_group$Temperate_chu_Lytic[rongyuan_TT2_group$Group2==Group]),]
  rongyuan_TT2 <-rbind(rongyuan_TT2,rongyuan_TT2_group)
}
rongyuan_TT<- rongyuan_TT2
aggregate(1:nrow(rongyuan_TT),by=list(rongyuan_TT$Diet),function(x) wilcox.test(rongyuan_TT$Temperate_chu_Lytic[x]~rongyuan_TT$Group[x])[['p.value']])
rongyuan_TT$Diet <- factor(rongyuan_TT$Diet,levels = c("CON","HFD","FUC"),labels = c("CON","HFD","FUC")  )
rongyuan_TT$Group <- paste(rongyuan_TT$Diet,rongyuan_TT$Time,sep="_")
rongyuan_TT$Group <- factor(rongyuan_TT$Group,levels = c("CON_D0", "CON_D10","HFD_D0", "HFD_D10","FUC_D0", "FUC_D10"),
                            labels = c("CON_D0", "CON_D10","HFD_D0", "HFD_D10","FUC_D0", "FUC_D10")  )

############################## [Fig3-C] Painting #####################################
color_Group <- c("#DDDDDD","#ACD26A",
                 "#DDDDDD","#DD7F88",
                 "#DDDDDD","#8BBCD6")
names(color_Group)<- unique(rongyuan_TT$Group)

ggplot(rongyuan_TT ,aes(x=Group,y=Temperate_chu_Lytic,fill=Group))+
  geom_boxplot()+
  theme_light()+ylab("Ratio of lytic phages to lysogenic phages")+
  scale_fill_manual(values =color_Group  )+
  xlab("")+
  # facet_wrap(.~chayi,scales = "free_y")+
  geom_signif(comparisons = list(c("CON_D0","CON_D10"),
                                 c("HFD_D0","HFD_D10"),
                                 c("FUC_D0","FUC_D10")),
              y_position =3.6,
              map_signif_level = TRUE, #转换p值为星号
              #xmin = c(0.8, 1.8, 2.8),
              #xmax = c(1.2, 2.2, 3.18),
              # ,tip_length = c(0.03,0.45,0.6,0.05,0.55,0.35)
              tip_length = 0.03,
              color="black")
ggsave("../fig3-C.pdf")


