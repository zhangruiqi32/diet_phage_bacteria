##Script to Figure 2-DE and Supplementary fig6-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","ggridges","viridisLite","viridis","circlize","RColorBrewer","ggpubr","ggplot2")
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
data_g_B <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2DE_S6C_gB.csv", row.names = 1)
data_g_V <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2DE_S6C_gV.csv", row.names = 1)

df_BV_name_duiyings <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2DE_S6C_nameBV.csv")
meta <- read.csv("~/cooperation/202409zhaofengxiang/data/fig2_sample.csv")
x <- read.csv("~/cooperation/202409zhaofengxiang/data/DABs_lefse.csv",sep=",",header = T )

############################## Correlatin analysis ##############################
BV_p_value<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],]),
           t(data_g_V[x[3],]))[["p.value"]]})

BV_R2<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],]),
           t(data_g_V[x[3],]))[["estimate"]]})

BV_p_value_Control<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="CON"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="CON"]]))[["p.value"]]})

BV_R2_Control<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="CON"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="CON"]]))[["estimate"]]})

BV_p_value_Fucoidan<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="FUC"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="FUC"]]))[["p.value"]]})

BV_R2_Fucoidan<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="FUC"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="FUC"]]))[["estimate"]]})

BV_p_value_HFD<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="HFD"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="HFD"]]))[["p.value"]]})

BV_R2_HFD<- apply(df_BV_name_duiyings,1,function(x){
  cor.test(t(data_g_B[x[1],colnames(data_g_B)%in%meta$sample[meta$treatment=="HFD"]]),
           t(data_g_V[x[3],colnames(data_g_V)%in%meta$sample[meta$treatment=="HFD"]]))[["estimate"]]})


df_BV_name_xiangguanxings_serach <- data.frame(df_BV_name_duiyings,BV_p_value=BV_p_value,BV_R2=BV_R2,
                                               BV_p_value_Control=BV_p_value_Control,BV_R2_Control,
                                               BV_p_value_Fucoidan=BV_p_value_Fucoidan,BV_R2_Fucoidan,
                                               BV_p_value_HFD,BV_R2_HFD)


df_BV_name_xiangguanxings_serach$Group_forp <-apply(df_BV_name_xiangguanxings_serach[,c("Bacteria_name","Viruses_name")],1,function(x) 
  paste(  sapply(strsplit(x[1],split = "\\|"),"[",length(strsplit(x[1],split = "\\|")[[1]])),
          sapply(strsplit(x[2],split = "\\|"),"[",length(strsplit(x[2],split = "\\|")[[1]])),sep = "_"))
df_BV_name_xiangguanxings_chayi <- df_BV_name_xiangguanxings_serach[df_BV_name_xiangguanxings_serach$Bacteria_name%in%x$Spe,]

#相关性峰峦图
Groups <- colnames(df_BV_name_xiangguanxings_serach)
Groups[which(Groups=="BV_p_value")] <- "BV_p_value_ALL"
ii=5

df_BV_name_xiangguanxings_serach_fengluans <- data.frame(matrix(ncol=7,nrow=0))
for (ii in grep("p_value",colnames(df_BV_name_xiangguanxings_serach))){
  df_BV_name_xiangguanxings_serach_fengluan <- data.frame(df_BV_name_xiangguanxings_serach[,-grep("p_value|R2",colnames(df_BV_name_xiangguanxings_serach))],
                                                          p_value= df_BV_name_xiangguanxings_serach[,ii],
                                                          R2=df_BV_name_xiangguanxings_serach[,ii+1],
                                                          Group=unlist(strsplit(Groups[ii],split = "value_"))[2])
  df_BV_name_xiangguanxings_serach_fengluans <- rbind(df_BV_name_xiangguanxings_serach_fengluans,df_BV_name_xiangguanxings_serach_fengluan)
}
df_BV_name_xiangguanxings_serach_fengluans <- df_BV_name_xiangguanxings_serach_fengluans[df_BV_name_xiangguanxings_serach_fengluans$Group%in%c("ALL","Disturbance","raodong")==F ,]
unique(df_BV_name_xiangguanxings_serach_fengluan$Viruses_name)
df_BV_name_xiangguanxings_serach_fengluans_quchu1<- substr(df_BV_name_xiangguanxings_serach_fengluans$Group_forp,1,nchar(df_BV_name_xiangguanxings_serach_fengluans$Group_forp)-1)
df_BV_name_xiangguanxings_serach_fengluans <-df_BV_name_xiangguanxings_serach_fengluans[df_BV_name_xiangguanxings_serach_fengluans_quchu1%in%df_BV_name_xiangguanxings_serach_fengluans$Group_forp==F,]


############################## [Fig2-D] Processing and Painting #####################################
#p值小于0.05
df_BV_name_xiangguanxings_serach_fengluans_forp <- df_BV_name_xiangguanxings_serach_fengluans[df_BV_name_xiangguanxings_serach_fengluans$p_value<0.05&is.na(df_BV_name_xiangguanxings_serach_fengluans$p_value)==F,]
title="#All bacteria are associated with their viruses(p<0.05)"

df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Bacteria_name,split = "\\."),"[",2)
df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank )
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Viruses_name,split = "\\."),"[",6)
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank)
df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <- as.character(df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank)


df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp[df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank%in%c("p__Cyanobacteria","p__Fusobacteria")==F,]
df_BV_name_xiangguanxings_serach_fengluans_forp2<- data.frame(table(df_BV_name_xiangguanxings_serach_fengluans_forp2[,c("Bact_p_rank","Group")]))
df_BV_name_xiangguanxings_serach_fengluans_forp2<- spread(df_BV_name_xiangguanxings_serach_fengluans_forp2,key ="Group",value ="Freq")
rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2) <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,1]
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,-1]
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,c("Control","HFD","Fucoidan")]
# colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)[which(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)=="Fucoidan")] <- "HFD+Fucoidan"
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[order(rowSums(df_BV_name_xiangguanxings_serach_fengluans_forp2),decreasing = T),]
color_name=c(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2),rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2))
color=c("#8BBCD6","#DD7F88","#A0CC58","#85C3B9","#7CA5C2","#E9C4DA","#DDEAC4")
names(color) <- color_name

chordDiagram(as.matrix(df_BV_name_xiangguanxings_serach_fengluans_forp2 ),
             grid.col = color,transparency = 0.5,grid.border = NA,
             # link.border = 'gray', transparency = 0.5,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.04)
)
legend("bottomright", 
       legend=as.character(rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2)),
       fill = color[4:7], border = F,bty="n",
       title = "")
circos.clear()
ggsave("../fig2-D.pdf")


############################## [FigS6-C] Processing and Painting #####################################
#差异菌,p<0.05
x$Spe2 <- gsub("\\|","\\.",x$Spe)
df_BV_name_xiangguanxings_serach_fengluans_forp <- df_BV_name_xiangguanxings_serach_fengluans[df_BV_name_xiangguanxings_serach_fengluans$p_value<0.05&is.na(df_BV_name_xiangguanxings_serach_fengluans$p_value)==F&df_BV_name_xiangguanxings_serach_fengluans$Bacteria_name%in%x$Spe2,]
title="#Differential bacteria are associated with their viruses(p<0.05)"

df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Bacteria_name,split = "\\."),"[",2)
df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank )
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Viruses_name,split = "\\."),"[",6)
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank)
df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <- as.character(df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank)


df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp[df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank%in%c("p__Cyanobacteria","p__Fusobacteria")==F,]
df_BV_name_xiangguanxings_serach_fengluans_forp2<- data.frame(table(df_BV_name_xiangguanxings_serach_fengluans_forp2[,c("Bact_p_rank","Group")]))
df_BV_name_xiangguanxings_serach_fengluans_forp2<- spread(df_BV_name_xiangguanxings_serach_fengluans_forp2,key ="Group",value ="Freq")
rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2) <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,1]
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,-1]
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[,c("Control","HFD","Fucoidan")]
# colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)[which(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)=="Fucoidan")] <- "HFD+Fucoidan"
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[order(rowSums(df_BV_name_xiangguanxings_serach_fengluans_forp2),decreasing = T),]
color_name=c(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2),rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2))
color=c("#8BBCD6","#DD7F88","#A0CC58","#85C3B9","#DDEAC4","#7CA5C2","#E9C4DA")
names(color) <- color_name


chordDiagram(as.matrix(df_BV_name_xiangguanxings_serach_fengluans_forp2 ),
             grid.col = color,transparency = 0.5,grid.border = NA,
             # link.border = 'gray', transparency = 0.5,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.04)
)
legend("bottomright", 
       legend=as.character(rownames(df_BV_name_xiangguanxings_serach_fengluans_forp2)),
       fill = color[4:7], border = F,bty="n",
       title = "")
circos.clear()
ggsave("../figS6-C.pdf")


############################## [Fig2-E] Processing and Painting #####################################
#差异菌,p<0.05
x$Spe2 <- gsub("\\|","\\.",x$Spe)
df_BV_name_xiangguanxings_serach_fengluans_forp <- df_BV_name_xiangguanxings_serach_fengluans[df_BV_name_xiangguanxings_serach_fengluans$p_value<0.05&is.na(df_BV_name_xiangguanxings_serach_fengluans$p_value)==F&df_BV_name_xiangguanxings_serach_fengluans$Bacteria_name%in%x$Spe2,]
title="#Differential bacteria are associated with their viruses(p<0.05)"

df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Bacteria_name,split = "\\."),"[",2)
df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Bact_p_rank )
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <- sapply(strsplit(df_BV_name_xiangguanxings_serach_fengluans_forp$Viruses_name,split = "\\."),"[",6)
df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank <-factor(df_BV_name_xiangguanxings_serach_fengluans_forp$Viru_p_rank)

df_BV_name_xiangguanxings_serach_fengluans_forp2<- df_BV_name_xiangguanxings_serach_fengluans_forp
colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)[which(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)=="Bact_p_rank")] <-  "Bacteria"
df_BV_name_xiangguanxings_serach_fengluans_forp2$Group <- factor(df_BV_name_xiangguanxings_serach_fengluans_forp2$Group,
                                                                 levels =  c("Control","HFD","Fucoidan"),labels = c("CON","HFD","FUC") )
df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[order(df_BV_name_xiangguanxings_serach_fengluans_forp2$Group,-df_BV_name_xiangguanxings_serach_fengluans_forp2$R2),]
# df_BV_name_xiangguanxings_serach_fengluans_forp2 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[order(df_BV_name_xiangguanxings_serach_fengluans_forp2$R2,decreasing = T),]
df_BV_name_xiangguanxings_serach_fengluans_forp2$Group_forp <- factor(df_BV_name_xiangguanxings_serach_fengluans_forp2$Group_forp,levels = unique(df_BV_name_xiangguanxings_serach_fengluans_forp2$Group_forp))
df_BV_name_xiangguanxings_serach_fengluans_forp2$Viral_abundance<- apply(df_BV_name_xiangguanxings_serach_fengluans_forp2,1,function(x) mean(t(data_g_V[x[which(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)=="Viruses_name")],
                                                                                                                                                        meta$Sample_id[(meta$GD=="D0"|meta$GD=="D10")&
                                                                                                                                                                     meta$Diet==x[which(colnames(df_BV_name_xiangguanxings_serach_fengluans_forp2)=="Group")]]])))
color <- c("#8BBCD6","#DD7F88","#A0CC58")
names(color) <- c("CON","HFD","FUC")
p1<- ggplot(df_BV_name_xiangguanxings_serach_fengluans_forp2,aes(x=Group_forp,y=R2,colour=Group,size=Viral_abundance))+
  geom_point()+theme_bw()+
  theme(axis.text.x = element_blank(),axis.ticks=element_blank())+
  scale_size_continuous(range = c(0.8,2))+
  scale_color_manual(values = color)+theme(legend.position="bottom")+xlab("#The correlations between Viruses and Bacteria")
df_BV_name_xiangguanxings_serach_fengluans_forp3 <- df_BV_name_xiangguanxings_serach_fengluans_forp2[!duplicated(df_BV_name_xiangguanxings_serach_fengluans_forp2$Group_forp),]
df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2 <- sapply(strsplit(as.character(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp),split = "s__"),"[",2)
df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses <- sapply(strsplit(as.character(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2),split = "_"),"[",1)
df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses <-as.character(factor(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses,levels = c("Acinetobacter","Bacteroides","crAssphage","Faecalibacterium","Gordonia"),
                                                                               labels =c("Acinetobacter_virus","Bacteroides_phage","crAssphage","Faecalibacterium_virus","Gordonia_phage") ))
df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family <- sapply(strsplit(as.character(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses_name),split = "\\."),"[",6)
df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family <- sapply(strsplit(as.character(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family),split = "__"),"[",2)
df_BV_name_xiangguanxings_serach_fengluans_forp3 <- df_BV_name_xiangguanxings_serach_fengluans_forp3[order(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses),]
annotation_colors2 <- c("#66C2A5","#FC8D62",brewer.pal(7,'Blues')[1:4],"#ABD9E9","#8DA0CB",brewer.pal(7,'RdPu')[1:4],"#E78AC3",brewer.pal(7,'YlGn')[2:4])
names(annotation_colors2) <- unique(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2)

aggregate(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2,by=list(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2),function(x)length(unique(x)))

df_BV_name_xiangguanxings_serach_fengluans_forp4 <- pivot_longer(df_BV_name_xiangguanxings_serach_fengluans_forp3,names_to="Label_Title",values_to ="Label_Content",cols = c("Bacteria","Viruses","Viruse_Family"))
df_BV_name_xiangguanxings_serach_fengluans_forp4 <- df_BV_name_xiangguanxings_serach_fengluans_forp4[order(df_BV_name_xiangguanxings_serach_fengluans_forp4$Label_Title,df_BV_name_xiangguanxings_serach_fengluans_forp4$Group_forp),]
table(df_BV_name_xiangguanxings_serach_fengluans_forp3[c("Bacteria","Viruses")])
df_BV_name_xiangguanxings_serach_fengluans_forp3$add <- paste(df_BV_name_xiangguanxings_serach_fengluans_forp3$Bacteria,df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses,df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family,sep="_")
p2 <- ggplot(df_BV_name_xiangguanxings_serach_fengluans_forp4,aes(x=Label_Title,y=Label_Content,fill=Label_Content))+geom_bar(stat="identity",position = "stack",width=1)+
  # scale_fill_manual(values = annotation_colors2)+theme_bw()+
  # scale_color_manual(values = annotation_colors2)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    panel.border = element_blank (),
    legend.position="bottom"
  )+xlab("")+ylab("")
annotation_colors <- brewer.pal(8,'Set2')[c(2,4,6:8)]
names(annotation_colors) <- unique(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruses)
p3<- ggplot(df_BV_name_xiangguanxings_serach_fengluans_forp3,aes(x=Group_forp,y=1,fill=Viruses))+geom_bar(stat="identity",position="dodge",width=1)+
  scale_fill_manual(values = annotation_colors)+theme_bw()+
  scale_color_manual(values = annotation_colors)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    panel.border = element_blank (),
    legend.position="bottom"
  )+xlab("")+ylab("")
annotation_colors3 <- brewer.pal(9,'Set3')[6:9]

names(annotation_colors3) <- as.character(unique(df_BV_name_xiangguanxings_serach_fengluans_forp3$Bacteria))
p4<- ggplot(df_BV_name_xiangguanxings_serach_fengluans_forp3,aes(x=Group_forp,y=1,fill=Bacteria))+geom_bar(stat="identity",position="dodge",width=1)+
  scale_fill_manual(values = annotation_colors3)+theme_bw()+
  scale_color_manual(values = annotation_colors3)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    panel.border = element_blank (),
    legend.position="bottom"
  )+xlab("")+ylab("")
annotation_colors4 <- brewer.pal(6,'Set2')[c(1,3,5)]
brewer.pal(11,'RdYlBu')
names(annotation_colors4) <- unique(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family)
aggregate(df_BV_name_xiangguanxings_serach_fengluans_forp3$Group_forp2,by=list(df_BV_name_xiangguanxings_serach_fengluans_forp3$Viruse_Family),function(x)length(unique(x)))
p5<- ggplot(df_BV_name_xiangguanxings_serach_fengluans_forp3)+geom_bar(aes(x=Group_forp,y=1,fill=Viruse_Family),stat="identity",position="dodge",width=1)+
  scale_fill_manual(values = annotation_colors4)+theme_bw()+
  scale_color_manual(values = annotation_colors4)+
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text=element_blank(),
    panel.border = element_blank (),
    legend.position="bottom"
  )+xlab("")+ylab("")

ggarrange(p4,p5,p3,p1,ncol = 1,heights = c(0.3,0.3,0.3,1) )

ggsave("../fig2-E.pdf")

