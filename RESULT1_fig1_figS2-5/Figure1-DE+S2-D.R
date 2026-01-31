##Script to Figure 1-DE and Supplementary fig2-D.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","treemap","treemapify","ggalluvial","RColorBrewer","ggplot2")
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
Diet1 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1DE_S2D_1.csv", row.names = 1, header = T)
As<- read.table("~/cooperation/202409zhaofengxiang/data/fig1_As")
diet1 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1DE_S2D_diet.csv")

data_g_V <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_V")
data_g_B <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_B")
prj2 <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_prj_data")

############################## [Fig1-D] Processing and Painting #####################################
x <- Diet1
x <- x[rowSums(is.na(x)==F)>0,]
x[is.na(x)] <- 0
x2 <- x
x2<- x2[rowSums(x2!=0)>=2,]

x4 <- data.frame(species=rownames(data_g_V))
x4$species_host <- sapply(strsplit(x4$species,split = "s__",fixed = T),"[",2)
x4$species_host[grep("crass",x4$species_host ,ignore.case = T)]<- "crAssphage"
x4$species_host <- sapply(strsplit(x4$species_host,split = "_",fixed = T),"[",1)
x4<- data.frame(table(x4$species_host))
colnames(x4) <- c("Type","Number_all")

x2 <- data.frame(species=rownames(x2))
x2$species_host <- sapply(strsplit(x2$species,split = "s__",fixed = T),"[",2)
x2$species_host[grep("crass",x2$species_host ,ignore.case = T)]<- "crAssphage"
x2$species_host <- sapply(strsplit(x2$species_host,split = "_",fixed = T),"[",1)
x2<- data.frame(table(x2$species_host))
colnames(x2) <- c("Type","Number_diff")

x2 <- merge(x2,x4,by="Type")
x2 <- x2[x2$Number_all>=5,]
x2 <- x2[x2$Number_diff>=2,]
x2$percent <- x2$Number_diff/x2$Number_all*100

a <- data.frame(Bact=rownames(data_g_B))
a$phylum <- sapply(strsplit(a$Bact,split = "p__") ,"[",2)
a$phylum <- sapply(strsplit(a$phylum,split = "|",fixed = T) ,"[",1)
a$genus <- sapply(strsplit(a$Bact,split = "g__") ,"[",2)
a$genus <- sapply(strsplit(a$genus,split = "|",fixed = T) ,"[",1)
x2<- merge(x2,unique(a[2:3]),by.x = "Type",by.y="genus",all.x = T)
x2$phylum[x2$Type=="crAssphage"] <- "Bacteroidetes"
x2[is.na(x2)] <- "Others"
colnames(x2)[colnames(x2)=="Type"] <- "Host"

a <- data.frame(table(x2$phylum))
a <- a[order(a$Freq,decreasing = T),]
x2$phylum <- factor(x2$phylum,levels = a$Var1,labels = a$Var1)
x2 <- x2[order(x2$phylum,-x2$percent),]

#### Painting
treemap(x2,index = colnames(x2)[c(5,1)],vSize="percent" ,palette = "Set2",border.lwds = 0.8,border.col = "white")

ggplot(df, aes(area = percent, fill = phylum, label = Host,
               subgroup = phylum)) + #子分组映射
  scale_fill_manual(values=c("#7CA5C2","#85C3B9","#DDEAC4", "#B7B3D1","#E9C4DA"))+ #Acti\Bacter\Firmi\Other\Prote
  geom_treemap(alpha=0.5) +
  geom_treemap_subgroup_border(color = 'white',size = 4) + #子分组边界线绘制
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.3, colour = "black", fontface = "italic", min.size = 0) +
  geom_treemap_text(colour = "black", place = "topleft",fontface = "bold",size = 6)+
  theme(legend.position = "none")
ggsave("../fig1-D.pdf") 


############################## [figS2-D] Processing and Painting #####################################
x <- Diet1
x <- x[rowSums(is.na(x)==F)>0,]
x[is.na(x)] <- 0
x2 <- x
x2<- x2[grepl("s_",rownames(x2)),]
x2<- x2[rowSums(x2!=0)>=2,]
x2<- data.frame(diff=colSums(x2!=0),spes=colnames(x2))

x4<- data.frame(all=colSums(is.na(As[2:ncol(As)])==F) ,spes=colnames(As)[2:ncol(As)] )
fenleis<-c("Run","BioProject","Diet1","Diet","Nutrition1","Nutrition","Antibiotic","Continent","Host","Age","Age2","Sex" )
i=3
prj_index <- i
x4$spes <- paste0(x4$spes,".",fenleis[prj_index])
x2$spes <- gsub("Abundance.","",x2$spes)
x2$spes
x2 <- merge(x2,x4,by="spes")
x2$percent <- x2$diff/x2$all*100
x2 <- x2[x2$percent>=5,]

x[is.na(x)==T] <- 0
i <- prj_index
diets<- unique(prj2[,fenleis[i]])
diets_index <- 1
x4 <- data.frame(matrix(ncol = 1,nrow=0))
for (diets_index in 1:length(diets)){
  diet <- diets[diets_index]
  prjs_clean_x <-paste0(unique(prj2[prj2[,fenleis[i]]==diet,"BioProject"]),".",fenleis[i])
  x4 <- rbind(x4,data.frame(diet=diet,
                            all=length(prjs_clean_x ),
                            diff=nrow(x2[x2$spes%in%prjs_clean_x,  ]),
                            percent=nrow(x2[x2$spes%in%prjs_clean_x,  ])/length(prjs_clean_x ))  )
  
}
x4 <- x4[x4$all>=5,]
x4 <- merge(x4,diet1,by.x = "diet",by.y = "Diet",all=T)

x4 <- x4[is.na(x4$Type)==F,]
x4$diet <- gsub("_"," ",x4$diet)
 
#### Painting
treemap(x4,index = colnames(x4)[c(5,1)],vSize="percent" ,
        # palette = "Set2",
        border.lwds = 0.8,border.col = "white")

x4$Type <- as.factor(x4$Type)
ggplot(x4, aes(area = percent, fill = Type, label = diet,
               subgroup = Type)) + #子分组映射
  scale_fill_manual(values=c("#F3A59A",NA,"#80D0C3",NA,"#9EAAC4","#A6DDEA",NA))+ 
  geom_treemap(alpha=0.5) +
  geom_treemap_subgroup_border(color = 'white',size = 4) + #子分组边界线绘制
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.3, colour = "black", fontface = "italic", min.size = 0) +
  geom_treemap_text(colour = "black", place = "topleft",fontface = "bold",size = 6)+
  theme(legend.position = "none")
ggsave("../figS2-D.pdf") 


############################## [Fig1-E] Processing and Painting #####################################
x <- Diet1
x <- x[rowSums(is.na(x)==F)>0,]
x[is.na(x)] <- 0
i <- prj_index
diets<- unique(prj2[,fenleis[i]])
diets_index <- 1
colnames(x) <- gsub("Abundance.","",colnames(x) )

spe_x_prj_diets_tongji <-  data.frame(matrix(ncol = 1,nrow=0))
colnames(spe_x_prj_diets_tongji) <-"spe"
spe_x_prj_diets_zong<-  data.frame(matrix(ncol = 1,nrow=0))
colnames(spe_x_prj_diets_zong) <-"spe"

for (diets_index in 1:length(diets)){
  diet <- diets[diets_index]
  prjs_clean_x <-paste0(unique(prj2[prj2[,fenleis[i]]==diet,"BioProject"]),".",fenleis[i])
  prjs_clean <-unique(prj2[prj2[,fenleis[i]]==diet,"BioProject"])
  prjs_reduce_x <-paste0(prjs_clean,"_",fenleis[i])
  
  x_prj_diet <- x[,intersect(prjs_clean_x,colnames(x)  )]
  if (class(x_prj_diet)=="data.frame"){
    x2 <- x_prj_diet[rowSums(x_prj_diet!=0)/ncol(x_prj_diet)>=0.1,]
    # x2 <- x_prj_diet[rowSums(x_prj_diet!=0)>0,]
    x2 <- data.frame(species=rownames(x2))
  }else{
    x2 <- data.frame(species=rownames(x)[x_prj_diet!=0]  )
  }
  
  if(nrow(x2)!=0){
    x2$species_host <- sapply(strsplit(x2$species,split = "s__",fixed = T),"[",2)
    x2$species_host[grep("crass",x2$species_host ,ignore.case = T)]<- "crAssphage"
    x2$species_host <- sapply(strsplit(x2$species_host,split = "_",fixed = T),"[",1)
    x2<- data.frame(table(x2$species_host))
    colnames(x2) <- c("spe",diet)
    spe_x_prj_diets_tongji <- merge(spe_x_prj_diets_tongji,x2,by="spe",all = T)
  }
  
  x12<- As
  x12[is.na(x12)] <- 0
  x12 <- x12[rowSums(x12!=0)>0,colSums(x12!=0)>0]######筛选标准1（存在比例标准）,As筛选标准2（是否存在标准）
  rownames(x12) <- x12$spes
  x12 <- x12[,-1]
  colnames(x12) <-paste0(colnames(x12),".",fenleis[i])  
  x_prj_diet <- x12[,intersect(prjs_clean_x,colnames(x12))]
  if (class(x_prj_diet)=="data.frame"){
    x2 <- x_prj_diet[rowSums(x_prj_diet!=0)/ncol(x_prj_diet)>=0.1,]
    # x2 <- x_prj_diet[rowSums(x_prj_diet!=0)>0,]
    x2 <- data.frame(species=rownames(x2))
  }else{
    x2 <- data.frame(species=rownames(x12)[x_prj_diet!=0]  )
  }
  
  if(nrow(x2)!=0){
    x2$species_host <- sapply(strsplit(x2$species,split = "s__",fixed = T),"[",2)
    x2$species_host[grep("crass",x2$species_host ,ignore.case = T)]<- "crAssphage"
    x2$species_host <- sapply(strsplit(x2$species_host,split = "_",fixed = T),"[",1)
    x2<- data.frame(table(x2$species_host))
    colnames(x2) <- c("spe",diet)
    spe_x_prj_diets_zong <- merge(spe_x_prj_diets_zong,x2,by="spe",all = T)
  }
  
}
spe_x_prj_diets_zong2 <- pivot_longer(spe_x_prj_diets_zong ,cols = colnames(spe_x_prj_diets_zong)[2:ncol(spe_x_prj_diets_zong)],values_to = "diet_change_all" )
spe_x_prj_diets_zong2 <- spe_x_prj_diets_zong2[is.na(spe_x_prj_diets_zong2$diet_change_all)==F,]
spe_x_prj_diets_zong2$type <- paste(spe_x_prj_diets_zong2$spe,spe_x_prj_diets_zong2$name,sep="_")

spe_x_prj_diets_tongji2 <- pivot_longer(spe_x_prj_diets_tongji ,cols = colnames(spe_x_prj_diets_tongji)[2:ncol(spe_x_prj_diets_tongji)] ,values_to = "diet_change_diff")
spe_x_prj_diets_tongji2 <- spe_x_prj_diets_tongji2[is.na(spe_x_prj_diets_tongji2$diet_change_diff)==F,]
spe_x_prj_diets_tongji2$type <- paste(spe_x_prj_diets_tongji2$spe,spe_x_prj_diets_tongji2$name,sep="_")

spe_x_prj_diets_tongji2 <- merge(spe_x_prj_diets_tongji2,spe_x_prj_diets_zong2[3:4],by="type")
a <- data.frame(Bact=rownames(data_g_B))
a$phylum <- sapply(strsplit(a$Bact,split = "p__") ,"[",2)
a$phylum <- sapply(strsplit(a$phylum,split = "|",fixed = T) ,"[",1)
a$genus <- sapply(strsplit(a$Bact,split = "g__") ,"[",2)
a$genus <- sapply(strsplit(a$genus,split = "|",fixed = T) ,"[",1)
spe_x_prj_diets_tongji2 <- merge(spe_x_prj_diets_tongji2,unique(a[2:3]),by.x = "spe",by.y="genus",all.x = T)
spe_x_prj_diets_tongji2$phylum[grep("crAssphage",spe_x_prj_diets_tongji2$sp,ignore.case = T)] <- "Bacteroidetes"
spe_x_prj_diets_tongji2[is.na(spe_x_prj_diets_tongji2)] <- "Others"
colnames(spe_x_prj_diets_tongji2)[colnames(spe_x_prj_diets_tongji2)=="Type"] <- "Host"

spe_x_prj_diets_tongji3 <- aggregate(list(diet_change_diff=spe_x_prj_diets_tongji2$diet_change_diff,
                                          diet_change_all=spe_x_prj_diets_tongji2$diet_change_all)
                                     ,list(diet=spe_x_prj_diets_tongji2$name,
                                           phylum=spe_x_prj_diets_tongji2$phylum),
                                     sum)
spe_x_prj_diets_tongji3 <- spe_x_prj_diets_tongji3[spe_x_prj_diets_tongji3$diet_change_all>=5,]###筛选标准3，类型标准
spe_x_prj_diets_tongji3$percent <- spe_x_prj_diets_tongji3$diet_change_diff/spe_x_prj_diets_tongji3$diet_change_all*100
spe_x_prj_diets_tongji3 <- spe_x_prj_diets_tongji3[spe_x_prj_diets_tongji3$diet%in%
                                                     c("Baseline","Chow_diet","Control","Omnivores","Low_urbanization","Low_fiber_diet","Low_gluten_diet","Habitual")==F,]
spe_x_prj_diets_tongji3 <- spe_x_prj_diets_tongji3[spe_x_prj_diets_tongji3$diet%in%
                                                     c("Table_Food","Formula","Breastmilk","MIX")==F,]

spe_x_prj_diets_tongji3$group <- 1:nrow(spe_x_prj_diets_tongji3)
spe_x_prj_diets_tongji3 <-  pivot_longer(spe_x_prj_diets_tongji3[c("diet","phylum","percent","group")],cols = c("diet","phylum")  )
spe_x_prj_diets_tongji3$value <- gsub("_"," ",spe_x_prj_diets_tongji3$value )
spe_x_prj_diets_tongji3$value <- gsub("diet","",spe_x_prj_diets_tongji3$value )
a <- aggregate(list(percent=spe_x_prj_diets_tongji3$percent),by=list(value=spe_x_prj_diets_tongji3$value,name=spe_x_prj_diets_tongji3$name),sum)
a <- a[order(a$name,a$percent,decreasing = T),]
spe_x_prj_diets_tongji3$value <- factor(spe_x_prj_diets_tongji3$value,levels = c(unique(a$value)[unique(a$value)!="Others"],"Others"),labels = c(unique(a$value)[unique(a$value)!="Others"],"Others")   )
spe_x_prj_diets_tongji3$name <- factor(spe_x_prj_diets_tongji3$name,levels = unique(a$name),labels = unique(a$name))
 
#### Painting
mycolor <-  c("#DDEAC4","#85C3B9", "#7CA5C2","#E9C4DA",colorRampPalette(brewer.pal(9,'OrRd'))(20)[17:3],"#B7B3D1")
ggplot(spe_x_prj_diets_tongji3,
       aes(x = name,y = percent, stratum = value, alluvium = group, fill =  value, label =  value) ) +
  geom_flow(width = 0.2, knot.pos = 0.1,curve_type = "arctangent",size = 0.1) + #间隔宽度,#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
  geom_stratum(alpha = 0.8,color=NA,width = 0.18) + #width图中方块的宽度
  geom_text(stat = "stratum", size =4,color="black") +
  scale_fill_manual(values = mycolor) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)
ggsave("../fig1-E.pdf") 



