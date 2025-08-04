##Script to Figure 6-A and Supplementary fig8-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("tidyverse","tidyr","treemap","treemapify","ggalluvial","RColorBrewer","ggplot2")
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
b <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6A_2.csv",row.names = 1)
ab <- read.csv("~/cooperation/202409zhaofengxiang/data/fig6A_1.csv")

prj <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv")
prj1 <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_prj_data")

As <- read.table("~/cooperation/202409zhaofengxiang/data/fig6_As")
prj2 <- read.table("~/cooperation/202409zhaofengxiang/data/fig6_prj_data")

data_g_B <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_B")

############################## Data processing ##############################
b_spe <- b
fenleis<-c("Run","BioProject","Diet1","Diet","Nutrition1","Nutrition","Antibiotic","Continent","Host","Age","Age2","Sex" )
i=5
prj_index <- i
assign(fenleis[i],b_spe[,grep(fenleis[i],colnames(b_spe),value = F)])
Diet1<- Nutrition1
colnames(Diet1) <- gsub("Abundance.","",colnames(Diet1))
x <- Diet1

ab<- unique(ab[c("prj_name","type3")])
ab <- ab[ab$type3=="Nutrient addition",]
search <- sapply(strsplit(colnames(x),split = ".",fixed = T),"[",1)
serach2 <- sapply(strsplit(colnames(x),split = ".",fixed = T),"[",2)
serach2[serach2%in%as.character(1:9)==F] <- ""
search <- paste(search,serach2,sep=".")
x <- x[rowSums(is.na(x)==F)>0,]


x[is.na(x)==T] <- 0
prj1[prj1==""] <- NA
prjs <- unique(prj1$BioProject)
prjs_clean <- paste0(prjs,".",fenleis[i])

colnames(x)<- gsub('*_4.20_Diet1_', '_', colnames(x))
x <- x[apply(x, 1, function(row) any(!is.na(row))),]
x <- x[rowSums(x!=0)>0,]

y <- apply(x, 2, function(x) table(x[x!=0]))
z <- apply(x, 1, function(x) table(t(x[x!=0])))


############################## [figS8-C] Processing and Painting #####################################
x <- Diet1
i<- prj_index
x <- x[rowSums(is.na(x)==F)>0,]
x[is.na(x)] <- 0
x2 <- x
x2<- x2[grepl("s_",rownames(x2)),]
x2<- x2[rowSums(x2!=0)>=2,]
x2<- data.frame(diff=colSums(x2!=0),spes=colnames(x2))

x4<- data.frame(all=colSums(is.na(As[2:ncol(As)])==F) ,spes=colnames(As)[2:ncol(As)] )
x4$spes <- paste0(x4$spes,".",fenleis[prj_index])
x2 <- merge(x2,x4,by="spes")
x2$percent <- x2$diff/x2$all*100
x2 <- x2[x2$percent>=0.1,]

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
x4 <- x4[x4$all>=2,]
x4 <- merge(x4,unique(prj[c("Nutrition1","Nutrition")]),by.x = "diet",by.y = "Nutrition1",all=T)
x4 <- x4[x4$Nutrition%in%c("Washout","Control","Baseline")==F,]
x4 <- x4[is.na(x4$all)==F,]
x4$diet <- gsub("_"," ",x4$diet)
df <- x4
df <- df[is.na(df$Nutrition)==F,]

treemap(df,index = colnames(x4)[c(5,1)],vSize="percent",
        palette = "Set2",border.lwds = 0.8,border.col = "white") 

###### Painting
ggplot(df, aes(area = percent, fill = Nutrition, label = diet,
               subgroup = Nutrition)) + #子分组映射
  scale_fill_manual(values = c(
    "Polyphenol"="#ABD1BC",
    "Polysaccharide"="#FCB6A5",
    "Oligosaccharide"="#E3BBED",
    "Others"="#F1F1F1",
    "Protein"="#BED0F9",
    "Probiotics"="#cCcc99",
    "Alkohol"="#9DA7CB",
    "Sweetener"="#F1F1F1"))+
  geom_treemap(alpha=0.5) +
  geom_treemap_subgroup_border(color = 'white',size = 4) + #子分组边界线绘制
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.3, colour = "black", fontface = "italic", min.size = 0) +
  geom_treemap_text(colour = "black", place = "topleft",fontface = "bold",size = 6)+
  theme(legend.position = "none")
ggsave("../figS8-C.pdf")


############################## [Fig6-A] Processing and Painting #####################################
x <- Diet1
x <- x[rowSums(is.na(x)==F)>0,]
x[is.na(x)] <- 0
i <- prj_index
diets<- unique(prj2[,fenleis[i]])
diets <- diets[is.na(diets)==F]
diets_index <- 1
colnames(x) <- gsub("Abundance.","",colnames(x) )

spe_x_prj_diets_tongji <-  data.frame(matrix(ncol = 1,nrow=0))
colnames(spe_x_prj_diets_tongji) <-"spe"
spe_x_prj_diets_zong<-  data.frame(matrix(ncol = 1,nrow=0))
colnames(spe_x_prj_diets_zong) <-"spe"

diets_index=1
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
  
  if(nrow(x2)!=0&length(grep("s__",x2$species))!=0   ){
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
                                                     c("Baseline","Chow_diet","Control","Omnivores","Low_urbanization","Low_fiber_diet","High_gluten_diet","Low_gluten_diet")==F,]
spe_x_prj_diets_tongji3 <- merge(spe_x_prj_diets_tongji3,unique(prj2[c("Nutrition1","Nutrition")]),by.x = "diet",by.y="Nutrition1")
spe_x_prj_diets_tongji3 <- spe_x_prj_diets_tongji3[spe_x_prj_diets_tongji3$Nutrition!="Water",]
spe_x_prj_diets_tongji3$group <- 1:nrow(spe_x_prj_diets_tongji3)
spe_x_prj_diets_tongji3 <- spe_x_prj_diets_tongji3[spe_x_prj_diets_tongji3$Nutrition%in%c("Control","Washout")==F,]
spe_x_prj_diets_tongji3 <-  pivot_longer(spe_x_prj_diets_tongji3[c("Nutrition","phylum","percent","group")],cols = c("Nutrition","phylum")  )

a <- aggregate(list(percent=spe_x_prj_diets_tongji3$percent),by=list(value=spe_x_prj_diets_tongji3$value,name=spe_x_prj_diets_tongji3$name),sum)
a <- a[order(a$name,a$percent,decreasing = T),]
spe_x_prj_diets_tongji3$value <- factor(spe_x_prj_diets_tongji3$value,levels = c(unique(a$value)[unique(a$value)!="Others"],"Others"),labels = c(unique(a$value)[unique(a$value)!="Others"],"Others")   )
spe_x_prj_diets_tongji3$name <- factor(spe_x_prj_diets_tongji3$name,levels = unique(a$name),labels = unique(a$name))
levels(spe_x_prj_diets_tongji3$value)

###### Painting
mycolor <-  c("#DDEAC4","#E9C4DA","#85C3B9", "#7CA5C2",
              colorRampPalette(brewer.pal(9,'Blues'))(9)[9:1],"#B7B3D1")
ggplot(spe_x_prj_diets_tongji3,
       aes(x = name,y = percent, stratum = value, alluvium = group, fill =  value, label =  value) ) +
  geom_flow(width = 0.2, knot.pos = 0.1,curve_type = "arctangent",size = 0.1) + #间隔宽度,#曲线形状，有linear、cubic、quintic、sine、arctangent、sigmoid几种类型可供调整
  geom_stratum(alpha = 0.8,color=NA,width = 0.18) + #width图中方块的宽度
  geom_text(stat = "stratum", size =4,color="black") +
  scale_fill_manual(values= mycolor)+ #values=mycolor)+
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)
ggsave("../fig6-A.pdf")





