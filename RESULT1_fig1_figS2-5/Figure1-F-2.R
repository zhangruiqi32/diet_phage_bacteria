##Script to Figure 1-F-2.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggplot2")
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
x_Viru <- read.csv("../data/fig1F-2_1.csv", row.names = 1)
x_Bact <- read.csv("../data/fig1F-2_2.csv", row.names = 1)

As <- read.table("../data/fig1_As")
prj2 <- read.table("../data/fig1_prj_data")
diet_controls <- read.csv("../data/fig1F_diet1.csv")

############################## Data processing ##############################
x_Bact2 <- x_Bact[grep("g__",rownames(x_Bact)),]
x_Bact2 <- x_Bact2[-grep("s__",rownames(x_Bact2)),]
x_Bact2 <-x_Bact2[rowSums(is.na(x_Bact2))!=0,]
a <- cbind(rownames(x_Bact2),x_Bact2)
x_Bact2 <- pivot_longer(a,names_to = "prjs",values_to = "LEFSE_result",
                        cols = colnames(x_Bact2))
x_Bact2 <- x_Bact2[is.na(x_Bact2$LEFSE_result)==F,]
colnames(x_Bact2)[1] <- "Spes"
x_Bact2$genus <- sapply(strsplit(x_Bact2$Spes,split = "g__",fixed = T),"[",2 )
x_Bact2$genus <- sapply(strsplit(x_Bact2$genus,split = "_",fixed = T),"[",1 )
x_Bact2$prj_diet_spes <- paste(x_Bact2$prjs,x_Bact2$genus,sep="_")
colnames(x_Bact2) <- paste(colnames(x_Bact2),"Bact",sep="_")
x_Bact2 <- unique(x_Bact2)


x_Viru2 <- x_Viru[grep("s__",rownames(x_Viru)),]
x_Viru2 <-x_Viru2[rowSums(is.na(x_Viru2))!=0,]
a <- cbind(rownames(x_Viru2),x_Viru2)
x_Viru2 <- pivot_longer(a,names_to = "prjs",values_to = "LEFSE_result",
                        cols = colnames(x_Viru2))
x_Viru2 <- x_Viru2[is.na(x_Viru2$LEFSE_result)==F,]
colnames(x_Viru2)[1] <- "Spes"
x_Viru2$Bact_host <- sapply(strsplit(x_Viru2$Spes,split = "s__",fixed = T),"[",2 )
x_Viru2$Bact_host <- sapply(strsplit(x_Viru2$Bact_host,split = "_",fixed = T),"[",1 )
x_Viru2$Bact_host <- sapply(strsplit(x_Viru2$Bact_host,split = ".",fixed = T),"[",1 )
x_Viru2$Bact_host[grep("crass",x_Viru2$Spes,ignore.case = T)] <-"Bacteroides"
x_Viru2$prj_diet_spes <- paste(x_Viru2$prjs,x_Viru2$Bact_host,sep="_")
x_Viru2 <- unique(x_Viru2)

x_Viru3 <- merge(x_Viru2,x_Bact2,by.x = "prj_diet_spes",by.y = "prj_diet_spes_Bact")
x_Viru3 <- data.frame(x_Viru3)

prjs_diet <- unique(prj2[,c("Diet1","BioProject")])
prjs_diet <- prjs_diet[prjs_diet$Diet1%in%diet_controls$case,]
prjs_diet$BioProject <- paste0(prjs_diet$BioProject,".Diet1")
x_Viru3 <- merge(x_Viru3,prjs_diet,by.x = "prjs",by.y = "BioProject")
x_Viru3$Change_type<- "Opposite direction"
x_Viru3$Change_type[x_Viru3$LEFSE_result==x_Viru3$LEFSE_result_Bact] <- "Same direction"
x_Viru3$prjs_spes <- paste0(x_Viru3$Spes,x_Viru3$prjs,sep="_")
x_Viru3 <- aggregate(list(number=1:nrow(x_Viru3)),
                     by=list(prjs=x_Viru3$prjs,Bact_host=x_Viru3$Bact_host,Change_type=x_Viru3$Change_type),length)
# x_Viru3 <- x_Viru3[x_Viru3$number>=2,]

As2 <- pivot_longer(As,names_to = "prjs",values_to = "LEFSE_result",
                    cols = colnames(As)[2:ncol(As)]   )
As2 <- As2[is.na(As2$LEFSE_result)==F,]
As2$Bact_host <- sapply(strsplit(As2$spes,split = "s__",fixed = T),"[",2 )
As2$Bact_host <- sapply(strsplit(As2$Bact_host,split = "_",fixed = T),"[",1 )
As2$Bact_host <- sapply(strsplit(As2$Bact_host,split = ".",fixed = T),"[",1 )
As2$Bact_host[grep("crass",As2$spes,ignore.case = T)] <-"Bacteroides"
As2$prjs <- paste0(As2$prjs,".Diet1")
As2$spes <- gsub("|",".",As2$spes,fixed = T)
As2$prjs_spes <- paste0(As2$spes,As2$prjs,sep="_")
As2 <- As2[As2$prjs_spes%in%x_Viru3$prjs_spes==F,]

As2 <- aggregate(list(number=1:nrow(As2)),
                 by=list(prjs=As2$prjs,Bact_host=As2$Bact_host),length )
As2<- As2[As2$Bact_host%in%x_Viru3$Bact_host,]


#prj2 <- prj_data
prjs_diet <- unique(prj2[,c("Diet1","BioProject")])
prjs_diet <- prjs_diet[prjs_diet$Diet1%in%diet_controls$case,]
x_Viru3 <- data.frame(x_Viru3)
prjs_diet$BioProject <- paste0(prjs_diet$BioProject,".Diet1")
x_Viru3 <- merge(x_Viru3,prjs_diet,by.x = "prjs",by.y = "BioProject")

x_Viru3$x <- paste(x_Viru3$Diet1,x_Viru3$prjs,x_Viru3$Bact_host,sep="_")


x_Viru3$Diet1 <- gsub("_"," ",x_Viru3$Diet1)
x_Viru3$Diet1 <- gsub("diet","",x_Viru3$Diet1)
x_Viru3$Diet1 <- factor(x_Viru3$Diet1,levels = c("High fat ","Vegatarian","EEN ","High fiber ","Mediterranean ","WTP ","High urbanization","MRE"),
                        labels = c("High fat ","Vegatarian","EEN ","High fiber ","Mediterranean ","WTP ","High urbanization","MRE")   )

a <- aggregate(list(number=x_Viru3$number),
               by=list(Change_type=x_Viru3$Change_type,x=x_Viru3$x,Diet1=x_Viru3$Diet1),
               sum)
a <- pivot_wider(a,names_from = "Change_type",values_from = "number")
a[is.na(a)] <- 0
a$percent <- a$`Same direction`/(a$`Same direction`+a$`Opposite direction`)
a$percent2 <- a$`Opposite direction`/(a$`Same direction`+a$`Opposite direction`)
a <- a[order(a$percent,-a$percent2),]
x_Viru3$x <- factor(x_Viru3$x,levels =a$x ,labels = a$x)

############################## [Fig1-F-2] Painting #####################################
ggplot(x_Viru3,aes(x=x,y=log10(number+1),fill=Change_type))+
  geom_bar(stat = "identity",position = "fill",width = 1)+
  theme_test()+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        legend.position="none")+ #输出pdf大小20X3
  facet_wrap(.~Diet1,scales = "free_x",nrow = 1)+
  scale_fill_manual(values = c("#91CAE8","#F48892","gray"))
ggsave("../fig1-F-2.pdf") 
