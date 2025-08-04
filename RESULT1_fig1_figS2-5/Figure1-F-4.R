##Script to Figure 1-F-4.

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
prj2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_prj.csv")
life_style_reads <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F-4.csv")

diet_controls <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_diet1.csv")
diet_controls2 <- read.csv("~/cooperation/202409zhaofengxiang/data/fig1F_diet2.csv")

crass_phage3 <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_crass_phage3")

############################## Data processing ##############################
#colnames(life_style_reads) <- c("Contig_number","Lytic","Temperate","Sample_id")
#life_style_reads <- life_style_reads[life_style_reads$Temperate!="Temperate",]

life_style_reads$Sample_id <- sapply(strsplit(life_style_reads$Sample_id,split = ".",fixed = T), "[",1)
life_style_reads$type <- "Temperate"
life_style_reads$type[life_style_reads$Lytic>0.5] <- "Lytic"
# life_style_reads2 <- aggregate(life_style_reads$type,by=list(life_style_reads$Sample_id),
#                                function(x){length(x[x=="Lytic"])/length(x[x=="Temperate"])})
life_style_reads2 <- aggregate(life_style_reads$type,by=list(life_style_reads$Sample_id),
                               function(x){length(x[x=="Lytic"])/length(x)})
colnames(life_style_reads2) <- c("Sample_id","Lifestyle_Bili")
life_style_reads2 <- merge(life_style_reads2,unique(prj2),by.x = "Sample_id",by.y = "Run")
life_style_reads2$Diet2 <- paste(life_style_reads2$Diet1,life_style_reads2$Nutrition1,sep="_")


unique(life_style_reads2$Diet1)
diet="Low_calorie_diet"
unique(prj2$Diet1[prj2$BioProject%in%unique(prj2$BioProject[prj2$Diet1==diet])  ])

a <- unique(life_style_reads2[,c("BioProject","Diet1")])
table(life_style_reads2$Diet1)
a<- aggregate(1:nrow(life_style_reads2),
              by=list(life_style_reads2$BioProject),function(x) 
                if(length(unique(life_style_reads2$Diet1[x]))>=2   )
                {kruskal.test(Lifestyle_Bili~Diet1,life_style_reads2[x,])[["p.value"]]}else{NA} )

df_rongyuanliejies <- data.frame(matrix(ncol=0,nrow = 0))
i=6
for (i in 1:nrow(diet_controls2) ){
  diet=diet_controls2$case[i]
  diet_control <- unlist(strsplit(diet_controls2$control[i],split = "/"))
  
  prjs <- unique(prj2$BioProject[prj2$Diet1==diet])
  j=1
  for (j in 1:length(prjs)){
    prj22 <- prj2[prj2$BioProject%in%prjs[j],]
    prj22 <- prj22[prj22$Diet1%in%c(diet_control,diet),]
    
    prj_case <- prj22[prj22$Diet1%in%diet,]
    if(diet_control[1]%in%prj22$Diet1){
      prj_control <- prj22[prj22$Diet1%in%diet_control[1],]
    }else{
      prj_control <- prj22[prj22$Diet1%in%diet_control[2],]
    }
    
    if(nrow(prj_control[prj_control$Run%in%life_style_reads2$Sample_id,])>0 & 
       nrow(prj_case[prj_case$Run%in%life_style_reads2$Sample_id,] )>0    ){
      median=median(life_style_reads2$Lifestyle_Bili[life_style_reads2$Sample_id%in%prj_case$Run])-
        median(life_style_reads2$Lifestyle_Bili[life_style_reads2$Sample_id%in%prj_control$Run])
      
      p_value=wilcox.test(x=life_style_reads2$Lifestyle_Bili[life_style_reads2$Sample_id%in%prj_case$Run],
                          y=life_style_reads2$Lifestyle_Bili[life_style_reads2$Sample_id%in%prj_control$Run])[["p.value"]]
      df_rongyuanliejie <- data.frame(diet=diet,diet_control=diet_control,life_style_change=median,p_value=p_value)
      df_rongyuanliejies <- rbind(df_rongyuanliejies,df_rongyuanliejie)
    }
    
  }
  
  
}


df_rongyuanliejies$life_style_change_Type <- "Reduce"
df_rongyuanliejies$life_style_change_Type[df_rongyuanliejies$life_style_change>0] <- "Increase"
df_rongyuanliejies$life_style_change_p_value <- "Sig"
df_rongyuanliejies$life_style_change_p_value[df_rongyuanliejies$p_value>0.05] <- "No_sig"
df_rongyuanliejies$Type <- paste(df_rongyuanliejies$life_style_change_p_value,df_rongyuanliejies$life_style_change_Type,sep="_")
lecels<- rev(c("Sig_Reduce","Sig_Increase","No_sig_Reduce","No_sig_Increase"))
df_rongyuanliejies$Type <- factor(df_rongyuanliejies$Type,levels = lecels,labels = lecels  )
colors<- rev(c("#3288BD","#D53E4F","#DEEBF7","#FEE0D2"))
names(colors) <- lecels

df_rongyuanliejies$diet <- gsub("_"," ",df_rongyuanliejies$diet)
df_rongyuanliejies$diet <- gsub("diet","",df_rongyuanliejies$diet)
diets<- c("High fat ","Vegatarian","EEN ","High fiber ","Mediterranean ","WTP ","High urbanization","MRE")
df_rongyuanliejies <- df_rongyuanliejies[df_rongyuanliejies$diet%in%diets,]

df_rongyuanliejies2 <- aggregate(list(Reduce=df_rongyuanliejies$life_style_change_Type),
                                 by=list(diet=df_rongyuanliejies$diet),
                                 function(x) length(x[x=="Reduce"])/length(x)   )
df_rongyuanliejies2$Increase <- 1-df_rongyuanliejies2$Reduce
df_rongyuanliejies2 <- pivot_longer(df_rongyuanliejies2,cols = colnames(df_rongyuanliejies2)[2:3],
                                    names_to = "life_style_change_Type",
                                    values_to = "percent")

df_rongyuanliejies2$diet <- factor(df_rongyuanliejies2$diet,levels =c("High fat ","Vegatarian","EEN ","High fiber ","Mediterranean ","WTP ","High urbanization","MRE") ,
                                   labels =c("High fat ","Vegatarian","EEN ","High fiber ","Mediterranean ","WTP ","High urbanization","MRE")  )

df_rongyuanliejies2$life_style_change_Type <- factor(df_rongyuanliejies2$life_style_change_Type,
                                                     levels = c("Reduce","Increase"),
                                                     labels = c("Lytic phage ratio reduces","Lytic phage ratio increases"))
#更改展示顺序
levels(df_rongyuanliejies2$life_style_change_Type) 
df_rongyuanliejies2$life_style_change_Type <- factor(df_rongyuanliejies2$life_style_change_Type,levels =c("Lytic phage ratio increases","Lytic phage ratio reduces"))
levels(df_rongyuanliejies2$life_style_change_Type) 

colors2<- c("#F48892","#91CAE8")
names(colors2) <- c("Lytic phage ratio increases","Lytic phage ratio reduces")

############################## [Fig1-F-4] Painting #####################################
ggplot(df_rongyuanliejies2, aes(x = "", y = percent, fill =life_style_change_Type)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) +
  facet_wrap(.~diet,nrow =1)+ 
  theme_test()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position="none")+ #输出pdf大小20X3
  scale_fill_manual(values = colors2)+
  xlab("")+ylab("")
ggsave("../fig1-F-4.pdf") 