##Script to Supplementary fig8-D.

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

diet_case <- read.csv("~/cooperation/202409zhaofengxiang/data/figS8D_dietcase.csv")

############################## Data processing ##############################
life_style_reads$Sample_id <- sapply(strsplit(life_style_reads$Sample_id,split = ".",fixed = T), "[",1)
life_style_reads$type <- "Temperate"
life_style_reads$type[life_style_reads$Lytic>0.5] <- "Lytic"

life_style_reads2 <- aggregate(life_style_reads$type,by=list(life_style_reads$Sample_id),
                               function(x){length(x[x=="Lytic"])/length(x)})
colnames(life_style_reads2) <- c("Sample_id","Lifestyle_Bili")
life_style_reads2 <- merge(life_style_reads2,unique(prj2),by.x = "Sample_id",by.y = "Run")


diet_case$Project <- gsub("_nutrition","",diet_case$Project)
diet_case$`PRJ-DIET` <- paste0(diet_case$Project,"-",diet_case$Diet)
life_style_reads2$`PRJ-DIET` <- paste0(life_style_reads2$BioProject,"","-",life_style_reads2$Nutrition1)
life_style_reads2 <- merge(life_style_reads2,diet_case[c("type","PRJ-DIET")],by="PRJ-DIET")

a<- diet_case[diet_case$type=="case",c("Project","Diet")]
colnames(a) <- c("Project","group")
life_style_reads2 <- merge(life_style_reads2,a,by.x="BioProject",by.y  ="Project")
life_style_reads2$group[life_style_reads2$type=="case"&life_style_reads2$group!=life_style_reads2$Nutrition1]<- NA
life_style_reads2 <- life_style_reads2[is.na(life_style_reads2$group)==F,]

df_rongyuanliejies<- aggregate(list(life_style_change_p_value=1:nrow(life_style_reads2)),
              by=list(group=life_style_reads2$group,BioProject=life_style_reads2$BioProject),function(x)
  if(  length(unique(life_style_reads2$type[x]))==2  )
    {wilcox.test(Lifestyle_Bili~type,life_style_reads2[x,])[["p.value"]]}else{NA}
  )
a <- aggregate(list(life_style_change=1:nrow(life_style_reads2)),
                by=list(life_style_reads2$group,life_style_reads2$BioProject),function(x)
  if(  length(unique(life_style_reads2$type[x]))==2  )
  { y <- life_style_reads2[x,]
    median(y$Lifestyle_Bili[y$type=="case"])-
      median(y$Lifestyle_Bili[y$type=="control"])  }else{NA}
)
df_rongyuanliejies$life_style_change <- a$life_style_change

 
prj<- prj2
df_rongyuanliejies<- merge(df_rongyuanliejies,unique(prj[c("Nutrition1","Nutrition")]),by.x = "group",by.y = "Nutrition1")
df_rongyuanliejies$diet<- df_rongyuanliejies$Nutrition

df_rongyuanliejies$life_style_change_Type <- "Reduce"
df_rongyuanliejies$life_style_change_Type[df_rongyuanliejies$life_style_change>0] <- "Increase"
df_rongyuanliejies$life_style_change_p_value <- "Sig"
df_rongyuanliejies$life_style_change_p_value[df_rongyuanliejies$p_value>0.05] <- "No_sig"
df_rongyuanliejies$Type <- paste(df_rongyuanliejies$life_style_change_p_value,df_rongyuanliejies$life_style_change_Type,sep="_")
lecels<- rev(c("Sig_Reduce","Sig_Increase","No_sig_Reduce","No_sig_Increase"))
df_rongyuanliejies$Type <- factor(df_rongyuanliejies$Type,levels = lecels,labels = lecels  )

df_rongyuanliejies2 <- aggregate(list(Reduce=df_rongyuanliejies$life_style_change_Type),
                                 by=list(diet=df_rongyuanliejies$diet),
                                 function(x) length(x[x=="Reduce"])/length(x)   )
df_rongyuanliejies2$Increase <- 1-df_rongyuanliejies2$Reduce
df_rongyuanliejies2 <- pivot_longer(df_rongyuanliejies2,cols = colnames(df_rongyuanliejies2)[2:3],
                                    names_to = "life_style_change_Type",
                                    values_to = "percent")

df_rongyuanliejies2$life_style_change_Type <- factor(df_rongyuanliejies2$life_style_change_Type,
                                                     levels = c("Reduce","Increase"),
                                                     labels = c("Lytic phage ratio reduces","Lytic phage ratio increases"))

############################## [FigS8-D] Painting #####################################
colors2<- c("#97B3D1","#DD7F88")
names(colors2) <- c("Lytic phage ratio reduces","Lytic phage ratio increases")

ggplot(df_rongyuanliejies2, aes(x = "", y = percent, fill =life_style_change_Type)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) +
  facet_wrap(.~diet,nrow =1)+
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = colors2)+
  xlab("")+ylab("")
write.csv(df_rongyuanliejies2,"~/cooperation/202409zhaofengxiang/code/figS8D_paintdata.csv",row.names = F)

ggsave("../figS8-D.pdf") 