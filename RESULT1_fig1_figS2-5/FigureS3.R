##Script to Supplementary fig3.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggforce","lemon","RColorBrewer","ggplot2")
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
data_g_B <- read.table("~/cooperation/202409zhaofengxiang/data/fig1_data_g_B")

spes_tongji_BACT <- read.csv("~/cooperation/202409zhaofengxiang/data/figS3_1.csv")
spes_tongji_VIRU <- read.csv("~/cooperation/202409zhaofengxiang/data/figS3_2.csv")

diet1 <- read.csv("~/cooperation/202409zhaofengxiang/data/figS3_diet.csv",sep=",", header = F)

############################## Data processing ##############################
# spes_tongji_VIRU$percent <- spes_tongji_VIRU$percent*100
spes_tongji_BACT <- spes_tongji_BACT[grep("s__",spes_tongji_BACT$spes,fixed = T),]
spes_tongjis <- rbind(spes_tongji_VIRU,spes_tongji_BACT)
spes_tongjis <- spes_tongjis[spes_tongjis$diet%in%c("Table_Food","Formula","Breastmilk","MIX","High_fermented_diet","High_gluten_diet","Low_gluten_diet","Whole_grain_diet")==F,]
spes_tongjis <- spes_tongjis[spes_tongjis$diet%in%c("Baseline","Chow_diet","Control","Omnivores","Low_urbanization","Low_fiber_diet","High_protein_diet")==F,]
spes_tongjis <- spes_tongjis[spes_tongjis$diet%in%c("Low_carbon_water_diet","Low_calorie_diet","HFHS_diet")==F,]
spes_tongjis <- spes_tongjis[spes_tongjis$Bact_host%in%c("Bacillus")==F,]
spes_tongjis$Bact_host2 <- as.character(spes_tongjis$Bact_host)
spes_tongjis$Bact_host2[grep("crass",spes_tongjis$Bact_host2,ignore.case = T)] <-"Bacteroides"
a<- data.frame(table(spes_tongjis[c("type","Bact_host2")]))
a <- a[a$Freq!=0,]
a <- a[duplicated(a$Bact_host2),]
spes_tongjis <- spes_tongjis[spes_tongjis$Bact_host2%in%a$Bact_host2,]


a <- data.frame(Bact=rownames(data_g_B))
a$phylum <- sapply(strsplit(a$Bact,split = "p__") ,"[",2)
a$phylum <- sapply(strsplit(a$phylum,split = "|",fixed = T) ,"[",1)
a$genus <- sapply(strsplit(a$Bact,split = "g__") ,"[",2)
a$genus <- sapply(strsplit(a$genus,split = "|",fixed = T) ,"[",1)
spes_tongjis<- merge(spes_tongjis,unique(a[2:3]),by.x = "Bact_host",by.y="genus",all.x = T)
spes_tongjis$phylum[spes_tongjis$Bact_host=="crAssphage"] <- "Bacteroidetes"
spes_tongjis$phylum[is.na(spes_tongjis$phylum)] <- "Others"

# spes_tongjis <- spes_tongjis[order(spes_tongjis$phylum,spes_tongjis$Bact_host,),]
colnames(diet1) <- c("Diet","Type")
spes_tongjis <- merge(spes_tongjis,diet1,by.x = "diet",by.y="Diet")


spes_tongjis$group2 <- paste(spes_tongjis$phylum,spes_tongjis$type,spes_tongjis$Bact_host2,spes_tongjis$Type,spes_tongjis$diet,sep="_")
spes_tongjis <- spes_tongjis [order(spes_tongjis$group2,spes_tongjis$percent),]
spes_tongjis$Bact_host2 <- factor(spes_tongjis$Bact_host2,levels = unique(spes_tongjis$Bact_host2),labels = unique(spes_tongjis$Bact_host2) )
spes_tongjis$die2 <- spes_tongjis$diet
spes_tongjis$diet <- gsub("_"," ",spes_tongjis$diet)
spes_tongjis$diet <- gsub("diet","",spes_tongjis$diet)
spes_tongjis$diet <- factor(spes_tongjis$diet,levels = unique(spes_tongjis$diet),labels = unique(spes_tongjis$diet) )


a<- data.frame(unlist(aggregate(1:nrow(spes_tongjis),by=list(spes_tongjis$group2),function(x) 1:length(x) )[2] ) )
spes_tongjis$number <- unlist(a[1])
xx<- table(spes_tongjis$group2)
spes_tongjis <- spes_tongjis[spes_tongjis$group2%in%names(xx[xx>=2]),]
spes_tongjis$number2 <- factor(spes_tongjis$number,
                               levels = unique(spes_tongjis$number), 
                               labels = unique(spes_tongjis$number)  )


############################## [FigS3] Painting #####################################
ggplot(spes_tongjis,aes(x=number,y=percent,fill=type )    )+
  # facet_grid(Bact_host2~diet,scales = "free",space = "free_x",drop = T,shrink = T)+
  # facet_rep_grid(Bact_host2~diet, scales='free')+
  facet_grid(Bact_host2~diet, scales='free',space = "free")+
  # facet_rep_wrap(Bact_host2~diet, scales='free',repeat.tick.labels = 'all' )+
  # facet_col(Bact_host2~diet,scales = "free_x")+
  geom_bar(
    stat="identity",linewidth=0.8,width = 1.05)+
  geom_hline(yintercept = 0,linewidth=0.1)+
  theme_bw()+
  scale_fill_manual(values = c(brewer.pal(11,"Spectral")[10] ,brewer.pal(11,"Spectral")[2]  ))+
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        # strip.text.x  = element_text(angle = 90),
        # strip.text.y  = element_text(angle =0),
        panel.grid=element_blank(),
        axis.ticks = element_blank() ,
        strip.text.y.right = element_text(angle = 0,hjust = 0),
        strip.background = element_rect(fill="white",color = "white"),
        strip.text.x  = element_text(#angle = 90,
          hjust = 0.5,vjust = 0.5)
        # axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
        # panel.background=element_blank(),
        # panel.grid = element_blank()
        # panel.spacing = unit(0.01, "cm", data = NULL)
  )+xlab("")+ylab("")
ggsave("../figS3.pdf") 
