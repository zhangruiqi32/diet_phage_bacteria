##Script to Figure 6-B.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("RColorBrewer","graphics","ggplot2")
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
dfs_crass2 <- read.csv("../data/fig6B_1.csv")
crass_phage2 <- read.csv("../data/fig6B_2.csv")
diet_case <- read.csv("../data/fig6B_diet.csv")
prj <- read.csv("../data/fig1F_prj.csv")

############################## Data processing ##############################
dfs_xiangguanxinngs <- dfs_crass2
dfs_xiangguanxinngs <- dfs_xiangguanxinngs[dfs_xiangguanxinngs$p_value<0.05,]

dfs_xiangguanxinngs$Diet <- sapply(strsplit(dfs_xiangguanxinngs$`PRJ_DIET`,split="-"),"[",2)
dfs_xiangguanxinngs$Diet <- gsub("/","",dfs_xiangguanxinngs$Diet,fixed = T)
dfs_xiangguanxinngs$Project <- sapply(strsplit(dfs_xiangguanxinngs$`PRJ_DIET`,split="-"),"[",1)

diet_case$`PRJ_DIET` <- paste0(diet_case$Project,"-",diet_case$Diet,"/")
dfs_xiangguanxinngs <- merge(dfs_xiangguanxinngs,diet_case[c("type","PRJ_DIET")],by="PRJ_DIET")

#####共享相关性-饮食分组
dfs_xiangguanxinngs2 <- dfs_xiangguanxinngs
dfs_xiangguanxinngs2$phylum[grep("crass",dfs_xiangguanxinngs2$Viruse,ignore.case = T)] <- "crAssphage"
dfs_xiangguanxinngs2 <- dfs_xiangguanxinngs2[dfs_xiangguanxinngs2$phylum%in%crass_phage2$phylum,]
dfs_xiangguanxinngs2 <- dfs_xiangguanxinngs2[order(dfs_xiangguanxinngs2$type),]

dfs_xiangguanxinngs2 <- unique(dfs_xiangguanxinngs2)
dfs_xiangguanxinngs2 <- dfs_xiangguanxinngs2[dfs_xiangguanxinngs2$type=="case",]
dfs_xiangguanxinngs2$BV <- paste(dfs_xiangguanxinngs2$Bact,dfs_xiangguanxinngs2$Viruse,sep="_")

dfs_xiangguanxinngs2 <- merge(dfs_xiangguanxinngs2,unique(prj[c("Nutrition1","Nutrition")]),by.x = "Diet",by.y = "Nutrition1")
dfs_xiangguanxinngs2$BV_type <- paste(dfs_xiangguanxinngs2$BV ,dfs_xiangguanxinngs2$Nutrition,sep="_")
dfs_xiangguanxinngs2 <- dfs_xiangguanxinngs2[!duplicated(dfs_xiangguanxinngs2$BV_type),]

dfs_xiangguanxinngs3 <- pivot_wider(dfs_xiangguanxinngs2[c("BV","R2","Nutrition","phylum")],
                                    names_from  = "Nutrition",values_from = "R2")
dfs_xiangguanxinngs3 <- dfs_xiangguanxinngs3[is.na(dfs_xiangguanxinngs3$phylum)==F,]
dfs_quan2 <- data.frame(dfs_xiangguanxinngs3,Bact_p=dfs_xiangguanxinngs3$phylum  )


j=1
as <- data.frame(matrix(ncol=1,nrow=0))
for (i in 1:length(unique(dfs_quan2$Bact_p))){
  Bact_p=unique(dfs_quan2$Bact_p)[i]
  ab <- dfs_quan2[dfs_quan2$Bact_p==Bact_p,3:(ncol(dfs_quan2)-1)]
  for (j in 1:ncol(ab) ){
    ac <- ab[is.na(ab[,j])==F,]
    a <- data.frame(shared_correlation=apply(ac,2, function(x) length(x[is.na(x)==F])),
                    diet=colnames(ab)[j],
                    Bact_p=Bact_p,
                    shared_diets=colnames(ac))
    as <- rbind(as,a) 
  }
}

as <- pivot_wider(as,names_from= "Bact_p",values_from="shared_correlation")
as <- as[,colSums(as!=0)>0]
diets <- gsub("+",".",x = unique(dfs_xiangguanxinngs2$`PRJ_DIET`),fixed = T)
diets <- gsub("/",".",x = diets,fixed = T)
diets <- gsub("-",".",x = diets,fixed = T)
unique(as$diet)

############################## [fig6-B] Processing and Painting #####################################
stars(log10(as[3:ncol(as)]+1),key.loc = c(-1,2),
      draw.segments = T, col.segments=colorRampPalette(c("#E2D4D4", "#5D0E12"))(13)[13:1])

ggsave("../fig6-B.pdf")

