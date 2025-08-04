##Script to Figure 1-C.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("vegan","sampling", "ggplot2")
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


############################## Data loading  ##############################
adoniss <- read.csv("../data/fig1C_1.csv")
adoniss2 <- read.csv("../data/fig1C_2.csv")

############################## [Fig1-C] Processing and Painting #####################################
adoniss2 <- merge(adoniss ,unique(adoniss2[c("prj_name","type2","type3")]),by="prj_name")
adoniss2$x <- paste(adoniss2$diets,adoniss2$prj_name,sep = "_")
 
  i=2
  adoniss3 <- adoniss2[adoniss2$type3==unique(adoniss2$type3)[i],]
  a <- aggregate(adoniss3$R2,by=list(adoniss3$type2),max)
  a <- a[order(a$x,decreasing = T),]
  adoniss3$type2 <- factor(adoniss3$type2,levels=c(a$Group.1),labels=c(  a$Group.1 )   )
  adoniss3 <- adoniss3[order(adoniss3$type2,-adoniss3$R2),]
  adoniss3$x <- factor(adoniss3$x,levels=adoniss3$x,labels=adoniss3$x   )
  adoniss3$p_value <- "sig"
  adoniss3$p_value[adoniss3$Pr..F.>=0.05] <- "sig_no"
  ggplot(adoniss3,aes(x=x,y=R2,fill=p_value))+
        geom_bar(stat = "identity")+
        facet_grid(.~type2,space = "free",scales = "free_x")+
        theme_classic()+
        scale_fill_manual(values = c("#2C7FB8" ,"gray"))+
        theme(axis.text.x= element_blank(),
              axis.ticks.x = element_blank())


ggsave("../fig1-C.pdf") 

