##Script to Figure 4-A.

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("ggraph","igraph","graphics","dplyr","RColorBrewer","ggplot2")
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
HGT_viruse_content_chayi2 <- read.table("~/cooperation/202409zhaofengxiang/data/fig4A_HGT_viruse_content_chayi2")
HOST_bact_ABUNDANCE2 <-  read.table("~/cooperation/202409zhaofengxiang/data/fig4A_HostBact_Abundance2")
DR_paintdata <- read.table("~/cooperation/202409zhaofengxiang/data/fig4A_HostBact_Donor_Recipient")


############################## Data processing ##############################
HGT_viruse_content_chayi2_forp <- HGT_viruse_content_chayi2[grep("D0|D10",HGT_viruse_content_chayi2$GD),]
HGT_viruse_content_chayi_lineage <- data.frame(HGT_viruse_content_chayi_lineage=c(HGT_viruse_content_chayi2_forp$recipient_lineage_genus,HGT_viruse_content_chayi2_forp$donor_lineage_genus))
HGT_viruse_content_chayi_lineage <- unique(HGT_viruse_content_chayi_lineage )
HGT_viruse_content_chayi_lineage$Kingdom <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",1)
HGT_viruse_content_chayi_lineage$phylum <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",2)
HGT_viruse_content_chayi_lineage$class <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",3)
HGT_viruse_content_chayi_lineage$order <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",4)
HGT_viruse_content_chayi_lineage$family <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",5)
HGT_viruse_content_chayi_lineage$genus <- sapply(strsplit(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage,split = ";"), "[",6)
HGT_viruse_content_chayi_lineage <- HGT_viruse_content_chayi_lineage[order(HGT_viruse_content_chayi_lineage$HGT_viruse_content_chayi_lineage),]
HGT_viruse_content_chayi_lineage2 <- merge(HGT_viruse_content_chayi_lineage,HOST_bact_ABUNDANCE2,by.x="genus",by.y="bact_g",all.x=T) 

unique_rank<- unique(HGT_viruse_content_chayi_lineage[c("Kingdom","phylum")])
unique_rank <- unique_rank[order(unique_rank[,2]),]
a <- aggregate(list(abundance=HGT_viruse_content_chayi_lineage2$abundance),by=list(HGT_viruse_content_chayi_lineage2[,c("phylum")]),sum )
a <- a[order(a$Group.1) ,]
edges <- data.frame(from=unique_rank[,1],to=unique_rank[,2],abundance=a$abundance)
unique_rank<- unique(HGT_viruse_content_chayi_lineage[c("phylum","class")])
a <- aggregate(list(abundance=HGT_viruse_content_chayi_lineage2$abundance),by=list(HGT_viruse_content_chayi_lineage2[,c("class")]),sum )
unique_rank <- unique_rank[order(unique_rank[,2]),]
a <- a[order(a$Group.1) ,]
edges <- rbind(edges,data.frame(from=unique_rank[,1],to=unique_rank[,2], abundance=a$abundance))
unique_rank<- unique(HGT_viruse_content_chayi_lineage[c( "class","order" )])
a <- aggregate(list(abundance=HGT_viruse_content_chayi_lineage2$abundance),by=list(HGT_viruse_content_chayi_lineage2[,c("order" )]),sum )
unique_rank <- unique_rank[order(unique_rank[,2]),]
a <- a[order(a$Group.1) ,]
edges <- rbind(edges,data.frame(from=unique_rank[,1],to=unique_rank[,2], abundance=a$abundance))
unique_rank<- unique(HGT_viruse_content_chayi_lineage[c("order", "family")])
a <- aggregate(list(abundance=HGT_viruse_content_chayi_lineage2$abundance),by=list(HGT_viruse_content_chayi_lineage2[,c("family")]),sum )
unique_rank <- unique_rank[order(unique_rank[,2]),]
a <- a[order(a$Group.1) ,]
edges <- rbind(edges,data.frame(from=unique_rank[,1],to=unique_rank[,2], abundance=a$abundance))
unique_rank<- unique(HGT_viruse_content_chayi_lineage[c("family" ,"genus")])
a <- aggregate(list(abundance=HGT_viruse_content_chayi_lineage2$abundance),by=list(HGT_viruse_content_chayi_lineage2[,c("genus")]),sum )
unique_rank <- unique_rank[order(unique_rank[,2]),]
a <- a[order(a$Group.1) ,]
edges <- rbind(edges,data.frame(from=unique_rank[,1],to=unique_rank[,2], abundance=a$abundance))

edges$from <-sapply(strsplit(edges$from,split = "__",fixed = T),"[",2)
edges$to <-sapply(strsplit(edges$to,split = "__",fixed = T),"[",2)

nodes <- data.frame(row.names =unique(c(edges$from,edges$to)) ,
                    name=unique(c(edges$from,edges$to))
)
nodes <- merge(nodes,edges,by.x="name",by.y="to",all=T)
colnames(nodes)[2] <- "group"

net<-graph_from_data_frame(d=edges,vertices = nodes,directed =T)


############################## [Fig4-A] Painting #####################################
###树状图+HGT事件方向
jiantous <- HGT_viruse_content_chayi2_forp[c("donor_bin","recipient_bin","recipient_lineage_genus","donor_lineage_genus","donor_genus","recipient_genus","GD","Diet","Group")]
jiantous <- unique(jiantous)
HGT_viruse_content_chayi_lineage$donor_number <- 0:(nrow(HGT_viruse_content_chayi_lineage)-1)
jiantous <- merge(jiantous,HGT_viruse_content_chayi_lineage[c("genus","donor_number")],by.x="donor_genus",by.y="genus")
HGT_viruse_content_chayi_lineage$recipient_number <- 0:(nrow(HGT_viruse_content_chayi_lineage)-1)
jiantous <- merge(jiantous,HGT_viruse_content_chayi_lineage[c("genus","recipient_number")],by.x="recipient_genus",by.y="genus")
jiantous <- jiantous[order(jiantous$donor_genus,jiantous$recipient_genus),]
jiantous$color[jiantous$GD=="D0"] <- "gray"
jiantous$color[jiantous$GD=="D10"] <- as.character(factor(jiantous$Diet[jiantous$GD=="D10"],levels = c("CON","HFD","FUC"),labels = c(NA,"#FB8073","#B2DE69") ))
jiantous$color[jiantous$GD=="D14"] <- as.character(factor(jiantous$Diet[jiantous$GD=="D14"],levels = c("CON","HFD","FUC"),labels = c(NA,"#FB8073","#B2DE69") ))
jiantous$color[jiantous$GD=="D18"] <- as.character(factor(jiantous$Diet[jiantous$GD=="D18"],levels = c("CON","HFD","FUC"),labels = c(NA,"#FB8073","#B2DE69") ))
df <- data.frame(x=c(),y=c())
df_names <- paste("df",1:nrow(jiantous),sep = "_")

p <- ggraph(net,layout = "dendrogram")+
  # geom_edge_elbow() +
  theme_graph()+
  # geom_edge_diagonal(aes(color=..index..))+
  # scale_edge_color_distiller(palette = "PuBu")+
  geom_edge_diagonal(color="#66C2A5")+
  geom_node_text(aes(label=name,x=x-0.25),angle=270,size=2.5,color="black")+
  # geom_node_point(color="#66C2A5")
  geom_node_point(aes(size=log10(abundance+1)),color="#66C2A5"  )+
  scale_size_continuous(range = c(0.8,3.5) )
for (i in 1:nrow(jiantous)){
  df_name <- df_names[i]
  color=jiantous$color[i]
  df2 <- data.frame(x=c(jiantous$donor_number[i],jiantous$donor_number[i],jiantous$recipient_number[i],jiantous$recipient_number[i]   ),
                    y=c(-0.5-0.1*i,-0.5-0.1*i-0.25,-0.5-0.1*i-0.25,-0.5-0.1*i))
  # df2 <- data.frame(x=c(jiantous$donor_number[i],jiantous$donor_number[i],jiantous$recipient_number[i],jiantous$recipient_number[i]   ),
  #                   y=c(0,0-0.2*i-0.25,0-0.2*i-0.25,0))
  assign(df_name,df2)
  p <- p+geom_path(data = get(df_name),aes(x,y),arrow = arrow(length  = unit(0.01, "native")),
                   color=color)
}
p


###供体受体数量
p2 <- stars(log10(DR_paintdata[2:3]+1), key.loc = c(0,3), scale = TRUE, 
      locations = NULL, len =1, radius = TRUE,
      full = TRUE, labels = NULL,draw.segments = TRUE,nrow = 1,
      col.segments=c("#BDB9B8","#CD776D"))

