#----------------------------#

setwd("/Users/maynefan/Desktop/Project/CRC")

library(dplyr)
library(tidyverse)
library(readxl)

#--------- metabolite bioproject
metabolite <- read_xlsx("PRJDB4176_tableS13_代谢组_有饮酒信息.xlsx",sheet = 18) %>%
  slice(-c(1,2)) 
colnames(metabolite) <- metabolite[1,]
metabolite <- metabolite[-1,]
names(metabolite)[1] <- "Metabolite"

metabolite_sample <- colnames(metabolite) %>%
  as.data.frame() %>%
  slice(-1) %>%
  setNames("Sample_name")

#--------- metadata bioproject
metadata <- read.table("SraRunTable.csv",sep = ",",header = T) %>%
  mutate(Sample_name = as.character(Sample_name)) %>%
  right_join(metabolite_sample , by = "Sample_name")
  
  
#--------- metadata gwas

