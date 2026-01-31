##Script to Supplementary fig2-E

############################## Packages loading ############################## 
#install if they are not installed yet)
cran_packages=c("eulerr")
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

############################## Data processing ##############################
process_data <- function(data, diet) {
  data[is.na(data)] <- 0
  data <- data[, colSums(data == diet) != 0, drop = FALSE]
  data <- data[rowSums(data == diet) != 0, , drop = FALSE]
  data <- data[grep("s__", rownames(data)), , drop = FALSE]
  rownames(data)
}

Diet1[is.na(Diet1)] <- 0
high_fiber <- process_data(Diet1, "High_fiber_diet")
vegetarian <- process_data(Diet1, "Vegatarian")

venn_data <- list(High_fiber = high_fiber, Vegetarian = vegetarian)
################################ [FigS2-E] Painting ################################
library(eulerr)
plot(euler(venn_data), 
     fills = list(fill = c("#DDFFDD","#B8E5FA","#DDDDFF"), 
                  alpha = 0.6),
     quantities = TRUE, 
     main = "Species overlap between diets")


