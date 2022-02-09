paks<-c("dplyr", "ggplot2", "bio3d", "IRkernel", "pheatmap", "jsonlite", "reshape2")

install.packages("BiocManager")

BiocManager::install(version = "3.14")
BiocManager::install(paks)

IRkernel::installspec()