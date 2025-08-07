# Library ----
# --------------------- #
library(dplyr)
# --------------------- #

# rLN GC markers----
Public_GCmarkers <- function(){
  jem_20200483_tables2 <- read.csv("./src/scRNAseq/JEM2020/article/jem_20200483_tables2.csv", skip = 2)[,1:26]
  jem_20200483_tables2[jem_20200483_tables2 == ""] <- NA
  
  signatures <- colnames(jem_20200483_tables2) %>%
    gsub("[.]", "_" , .) %>%
    sub("_", "", .)
  
  markers_UP <- list()
  for(i in seq(1,length(signatures),2)){
    D <- jem_20200483_tables2[-1,][,i]
    markers_UP[[i]] <- D[is.na(D)==F]
    names(markers_UP)[i] <- signatures[i]
    markers_UP <- markers_UP[!sapply(markers_UP, is.null)]
  }
  
  names(markers_UP)[11] <- "PreM_UP"
  
  saveRDS(markers_UP, "./src/scRNAseq/JEM2020/article/markers_UP.rds")
  
  # signatures <- reshape2::melt(markers_UP)
  # signatures <- signatures[,c(2,1)]
  # colnames(signatures) <- c("TERM","GENE")
  # dput(unique(signatures$TERM))
  # signatures$TERM <- factor(signatures$TERM,
  #                           levels = c("DZa_UP", "DZb_UP", "DZc_UP", 
  #                                      "INTa_UP", "INTb_UP", "INTc_UP","INTd_UP", "INTe_UP", 
  #                                      "LZa_UP", "LZb_UP", 
  #                                      "PreM_UP", 
  #                                      "PBLa_UP", "PBLb_UP"))
  # signatures <- signatures[order(signatures$TERM),]
  # signatures -> GC.gene.sig.detail
  # 
  # return(GC.gene.sig.detail)
}
