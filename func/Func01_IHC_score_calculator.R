# --------------------- #
library(data.table)
library(pbapply)
library(stringr)
library(doParallel)
library(tidyr)
library(dplyr)
library(tictoc)
library(openxlsx)
library(reshape2)
# --------------------- #

Calculate_Scores <- 
  function(target, pattern, Detect.prob, multicore, cohort){
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("=== Shaping Data ===")
    
    files <- list.files(path = paste0("../QuPath/", target, "/measure"), 
                        pattern = "detect.txt",
                        full.names = TRUE)
    
    exst_cols <- c("Image","TMA core", "Class", "DAB: Mean", "Detection probability", 
                   "Area µm^2", "Circularity", "Solidity", "Centroid X µm", "Centroid Y µm")
    
    
    tictoc::tic()
    
    cl <- makeCluster(multicore) 
    registerDoParallel(cl)
    
    dat <- foreach(i=1:length(files), .combine=rbind, .packages="foreach") %dopar% {
      data.table::fread(file = files[i], encoding = "UTF-8", sep="\t",
                        select = exst_cols)
      
    }
    
    stopCluster(cl)
    
    cat("\n")
    tictoc::toc() 
    
    dat <- data.frame(dat)
    dim(dat) 
    cat("\n- Detection probability before filtering\n"); print(summary(dat$Detection.probability))
    Detect.prob <- 0.05
    dat <- dat[dat$Detection.probability > Detect.prob, ]
    cat("\n- Detection probability after filtering\n"); print(summary(dat$Detection.probability))
    dim(dat) 
    
    if(cohort=="OKU"){
      cat("\n- cohort =", cohort, "-")
      dat <- data.frame(dat)
      dat$core_ID <- paste0("No.",str_extract(dat$Image, "0\\d|1\\d"))
      dat$core_ID <- paste0(paste(dat$core_ID, dat$TMA.core, sep = "_"))
    }
    
    if(cohort=="BCC"){
      cat("\n- cohort =", cohort, "-")
      dat$core_ID <- dat$Image
      dat <- dat %>% separate(core_ID, c("core_ID","x"), sep="_Sector")
      dat$x <- NULL
      check <- data.frame(unique(dat$core_ID))
      cat("\n- n=",nrow(check),"unique core_IDs -\n") #[1] "n=3unique core_IDs"
      
      dat$core_ID <- paste0(dat$Image, "_", dat$TMA.core)
      check <- data.frame(unique(dat$core_ID))
      cat("\n- n=",nrow(check),"unique cell_IDs -\n") #[1] "n=301unique cell_IDs"
    }
    
    dat <- dat %>%
      arrange(core_ID) %>%
      group_by(core_ID) %>%
      mutate(cell_ID = row_number())
    
    dat$cell_ID <- formatC(dat$cell_ID, width = max(nchar(dat$cell_ID)), flag = "0")
    
    dat$cell_ID <- paste0(dat$core_ID, "_", dat$cell_ID)
    exst <- c("cell_ID", "core_ID", "core_ID", "Centroid.X.µm", "Centroid.Y.µm", "Detection.probability", 
              "Area.µm.2", "Circularity", "Solidity", "DAB..Mean", "Class")
    
    dat.exst <- dat[ ,exst]
    dat.exst <- data.frame(dat.exst)
    colnames(dat.exst) <- c("cell_ID","core_ID", "core_ID", "centroid_x", "centroid_y", 
                            "detection_prob", "area", "circularity", 
                            "solidity", "DAB_mean", "DAB_class")
    
    rownames(dat.exst) <- dat.exst$cell_ID
    dat.exst <- dat.exst[dat.exst$DAB_class %in% c("1+","2+","3+","Negative"),]
    
    dim(dat.exst) 
    
    if(cohort=="OKU"){
      dat.exst <- dat.exst[nchar(dat.exst$core_ID) == 9, ]
      dim(dat.exst)
    }
    
    cat("\n == Slides: n =", length(unique(dat.exst$core_ID)), "==")
    cat("\n == TMA cores: n =", length(unique(dat.exst$core_ID)), "==")
    
    Scores <- data.frame(prop.table(table(dat.exst$core_ID, dat.exst$DAB_class), margin = 1))
    Scores <- Scores[order(Scores$Var1), ]
    Scores <- spread(data = Scores, key = Var2, value = Freq)
    
    Scores <- Scores %>% rename("Patient_ID" = "Var1")
    bind <- Scores[,1:2]
    bind$`1+` <- NULL
    
    Scores$n_total <- data.frame(table(dat.exst$core_ID))[,2]
    
    n_pos <- data.frame(table(dat.exst[dat.exst$DAB_class != "Negative", ]$core_ID))
    colnames(n_pos) <- c("Patient_ID", "n_pos")
    bind_pos <- left_join(bind, n_pos, by = "Patient_ID")
    bind_pos <- mutate_at(bind_pos, "n_pos", ~replace(., is.na(.), 0))
    
    Scores <- left_join(Scores, bind_pos, by = "Patient_ID")
    Scores$PosRate <- Scores$n_pos / Scores$n_total
    
    rownames(Scores) <- Scores$Patient_ID
    Scores$Patient_ID <- NULL
    
    Scores$Hscore <- apply(X = Scores, MARGIN = 1, FUN = function(x){
      (x[["Negative"]]*100)*0 + (x[["1+"]]*100)*1 + (x[["2+"]]*100)*2 + (x[["3+"]]*100)*3
    })
    
    cat("\n == Missing Hscore == "); print(table(is.na(Scores$Hscore)))
    cat("\n == Missing PosRate == "); print(table(is.na(Scores$PosRate)))
    
    Scores$core_ID <- rownames(Scores)
    
    write.csv(x = Scores, file = paste0("./saved/Scores/", target, "_Scores.csv"), row.names = F)
    
    tictoc::toc()
  }

Arrange_Files <-
  function(target, minimum_counts, cohort){
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== Arranging File ===\n")
    # read Hscore file
    Scores <- read.csv(paste0("./saved/Scores/", target, "_Scores.csv"))
    cat(paste("\n === Before arranging,", nrow(Scores), "TMA cores are included === \n"))
    
    if(cohort=="BCC"){
      Scores$core_ID <- gsub(paste0("_Sector","\\d"),"",Scores$core_ID)
      Scores$core_ID <- gsub(paste0("_0","\\d"),"",Scores$core_ID)
    }
    
    # convert slide number to patient unique IDs
    if(cohort=="OKU"){
      id_converter <- read.csv("./src/OKU/TMA_OHDID.csv")
      id_converter <- id_converter %>% dplyr::rename("Patient_ID" = "OHD_ID")
    }
    
    if(cohort=="BCC"){
      id_converter <- data.frame(fread("./src/BCC/BCC_TMA_DLC_to_Core.csv"))
      
      id_converter <- id_converter %>% unite(core_ID, Image, Core, remove = F, sep = "_")
      id_converter$core_ID <- paste0("DLBCL380_", id_converter$core_ID)
      id_converter <- id_converter %>% dplyr::rename("Patient_ID" = "DLC")
    }
    
    cat(paste("\n ==", length(setdiff(Scores$core_ID, id_converter$core_ID)), "missing core_IDs == \n"))
    
    foo <- left_join(Scores, id_converter, by = "core_ID")
    foo <- foo[,c("Patient_ID", "core_ID", "Hscore", "PosRate", "n_total")]
    colnames(foo) <- c("Patient_ID", "core_ID", paste0(target,"_Hscore"),paste0(target,"_PosRate"),paste0(target,"_TotalCells"))
    cat(paste("\n ==", length(unique(foo$Patient_ID)),  "unique samples in total == \n"))
    
    foo <- foo[foo[,paste0(target,"_TotalCells")]>=minimum_counts, ]
    cat(paste0("\n == ", length(unique(foo$Patient_ID)),  " samples have sufficient cells (>=", minimum_counts, ") == \n"))
    cat(paste0("\n == ref. Rimsza LM, et al. Blood. 2004. == \n"))
    
    if(nrow(foo[is.na(foo$Patient_ID),])>0){
      cat(paste("\n == ==", nrow(foo[is.na(foo$Patient_ID),]),  " Patient_IDs are NA, so remove == == \n"))
      foo <- foo[!is.na(foo$Patient_ID),]
      cat(paste0("\n == == After removing, n=", length(unique(foo$Patient_ID)),  " unique Patient_IDs == == \n"))
    }
    
    duplicated <- foo[duplicated(foo$Patient_ID)|duplicated(foo$Patient_ID, fromLast = TRUE),] 
    if(nrow(duplicated)>0){
      duplicated <- duplicated %>% arrange(Patient_ID)
      cat(paste("\n == ==", nrow(duplicated),  "IDs are duplicated == == \n"))
      
      foo <- foo[!(foo$Patient_ID %in% duplicated$Patient_ID), ]
      
      df<- data.frame(matrix(ncol = ncol(duplicated),nrow = nrow(duplicated)))
      colnames(df) <- colnames(duplicated)
      
      for (i in seq(1,nrow(duplicated),2)){
        # for duplicated cases, use maximum Hscore
        compare <- c(duplicated[i,paste0(target,"_Hscore")],duplicated[i+1,paste0(target,"_Hscore")])
        Maxi <- max(compare)
        MaxLocation <- which.max(compare)
        
        if (MaxLocation == 1){
          df[i,] <- duplicated[i,]
        }
        
        else{
          df[i,] <- duplicated[i+1,]
        }
      } 
      
      df <- na.omit(df)
      
      foo <- rbind(foo, df)
      
      cat(paste("\n == After removing duplicated IDs, n =", nrow(foo),  "IDs are unique == == \n"))
    }
    
    foo <- foo %>% arrange(Patient_ID)
    
    if(cohort=="OKU"){
      rr <- foo[grep("R",foo$Patient_ID, invert=F), ]
      cat(paste("\n ==", nrow(rr), "TMA cores are from relapse cases == \n"))
      
      denovo <- foo[grep("R",foo$Patient_ID, invert=T), ]
      cat(paste("\n ==", nrow(denovo), "TMA cores are from de novo cases == \n"))
      write.csv(denovo, paste0("./saved/Scores/", target, "_Scores_ID.csv"), row.names = FALSE)
    }
    
    if(cohort=="BCC"){
      write.csv(foo, paste0("./saved/Scores/", target, "_Scores_ID.csv"), row.names = FALSE)
      }
    cat("\n === Finished merging Scoring data with Patient IDs === \n")
    }

Combine_metadata <- 
  function(cohort){
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n === Start adding clinical data === \n")
    
    if(cohort=="OKU"){
      cutoff = 0.2
      cat("\n=== current study cohort ===\n")
      meta <- read.csv("./src/OKU/OKU_metadata_final.csv")
      meta <- meta %>% dplyr::rename("Patient_ID" = "OHD_ID")
      df1 <- read.csv("./saved/Scores/OKU_CD79b_Scores_ID.csv")
      df2 <- read.csv("./saved/Scores/OKU_CD79a_Scores_ID.csv")
      df3 <- read.csv("./saved/Scores/OKU_CD20_Scores_ID.csv")
      df <- full_join(df1, df2, by=c("Patient_ID","core_ID"))
      df <- full_join(df, df3, by=c("Patient_ID","core_ID"))
      cat("\n == CD79B H-score data has,", nrow(df[is.na(df$OKU_CD79b_Hscore),]), "NA == \n")
      cat("\n == CD79A H-score data has,", nrow(df[is.na(df$OKU_CD79a_Hscore),]), "NA == \n")
      cat("\n == CD20 H-score data has,", nrow(df[is.na(df$OKU_CD20_Hscore),]), "NA == \n")
      
      tmp <- na.omit(df)
      tmp$available_IHCdata <- rep("Available")
      
      nas <- df[!(df$Patient_ID %in% tmp$Patient_ID), ]
      nas$available_IHCdata <- rep("Unavailable")
      
      df <- rbind(tmp, nas)
      df <- df %>% arrange(Patient_ID)
      df$TumorContent <- ifelse(df$OKU_CD20_PosRate>=cutoff|df$OKU_CD79a_PosRate>=cutoff, "High", "Low")
      cat("\n == To determine tumor content, ref. PMID: 23887293 ==\n")
      
      
      merged <- left_join(meta, df, by="Patient_ID")
      merged <- mutate_at(merged, 
                          c("available_clinicaldata",
                            "available_IHCdata"),
                          ~replace(., is.na(.),"Unavailable"))
      cat("\n == ", nrow(merged[merged$available_clinicaldata=="Available" & merged$available_IHCdata=="Available", ]), 
          "TMA cores are available both clinical and IHC data == \n")
    
      write.csv(merged, paste0("./saved/Scores/", cohort, "_Final.csv"))
    }    
    
    if(cohort=="BCC"){
      cat("\n=== validation cohort ===\n")
      meta <- read.xlsx("./src/BCC/meta_jco_updated.xlsx")
      meta <- meta %>% dplyr::rename("Patient_ID" = "DLC-ID")
      meta$available_clinicaldata <- ifelse(is.na(meta$DLBCL90_Group_updated), "Unavailable", "Available")
      Scores <- read.csv(paste0("./saved/Scores/BCC_CD79b_Scores_ID.csv"))
      Scores$Patient_ID <- gsub("_","",Scores$Patient_ID)
      Scores$available_IHCdata <- rep("Available")
      merged <- left_join(Scores, meta, by="Patient_ID")
      merged <- mutate_at(merged, 
                          c("available_clinicaldata",
                            "available_IHCdata"),
                          ~replace(., is.na(.),"Unavailable"))
      write.xlsx(merged, paste0("./saved/Scores/", cohort, "_Final.xlsx"))
      cat("\n == ", nrow(merged[merged$available_clinicaldata=="Available" & merged$available_IHCdata=="Available", ]), 
          "TMA cores are available both clinical and IHC data == \n")
    }
    
    cat("\n === Finished adding clinical data and saved file === \n")
    
    }




































