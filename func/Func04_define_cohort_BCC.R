# --------------------- #
library(openxlsx)
library(dplyr)
# --------------------- #

define_cohort <- 
  function(exclude_UNC){
    exclude_UNC = exclude_UNC
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n ====== Start Analyizing BCC cohort ====== \n")
    
    cat("\n=== Shaping Data ===\n")
    df <- read.xlsx("./saved/Scores/BCC_Final.xlsx", rowNames = F)
    exst <- c("Patient_ID", "core_ID", "available_clinicaldata",
              "AGE", "Sex", "Stage", "Age_IPI", "PS_IPI", "LDH_IPI", 
              "Stage_IPI", "ES_IPI", "IPI", "CODE_OS", "CODE_PFS", 
              "OS", "PFS", "DLBCL90_Group_updated", 
              "BCC_CD79b_Hscore", "BCC_CD79b_PosRate", 
              "BCC_CD79b_TotalCells", "available_IHCdata", 
              "MYC-TR", "BCL2-TR", "BCL6-TR", "HGBL-DH/TH-BCL2", 
              "MYC-IHC", "BCL2-IHC", "DPE", 
              "CD10-IHC", "BCL6-IHC", "MUM1-IHC", "Hans_COO",
              "Ki67-IHC")
    df <- df[,exst]
    df <- df[df$available_clinicaldata=="Available" &
               df$available_IHCdata=="Available", ]
    cat("\n---- both clinical and IHC(CD79B) data are available in",
        length(unique(df$Patient_ID)), "samples ----\n")
    
    cat("\n---- Tumor Content is not calculated in BCCcohort ----\n")
    
    cat("\n---- Combine Staining Pattern Data ----\n")
    pattern <- read.csv("./src/BCC/BCC_CD79B_Staining_Pattern_Final_20240415.csv")
    pattern$Staining_pattern <- factor(pattern$Staining_pattern, levels = c("M","C","NEG","UD"))

    pattern$available_CD79B_Pattern <- rep("Available")
    pattern <- pattern[,c(1,2,11,18)]
    pattern <- pattern %>% dplyr::rename("BCC_CD79b_Pattern"="Staining_pattern")
    df <- left_join(df, pattern, by=c("core_ID","Patient_ID"))
    df <- mutate_at(df,
                    c("available_clinicaldata","available_CD79B_Pattern"),
                    ~replace(., is.na(.), "Unavailable"))
    print(table(df$BCC_CD79b_Pattern))
    # cat("\n---- n =", nrow(df[df$BCC_CD79b_Pattern=="UD",]),"are UD pattern, so exclude ----\n")
    # df <- df[df$available_CD79B_Pattern=="Available"&
    #            df$BCC_CD79b_Pattern!="UD", ]
    # df$BCC_CD79b_Pattern <- factor(df$BCC_CD79b_Pattern, levels = c("NEG","C","M"))
    
    cat("\n---- n =", length(unique(df$Patient_ID)), " are the Final cohort (with UNCLASS) ----\n")
    
    res <- summary(df$BCC_CD79b_Hscore)
    print(round(res, 1))
    cat("\n---- Define CD79B Hscore class using median Hscore", 
        paste0("(", round(res[[3]],1),")"), "as threshold ----\n")
    df$BCC_CD79b_Hscore_class <- ifelse(df$BCC_CD79b_Hscore>=res[[3]], "High", "Low")
    df$BCC_CD79b_Hscore_class <- factor(df$BCC_CD79b_Hscore_class, levels = c("Low","High"))
    print(table(df$BCC_CD79b_Hscore_class))
    
    cat("\n---- Define CD79B positive/negative using Staining pattern ----\n")
    cat("\n---- i.e. M,C:positive, NEG:negative ----\n")
    df$BCC_CD79b_Pattern_PosNeg <- ifelse(df$BCC_CD79b_Pattern=="NEG","NEG","POS")
    df$BCC_CD79b_Pattern_PosNeg <- factor(df$BCC_CD79b_Pattern_PosNeg, levels = c("NEG","POS"))
    print(table(df$BCC_CD79b_Pattern))
    print(round(prop.table(table(df$BCC_CD79b_Pattern))*100,1))
    print(table(df$BCC_CD79b_Pattern_PosNeg))
    print(round(prop.table(table(df$BCC_CD79b_Pattern_PosNeg))*100,1))
    
    if(exclude_UNC=="Yes"){
      df$DLBCL90_Group_updated <- factor(df$DLBCL90_Group_updated,
                                         levels = c("ABC","GCB","DZsig-pos","UNC"))
      print(table(df$DLBCL90_Group_updated))
      print(round(prop.table(table(df$DLBCL90_Group_updated))*100,2))
      cat("\n---- n =", nrow(df[df$DLBCL90_Group_updated=="UNC",]),"are UNCLASS COO, so exclude ----\n")
      df <- df[df$DLBCL90_Group_updated != "UNC", ]
      df$DLBCL90_Group_updated <- factor(df$DLBCL90_Group_updated, levels = c("ABC","GCB","DZsig-pos"))
      cat("\n---- n =", length(unique(df$Patient_ID)), " are the Final cohort (no UNCLASS) ----\n")
      print(table(df$DLBCL90_Group_updated))
      print(round(prop.table(table(df$DLBCL90_Group_updated))*100,2)) 
    }
    
    cat("\n---- Adding MHC IHC data ----\n")
    mhc <- read.xlsx("./src/BCC/CancerDiscovery2019/TableS1.xlsx", rowNames = F, startRow = 2)
    mhc <- mhc[,c("DLC_ID","MHC-I.expression","MHC-II.expression")]
    colnames(mhc) <- c("Patient_ID","MHCI_IHC","MHCII_IHC")
    df <- left_join(df, mhc, by = "Patient_ID")
    cat("\n MHCI \n")
    print(table(df$MHCI_IHC))
    cat("\n MHCII \n")
    print(table(df$MHCII_IHC))
    
    cat("\n---- Adding Mutation data ----\n")
    mut <- data.frame(fread("./src/BCC/DLC380_RNAseq/DLC_Mutation_20180313.csv"), row.names = 1)
    gene_order <- names(sort(colSums(mut), decreasing = T))
    donor_order <- names(sort(rowSums(mut), decreasing = T))
    mut$Patient_ID <- rownames(mut)
    mut$available_mut <- rep("Available")
    mut <- mut[,c(54,55,1:53)]
    genes <- colnames(mut[,3:55])
    
    df <- left_join(df, mut, by = "Patient_ID")
    df[,genes] <- lapply(df[,genes], as.factor)
    df <- mutate_at(df, 
                    c("available_mut"), 
                    ~replace(., is.na(.), "Unavailable"))
    print(table(df[!is.na(df$CD79B_M),]$CD79B_M))
    print(table(df[!is.na(df$CD79B_M),]$CD79B_M, 
                df[!is.na(df$CD79B_M),]$DLBCL90_Group_updated))
    
    cat("\n---- Adding LymphGen Subtypes ----\n")
    LG <- data.frame(fread("./src/BCC/GW_Wright_CancerCell2020SupTable2.csv", fill = T, encoding = "UTF-8", skip = 1, header = T))
    colnames(LG) <- gsub(pattern = "[.]", replacement = "_", colnames(LG))
    colnames(LG) <- gsub(pattern = "__", replacement = "_", colnames(LG))
    colnames(LG) <- gsub(pattern = "Years_", replacement = "Years", colnames(LG))
    LG <- LG %>% dplyr::rename("Patient_ID"="Donor_name")
    LG <- LG %>% dplyr::rename("LymphGen_orig"="LymphGen_call")
    LG$available_LymphGen <- rep("Available")
    LG$LymphGen_caption <- ifelse(grepl("/", LG$LymphGen_orig), "composite",
                                  ifelse(LG$LymphGen_orig == "Other", "Other", "primary"))
    LG$LymphGen_MCD_bin <- ifelse(grepl("MCD", LG$LymphGen_orig), "MCD", "Non-MCD")
    LG$LymphGen_MCD_tri <- ifelse(grepl("MCD", LG$LymphGen_orig) & LG$LymphGen_caption=="primary", "primary-MCD",
                                  ifelse(grepl("MCD", LG$LymphGen_orig) & LG$LymphGen_caption=="composite", "composite-MCD", "Non-MCD"))
    LG_BCC <- LG[LG$Cohort == "BCA", c("Patient_ID", "available_LymphGen","LymphGen_orig","LymphGen_caption",
                                       "LymphGen_MCD_bin","LymphGen_MCD_tri")]
    df <- left_join(df, LG_BCC, by = "Patient_ID")
    df <- mutate_at(df, 
                    c("available_LymphGen"), 
                    ~replace(., is.na(.), "Unavailable"))
    df$LymphGen_MCD_bin <- factor(df$LymphGen_MCD_bin, levels = c("Non-MCD","MCD"))
    
    print(table(df[df$LymphGen_caption=="primary",]$LymphGen_orig))
    print(table(df$LymphGen_MCD_bin))
    print(round(prop.table(table(df$LymphGen_MCD_bin))*100,1))
    print(table(df$LymphGen_MCD_bin, df$DLBCL90_Group_updated))
    
    tmp <- df[df$LymphGen_orig %in% c("EZB","EZB/ST2/A53","EZB/A53","EZB/ST2"), ]
    tmp$MYC_M_tmp <- tmp$MYC_M
    tmp$MYC_M_tmp <- gsub(1, "pos", tmp$MYC_M_tmp)
    tmp$MYC_M_tmp <- gsub(0, "neg", tmp$MYC_M_tmp)
    tmp$LymphGen_orig7 <- paste0("MYC", tmp$MYC_M_tmp, "EZB")
    tmp$MYC_M_tmp <- NULL
      
    df <- df[!(df$Patient_ID%in%tmp$Patient_ID), ]
    df$LymphGen_orig7 <- df$LymphGen_orig
    
    df <- rbind(df, tmp)
    
    table(df[df$LymphGen_caption=="primary",]$LymphGen_orig7, 
          df[df$LymphGen_caption=="primary",]$DLBCL90_Group_updated)
    
    if(exclude_UNC=="Yes"){
      write.csv(df, "./saved/Scores/BCC_defined_cohort_withoutUNC.csv", row.names = F)
    }
    
    if(exclude_UNC=="No"){
      write.csv(df, "./saved/Scores/BCC_defined_cohort_withUNC.csv", row.names = F)
    }
    
    }































