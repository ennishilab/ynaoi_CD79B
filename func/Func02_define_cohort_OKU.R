# --------------------- #
library(dplyr)
# --------------------- #
.myfunc.env = new.env()
sys.source( "./func/OKU_meta.R", envir = .myfunc.env )
attach( .myfunc.env )

define_cohort <- 
  function(){
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n ====== Start OKU cohort Difinition ====== \n")
    
    cat("\n=== Shaping Data ===\n")
    df <- read.csv("./saved/Scores/OKU_Final.csv", row.names = 1)
    df <- df[df$available_clinicaldata=="Available" &
               df$available_IHCdata=="Available", ]
    cat("\n---- both clinical and IHC(CD79B/CD20/CD79A) data are available in",
        length(unique(df$Patient_ID)), "samples ----\n")
    
    cat("\n---- Tumor Content (cutoff=20%) ----\n")
    df$TumorContent <- factor(df$TumorContent, levels = c("Low","High"))
    print(table(df$TumorContent))
    cat("\n---- Samples with Low tumor content", 
        paste0("(n=",nrow(df[df$TumorContent=="Low",]),")"), "are excluded ----\n")
    df <- df[df$TumorContent!="Low", ]
    
    cat("\n---- Combine Staining Pattern Data ----\n")
    pattern <- read.csv("./src/OKU/CD79B_staining_pattern_Final.csv")
    pattern$Pattern <- gsub("M ","M", pattern$Pattern)
    pattern$Dot.Pattern <- gsub("N ","N", pattern$Dot.Pattern)
    pattern$Pattern <- factor(pattern$Pattern, levels = c("M","MC","C","NEG","UD"))
    pattern$Dot.Pattern <- factor(pattern$Dot.Pattern, levels = c("N","Y","NEG","UD"))
    pattern$available_CD79B_Pattern <- rep("Available")
    pattern <- pattern[,c(1,2,6,3,4,5)]
    pattern <- pattern[,-6]
    colnames(pattern) <- c("Patient_ID", "core_ID", "available_CD79B_Pattern_OKU", 
                           "OKU_CD79b_Pattern", "OKU_CD79b_PosRate_Manual")
    df <- left_join(df, pattern, by=c("Patient_ID","core_ID"))
    
    df <- mutate_at(df, "available_CD79B_Pattern_OKU", ~replace(., is.na(.), "Unavailable"))
    cat("\n---- Pattern data is unavailable in n =", length(unique(df[df$available_CD79B_Pattern_OKU=="Unavailable",]$Patient_ID)),
        "samples, so exclude ----\n")
    df <- df[df$available_CD79B_Pattern_OKU=="Available",]
    print(table(df$OKU_CD79b_Pattern))
    cat("\n---- n =", length(unique(df[df$OKU_CD79b_Pattern=="UD",]$Patient_ID)),
        "are UD pattern, so exclude ----\n")
    df <- df[df$OKU_CD79b_Pattern!="UD",]
    
    df$OKU_CD79b_Pattern <- factor(df$OKU_CD79b_Pattern, levels = c("NEG","C","MC","M"))
    table(df$OKU_CD79b_Pattern)
    cat("\n---- n =", length(unique(df[df$OKU_CD79b_Pattern=="MC",]$Patient_ID)),
        "are MC pattern, so change to M pattern ----\n")
    df$OKU_CD79b_Pattern <- gsub("MC", "M", df$OKU_CD79b_Pattern)
    df$OKU_CD79b_Pattern <- factor(df$OKU_CD79b_Pattern, levels = c("NEG","C","M"))
    table(df$OKU_CD79b_Pattern)
    cat("\n---- n =", length(unique(df$Patient_ID)), " are the Final cohort (with UNCLASS) ----\n")
    
    res <- summary(df$OKU_CD79b_Hscore)
    print(round(res, 1))
    cat("\n---- Define CD79B Hscore class using median Hscore", 
        paste0("(", round(res[[3]],1),")"), "as threshold ----\n")
    df$OKU_CD79b_Hscore_class <- ifelse(df$OKU_CD79b_Hscore>=res[[3]], "High", "Low")
    df$OKU_CD79b_Hscore_class <- factor(df$OKU_CD79b_Hscore_class, levels = c("Low","High"))
    print(table(df$OKU_CD79b_Hscore_class))
    
    cat("\n---- Define CD79B positive/negative using Staining pattern ----\n")
    cat("\n---- i.e. M,C:positive, NEG:negative ----\n")
    df$OKU_CD79b_Pattern_PosNeg <- ifelse(df$OKU_CD79b_Pattern=="NEG","NEG","POS")
    df$OKU_CD79b_Pattern_PosNeg <- factor(df$OKU_CD79b_Pattern_PosNeg, levels = c("NEG","POS"))
    print(table(df$OKU_CD79b_Pattern))
    print(round(prop.table(table(df$OKU_CD79b_Pattern))*100,1))
    print(table(df$OKU_CD79b_Pattern_PosNeg))
    print(round(prop.table(table(df$OKU_CD79b_Pattern_PosNeg))*100,1))
    
    cat("\n---- Adding other IHC data ----\n")
    meta <- OKU_meta()
    meta <- meta[["meta0"]]
    meta <- meta %>% dplyr::rename("Patient_ID"="OHD_ID")
    other_IHC <- c("Patient_ID", 
                   "CD10.Positive", "MUM1.Positive", "BCL6.Positive",
                   "CD3.Positive", "CD4.Positive", "CD8.Positive", 
                   "CD31.Positive", 
                   "CD68.Positive", "CD68.H.score", 
                   "BCL2.Positive", "c.Myc.Positive", "MYC_FISH", "BCL2_FISH"
                   )
    meta <- meta[,other_IHC]
    colnames(meta) <- c("Patient_ID", 
                        "OKU_CD10_PosRate", "OKU_MUM1_PosRate", "OKU_BCL6_PosRate",
                        "OKU_CD3_PosRate", "OKU_CD4_PosRate", "OKU_CD8_PosRate", 
                        "OKU_CD31_PosRate", 
                        "OKU_CD68_PosRate", "OKU_CD68_Hscore", 
                        "OKU_BCL2_PosRate", "OKU_cMyc_PosRate", "OKU_MYC_FISH", "OKU_BCL2_FISH"
                        )
    df <- left_join(df, meta, by="Patient_ID")
    
    cat("\n---- Define Hans classification ----\n")
    df$Hans <- ifelse(df$OKU_CD10_PosRate>=0.3, "GCB",
                      ifelse(df$OKU_CD10_PosRate<0.3 & df$OKU_BCL6_PosRate<0.3, "NonGC",
                             ifelse(df$OKU_CD10_PosRate<0.3 & df$OKU_BCL6_PosRate>0.3 & df$OKU_MUM1_PosRate>0.3, "NonGC", "GCB")))
    df$Hans <- factor(df$Hans, levels = c("NonGC", "GCB"))
    print(table(df$Hans, df$DLBCL90))
    print(table(df$Hans))
    print(round(prop.table(table(df$Hans)),2))
    
    cat("\n---- Adding MHC data ----\n")
    mhc1 <- read.csv("./saved/Scores/OKU_MHC1_Scores_ID.csv")
    mhc2 <- read.csv("./saved/Scores/OKU_MHC2_Scores_ID.csv")
    df <- left_join(df, mhc1, by = c("Patient_ID","core_ID"))
    df <- left_join(df, mhc2, by = c("Patient_ID","core_ID"))
    
    write.csv(df, "./saved/Scores/OKU_defined_cohort_withUNC.csv", row.names = F)
    
    # COO=UNCLASS -> exclude
    print(table(df$DLBCL90))
    print(round(prop.table(table(df$DLBCL90))*100,2))
    cat("\n---- n =", nrow(df[df$DLBCL90=="UNCLASS",]),"are UNCLASS COO, so exclude ----\n")
    df <- df[df$DLBCL90 != "UNCLASS", ]
    df$DLBCL90 <- factor(df$DLBCL90, levels = c("ABC","GCB","DZsig-pos"))
    cat("\n---- n =", length(unique(df$Patient_ID)), " are the Final cohort (no UNCLASS) ----\n")
    print(table(df$DLBCL90))
    print(round(prop.table(table(df$DLBCL90))*100,2))
    
    print(table(df$DLBCL90, df$OKU_CD79b_Pattern))
    print(table(df$DLBCL90, df$OKU_CD79b_Pattern) %>% chisq.test())
    print(table(df$DLBCL90, df$OKU_CD79b_Pattern) %>% fisher.test(workspace = 1e+09))
    
    write.csv(df, "./saved/Scores/OKU_defined_cohort_withoutUNC.csv", row.names = F)
    cat("\n ====== Finish OKU cohort Difinition ====== \n")
    }































