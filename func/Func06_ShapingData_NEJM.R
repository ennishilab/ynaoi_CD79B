# --------------------- #
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(edgeR)
library(DESeq2)
library(tibble)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GSVA)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(tidyHeatmap)
library(openxlsx)
# --------------------- #

ShapingData <- 
  function(exclude_UNC){
    cat("\n=== Start Extracting Protein-coding Genes ===\n")
    edb <- EnsDb.Hsapiens.v86
    gene.coding <- ensembldb::genes(edb, filter = ~ gene_biotype == "protein_coding")
    cat("\n---- length(gene.coding)", length(gene.coding), "genes ----\n")
    coding.symbol <- gene.coding$symbol
    
    cnt.tpm <- data.frame(data.table::fread(file = "./src/NEJM/GDC/TPM_GeneSymbol_220213.csv"), row.names = 1)
    colnames(cnt.tpm) <- gsub(pattern = "\\.", replacement = "_", colnames(cnt.tpm))
    
    cat("\n---- NEJM data has", 
        nrow(cnt.tpm), "genes,",
        ncol(cnt.tpm), " samples ----\n")
    
    cnt.tpm <- cnt.tpm[rownames(cnt.tpm) %in% coding.symbol, ]
    cat("\n---- extract protein coding genes", 
        nrow(cnt.tpm), "genes,",
        ncol(cnt.tpm), " samples ----\n")
    
    cnt.tpm <- data.frame(cnt.tpm)
    cnt.tpm.df <- data.frame(t(cnt.tpm))
    write.csv(cnt.tpm.df, "./saved/bulkRNAseq/NEJM_CodingGENE_NrmCount(TPM).csv")
    cat("\n=== Finished Extracting Protein-coding Genes ===\n")
    
    cat("\n=== Provided couint data is already normalized (method=TPM), so check it ===\n")
    cnt.tpm <- data.frame(data.table::fread("./saved/bulkRNAseq/NEJM_CodingGENE_RowCount.csv"), row.names = 1)
    cnt.tpm.t <- t(cnt.tpm)
    
    cat("\n---- Filtering ----\n")
    cnt.tpm$TotalExpr <- rowSums(cnt.tpm)
    cat("\n---- n =", nrow(cnt.tpm[cnt.tpm$TotalExpr==0, ]), " are zero total count, so exclude ----\n")
    cnt.tpm <- cnt.tpm[cnt.tpm$TotalExpr>0, ]
    cat("\n---- Filtered ----\n")
    
    cat("\n---- samples are randomly devidec into 5 clusters for normalization ----\n")
    
    nums <- sample(c(rep(seq(10),68), seq(3)))
    cnt.tpm$rand <- nums
    cnt.tpm$rand <- as.factor(cnt.tpm$rand)
    cnt.tpm$Patient_id <- rownames(cnt.tpm)
    rand_group <- cnt.tpm[,c("Patient_id","rand", "TotalExpr")]
    cnt.tpm$TotalExpr_log2 <- log2(cnt.tpm$TotalExpr)
    cnt.tpm$TotalExpr_log2_scaled <- scale(cnt.tpm$TotalExpr_log2)
    print(table(cnt.tpm$rand))
    
    tpm <- ggboxplot(cnt.tpm, 
                     x = "rand", 
                     y = "TotalExpr_log2_scaled", 
                     outlier.shape = F,
                     ylim = c(-5,1.5)
    ) + 
      ylab("log2 total expression")+
      xlab("randomly assigned group")+
      labs(title = "TPM-normalized")
    print(tpm)
    
    cat("\n=== Clinical metadata for NEJM cohort ===\n")
    meta <- data.frame(fread("./src/NEJM/GW Wright Cancer Cell 2020 Sup.Table2.csv", fill = T, encoding = "UTF-8", skip = 1, header = T))
    colnames(meta) <- gsub(pattern = "[.]", replacement = "_", colnames(meta))
    colnames(meta) <- gsub(pattern = "__", replacement = "_", colnames(meta))
    colnames(meta) <- gsub(pattern = "Years_", replacement = "Years", colnames(meta))
    meta <- meta[meta$Cohort=="NCI",]
    meta$Donor_name <- gsub("-","_",meta$Donor_name)
    print(table(meta$Dbl_Hit_Call, meta$COO_Class))
    
    meta$DLBCL90_Group_updated <- ifelse(meta$COO_Class=="GCB"&meta$Dbl_Hit_Call=="POS","DZsig-pos",
                                         ifelse(meta$COO_Class=="Unclass"&meta$Dbl_Hit_Call=="POS","DZsig-pos",
                                                ifelse(meta$COO_Class=="GCB"&meta$Dbl_Hit_Call=="NEG","GCB",
                                                       ifelse(meta$COO_Class=="Unclass"&meta$Dbl_Hit_Call=="NEG","UNC","ABC")))) 
    meta <- meta[!is.na(meta$DLBCL90_Group_updated),]
    cat("\n---- Clinical data is available for n =", 
        length(unique(meta$Donor_name)), "----\n")
    print(table(meta$Dbl_Hit_Call, meta$DLBCL90_Group_updated))
    meta$DLBCL90 <- factor(meta$DLBCL90_Group_updated, levels = c("DZsig-pos","GCB","ABC","UNC"))
    meta$R_CHOp_like_Chemo
    
    exst <- c("Donor_name", "Age_at_Diagnosis", "Gender", 
              "DLBCL90_Group_updated", "LymphGen_call",
              "IPI_Age", "IPI_ECOG", "IPI_LDH", "IPI_Stage", "IPI_Extranodal_Sites", 
              "IPI_Score", "R_CHOp_like_Chemo", "Time_to_last_Follow_up_Years", "Status_at_last_Follow_up")
    meta <- meta[,exst]
    cols <- c("Patient_ID", "AGE", "Sex", 
              "DLBCL90_Group_updated", "LymphGen_call",
              "IPI_Age", "IPI_ECOG", "IPI_LDH", "IPI_Stage", "IPI_Extranodal_Sites", 
              "IPI_Score", "RCHOP_like_Chemo", "OS", "CODE_OS")
    colnames(meta) <- cols
    
    cat("\n---- Shaping LymphGen cluster ----\n")
    
    cat("\n---- Adding CD79B Mutation data ----\n")
    mut <- data.frame(data.table::fread(file = "./src/NEJM/GDC/MAF_NCICCR-DLBCL_phs001444.txt", sep = "\t"))
    mut <- mut[mut$SUBJECT_NAME %in% meta$Patient_ID,]
    mut_CD79B <- mut[mut$GENE.SYMBOL == "CD79B",]
    cat("\n----", 
        paste0("n=",length(unique(mut_CD79B$SUBJECT_NAME))),
        paste0("out of n=", length(unique(mut$SUBJECT_NAME))),
        paste0("(",round((length(unique(mut_CD79B$SUBJECT_NAME)))/length(unique(mut$SUBJECT_NAME))*100,1),"%)"),
        "has CD79B mutation (only with metadata) ----\n")
    
    tmp <- meta[meta$Patient_ID %in% unique(mut$SUBJECT_NAME),]
    tmp1 <- tmp[tmp$Patient_ID %in% unique(mut_CD79B$SUBJECT_NAME),]
    tmp1$CD79B_M <- rep(1)
    tmp2 <- tmp[!(tmp$Patient_ID %in% mut_CD79B$SUBJECT_NAME),]
    tmp2$CD79B_M <- rep(0)
    tmp <- rbind(tmp1,tmp2)
    
    meta <- left_join(meta, tmp[,c("Patient_ID","CD79B_M")], by="Patient_ID")
    print(table(meta$CD79B_M))
    print(round(prop.table(table(meta$CD79B_M))*100,1))
    
    cat("\n---- Adding MYC Mutation data ----\n")
    mut_MYC <- mut[mut$GENE.SYMBOL == "MYC",]
    cat("\n----", 
        paste0("n=",length(unique(mut_MYC$SUBJECT_NAME))),
        paste0("out of n=", length(unique(mut$SUBJECT_NAME))),
        paste0("(",round((length(unique(mut_MYC$SUBJECT_NAME)))/length(unique(mut$SUBJECT_NAME))*100,1),"%)"),
        "has MYC mutation (only with metadata) ----\n")
    
    tmp <- meta[meta$Patient_ID %in% unique(mut$SUBJECT_NAME),]
    tmp1 <- tmp[tmp$Patient_ID %in% unique(mut_MYC$SUBJECT_NAME),]
    tmp1$MYC_M <- rep(1)
    tmp2 <- tmp[!(tmp$Patient_ID %in% mut_MYC$SUBJECT_NAME),]
    tmp2$MYC_M <- rep(0)
    tmp <- rbind(tmp1,tmp2)
    
    meta <- left_join(meta, tmp[,c("Patient_ID","MYC_M")], by="Patient_ID")
    print(table(meta$MYC_M))
    print(round(prop.table(table(meta$MYC_M))*100,1))
    
    cat("\n---- Shaping LymphGen Subtypes ----\n")
    meta$LymphGen_orig <- meta$LymphGen_call
    meta$LymphGen_call <- NULL
    meta$LymphGen_caption <- ifelse(grepl("/", meta$LymphGen_orig), "composite",
                                  ifelse(meta$LymphGen_orig == "Other", "Other", "primary"))
    meta$LymphGen_MCD_bin <- ifelse(grepl("MCD", meta$LymphGen_orig), "MCD", "Non-MCD")
    meta$LymphGen_MCD_tri <- ifelse(grepl("MCD", meta$LymphGen_orig) & meta$LymphGen_caption=="primary", "primary-MCD",
                                  ifelse(grepl("MCD", meta$LymphGen_orig) & meta$LymphGen_caption=="composite", "composite-MCD", "Non-MCD"))
    meta$LymphGen_MCD_bin <- factor(meta$LymphGen_MCD_bin, levels = c("Non-MCD","MCD"))
    
    print(table(meta[meta$LymphGen_caption=="primary",]$LymphGen_orig))
    print(table(meta$LymphGen_MCD_bin))
    print(round(prop.table(table(meta$LymphGen_MCD_bin))*100,1))
    print(table(meta$LymphGen_MCD_bin, meta$DLBCL90_Group_updated))
    
    tmp <- meta[meta$LymphGen_orig %in% c("EZB","EZB/ST2/A53","EZB/A53","EZB/ST2"), ]
    tmp$MYC_M_tmp <- tmp$MYC_M
    tmp$MYC_M_tmp <- gsub(1, "pos", tmp$MYC_M_tmp)
    tmp$MYC_M_tmp <- gsub(0, "neg", tmp$MYC_M_tmp)
    tmp$LymphGen_orig7 <- paste0("MYC", tmp$MYC_M_tmp, "EZB")
    tmp$MYC_M_tmp <- NULL
    
    meta <- meta[!(meta$Patient_ID%in%tmp$Patient_ID), ]
    meta$LymphGen_orig7 <- meta$LymphGen_orig
    
    meta <- rbind(meta, tmp)
    
    table(meta[meta$LymphGen_caption=="primary",]$LymphGen_orig7, 
          meta[meta$LymphGen_caption=="primary",]$DLBCL90_Group_updated)
    
    if(exclude_UNC=="Yes"){
      write.csv(meta[meta$DLBCL90_Group_updated %in% c("DZsig-pos","GCB","ABC"), ], 
                "./saved/bulkRNAseq/NEJM_metadata_withoutUNC.csv", row.names = F)
    }
    
    if(exclude_UNC=="No"){
      write.csv(meta, "./saved/bulkRNAseq/NEJM_metadata_withUNC.csv", row.names = F)
    }
    
    cat("\n=== Finished Shaping metadata ===\n")
  }

scGCcluster <-
  function(){
    cat("\n---- scGC-cluster from Holmes AB et al. JEM 2020 ----\n")
    cat("\n---- load GEP data (only protein coding genes) ----\n")
    
    meta <- read.csv("./saved/bulkRNAseq/NEJM_metadata_withUNC.csv")
    
    cnt.nrm <- data.frame(data.table::fread("./saved/bulkRNAseq/NEJM_CodingGENE_NrmCount(TPM).csv"), row.names = 1)
    cnt.nrm <- cnt.nrm[rownames(cnt.nrm) %in% meta$Patient_ID, ]
    cnt.nrm$Patient_ID <- rownames(cnt.nrm)
    
    holmes <- read.xlsx("./src/BCC/JEM2020/jem_20200483_tables5.xlsx", rowNames = F, startRow = 3, sheet = "Table_S5")
    exst <- c("Sample.ID","sc-COO.Class","sc-COO.Score","sc-COO.Group")
    holmes <- holmes[holmes$Dataset=="NCI-DLBCL", exst]
    colnames(holmes) <- c("Patient_ID","JEMsupp_scGC_Class","JEMsupp_scGC_Score","JEMsupp_scGC_Group")

    gep <- left_join(cnt.nrm, holmes, by="Patient_ID")
    
    gep$JEMsupp_scGC_Class_Simple <- 
      ifelse(gep$JEMsupp_scGC_Class %in% c("DZ a","DZ b","DZ c"), "DZ-like",
             ifelse(gep$JEMsupp_scGC_Class %in% c("LZ a","LZ b"), "LZ-like",
                    ifelse(gep$JEMsupp_scGC_Class %in% c("PBL a","PBL b"), "PBL-like",
                           ifelse(gep$JEMsupp_scGC_Class == "PreM", "PreM-like", "INT-like"))))
    gep$JEMsupp_scGC_Class_Simple <- factor(gep$JEMsupp_scGC_Class_Simple,
                                               levels = c("DZ-like","INT-like","LZ-like","PBL-like","PreM-like"))
    
    
    cat("\n---- In NEJM cohort, scGC cluster data is available in n =", 
        nrow(gep[!is.na(gep$JEMsupp_scGC_Class_Simple), ]), "samples (with COO=UNCLASS) ----\n")
    print(table(gep$JEMsupp_scGC_Class_Simple))
    print(round(prop.table(table(gep$JEMsupp_scGC_Class_Simple))*100,2))
    
    write.csv(gep, "./saved/bulkRNAseq/NEJM_scGCcluster.csv")
  }
































