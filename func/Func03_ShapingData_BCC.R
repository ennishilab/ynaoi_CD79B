# --------------------- #
library(data.table)
library(dplyr)
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

Extract_CodingGenes <- 
  function(){
    cat("\n=== Start Extracting Protein-coding Genes ===\n")
    edb <- EnsDb.Hsapiens.v86
    gene.coding <- ensembldb::genes(edb, filter = ~ gene_biotype == "protein_coding")
    cat("\n---- length(gene.coding)", length(gene.coding), "genes ----\n")
    coding.id <- gene.coding$gene_id
    
    cnt.raw <- data.table::fread(file = "./src/BCC/DLC380_RNAseq/all_libs_raw_counts.txt", sep="\t")
    cnt.raw <- column_to_rownames(.data = cnt.raw, var = "gene_id")
    cnt.raw <- cnt.raw[, -1]
    cnt.raw <-  apply(X = cnt.raw, MARGIN = 2, FUN = function(x){ x * (10^6 / sum(x))})
    
    cat("\n---- BCC data has", 
        nrow(cnt.raw), "genes,",
        ncol(cnt.raw), " samples ----\n")
    
    cnt.raw <- cnt.raw[rownames(cnt.raw) %in% coding.id, ]
    cat("\n---- extract protein coding genes", 
        nrow(cnt.raw), "genes,",
        ncol(cnt.raw), " samples ----\n")
    
    Gene.name <- ensembldb::select(edb, rownames(cnt.raw), keytype = "GENEID", columns = c("GENEID","SYMBOL"))
    cnt.raw <- data.frame(cnt.raw)
    cnt.raw$GENEID <- rownames(cnt.raw)
    cnt.raw <- left_join(x = Gene.name, y = cnt.raw, by = "GENEID")
    cnt.raw <- cnt.raw[ !duplicated(cnt.raw$SYMBOL) , ]
    cnt.raw <- remove_rownames(cnt.raw) %>% column_to_rownames(var = "SYMBOL")
    cnt.raw$GENEID <- NULL
    
    cnt.raw.df <- data.frame(t(cnt.raw))
    
    write.csv(cnt.raw.df, "./saved/bulkRNAseq/BCC_CodingGENE_RowCount.csv")
    cat("\n=== Finished Extracting Protein-coding Genes ===\n")
  }

Normalize_GEP <- 
  function(){
    cat("\n=== Start Normalization ===\n")
    
    cnt.raw <- data.frame(data.table::fread("./saved/bulkRNAseq/BCC_CodingGENE_RowCount.csv"), row.names = 1)
    cnt.raw.t <- t(cnt.raw)
    
    cat("\n---- Filtering ----\n")
    cnt.raw$TotalExpr <- rowSums(cnt.raw)
    cat("\n---- n =", nrow(cnt.raw[cnt.raw$TotalExpr==0, ]), " are zero total count, so exclude ----\n")
    cnt.raw <- cnt.raw[cnt.raw$TotalExpr>0, ]
    cat("\n---- Filtered ----\n")
    
    cat("\n---- samples are randomly devidec into 5 clusters for normalization ----\n")
    
    nums <- sample(c(rep(seq(5),62), seq(4)))
    cnt.raw$rand <- nums
    cnt.raw$rand <- as.factor(cnt.raw$rand)
    cnt.raw$Patient_id <- rownames(cnt.raw)
    rand_group <- cnt.raw[,c("Patient_id","rand", "TotalExpr")]
    cnt.raw$TotalExpr_log2 <- log2(cnt.raw$TotalExpr)
    cnt.raw$TotalExpr_log2_scaled <- scale(cnt.raw$TotalExpr_log2)
    print(table(cnt.raw$rand))
    
    pre <- ggboxplot(cnt.raw, 
                     x = "rand", 
                     y = "TotalExpr_log2_scaled", 
                     outlier.shape = F,
                     ylim = c(-2,2)
    ) + 
      ylab("log2 total expression")+
      xlab("randomly assigned group")+
      labs(title = "pre-normalize")
    
    group <- data.frame(rand = factor(cnt.raw[ ,"rand"]))
    
    cat("\n---- Normalization (vst) ----\n")
    dds <- DESeqDataSetFromMatrix(countData = round(cnt.raw.t), colData = group, design = ~ rand)
    dds <- DESeq(dds)
    vsd <- vst(dds, blind=FALSE)
    
    cnt.nrm <- data.frame(t(assay(vsd)))
    cnt.nrm -> cnt.nrm.vst
    cnt.nrm$rand <- nums
    cnt.nrm$TotalExpr <- rowSums(cnt.nrm)
    cnt.nrm$TotalExpr_log2 <- log2(cnt.nrm$TotalExpr)
    cnt.nrm$TotalExpr_log2_scaled <- scale(cnt.nrm$TotalExpr_log2)
    
    post <- ggboxplot(cnt.nrm, 
                      x = "rand", 
                      y = "TotalExpr_log2_scaled", 
                      outlier.shape = F,
                      ylim = c(-2,2)
    ) + 
      ylab("log2 total expression")+
      xlab("randomly assigned group")+
      labs(title = "post-normalize (vst)")
    post -> post.vst
    
    cat("\n---- Normalization (EdgeR) ----\n")
    y <- DGEList(counts=cnt.raw.t, group=group$rand)
    y <- calcNormFactors(y, method = "TMM")
    cnt.nrm <- edgeR::cpm(y)
    cnt.nrm <- data.frame(t(cnt.nrm))
    cnt.nrm -> cnt.nrm.edgeR
    cnt.nrm$rand <- nums
    cnt.nrm$TotalExpr <- rowSums(cnt.nrm)
    cnt.nrm$TotalExpr_log2 <- log2(cnt.nrm$TotalExpr)
    cnt.nrm$TotalExpr_log2_scaled <- scale(cnt.nrm$TotalExpr_log2)
    
    post <- ggboxplot(cnt.nrm, 
                      x = "rand", 
                      y = "TotalExpr_log2_scaled", 
                      outlier.shape = F,
                      ylim = c(-2,2)
    ) + 
      ylab("log2 total expression")+
      xlab("randomly assigned group")+
      labs(title = "post-normalize (EdgeR)")
    post -> post.edger
    
    cat("\n==== visualize ====\n")
    plot(pre+post.vst+post.edger)
    
    write.csv(cnt.nrm.edgeR, "./saved/bulkRNAseq/BCC_CodingGENE_NrmCount(edgeR).csv")
    write.csv(cnt.nrm.vst, "./saved/bulkRNAseq/BCC_CodingGENE_NrmCount(vst).csv")
    cat("\n=== Finished Normalization ===\n")
  }

LME <-
  function(){
    cat("\n=== Calculate LME score (Kotolov N et al. Cancer Discovery. 2021.) ===\n")
    
    cat("\n---- load H-score data with clinical information (with UNCLASS) ----\n")
    hscore <- data.frame(data.table::fread("./saved/Scores/BCC_defined_cohort_withUNC.csv"))
    hscore <- hscore[hscore$available_clinicaldata=="Available"&hscore$available_IHCdata=="Available", ]
    ids <- hscore$Patient_ID
    
    cat("\n---- load GEP data (only protein coding genes) ----\n")
    cnt.nrm <- data.frame(data.table::fread("./saved/bulkRNAseq/BCC_CodingGENE_NrmCount(edgeR).csv"), row.names = 1)
    cnt.nrm <- cnt.nrm[rownames(cnt.nrm) %in% ids, ]
    
    cat("\n---- In this study cohort, GEP data is available in n =", nrow(cnt.nrm), "samples ----\n")
    cat("\n---- In n =", nrow(cnt.nrm), "samples, we're going to calculate LME scores ----\n")
    
    cat("\n---- Basically, LME score is GSVA score ----\n")
    
    cat("\n loading LME signatures ----\n")
    lme <- read.csv("./src/BCC/CancerDiscovery2021/CANCER_DISCOVERY_2021_LME_signatures.csv")
    LMEsig <- list()
    LMEsig <- apply(lme, MARGIN = 2, FUN = function(x){LMEsig <- c(LMEsig, x)})
    LMEsig <- lapply(LMEsig, FUN = function(x){x <- x[ x !=  ""] %>% unlist() })
    
    cat("\n---- calculate GSVA score as LME score ----\n")
    expr <- as.matrix(t(na.omit(cnt.nrm)))
    LMEscore <- gsva(expr = expr, 
                     gset.idx.list = LMEsig,
                     mx.diff = FALSE
                     )
    lme.targets <- rownames(LMEscore)
    LMEscore.df <- data.frame(t(LMEscore))
    LMEscore.df$Patient_ID <- rownames(LMEscore.df)
    
    exst <- c("DLBCL90_Group_updated", 
              "BCC_CD79b_Hscore", "BCC_CD79b_Hscore_class", 
              "BCC_CD79b_Pattern", "BCC_CD79b_Pattern_PosNeg")
    LMEscore.df <- left_join(LMEscore.df, hscore[,c("Patient_ID",exst)], by="Patient_ID")
    LMEscore.df$Patient_ID -> rownames(LMEscore.df)
    LMEscore.df$BCC_CD79b_Hscore_class <- factor(LMEscore.df$BCC_CD79b_Hscore_class, 
                                                 levels = c("Low","High"))
    
    LMEscore.df <- LMEscore.df %>% 
      arrange(BCC_CD79b_Pattern) %>%
      arrange(BCC_CD79b_Hscore_class)%>%
      arrange(BCC_CD79b_Pattern_PosNeg) %>%
      arrange(DLBCL90_Group_updated)
    anno <- LMEscore.df[,exst]
    mat <- LMEscore.df[,lme.targets]
    
    ph <- pheatmap::pheatmap(mat = t(mat), 
                             annotation_col = anno, 
                             show_colnames = F,
                             cluster_cols = T,
                             scale = "row"
    )
    print(ph)
    
    plot_cluster <- data.frame(x =matrix(nrow = nrow(LMEscore.df), ncol = 1))
    rownames(plot_cluster) <- rownames(LMEscore.df)
    for(i in 4){
      ph_cluster <- cutree(ph$tree_col, k = i)
      plot_cluster <- cbind(plot_cluster, ph_cluster)
    }
    plot_cluster$x <- NULL
    
    LMEscore.df <- cbind(LMEscore.df, plot_cluster)
    LMEscore.df$LME <- ifelse(LMEscore.df$ph_cluster == 1, "GC-like",
                              ifelse(LMEscore.df$ph_cluster == 2, "Inflammatory", 
                                     ifelse(LMEscore.df$ph_cluster == 3, "Depleted", "Mesenchymal")))
    LMEscore.df$LME <- factor(LMEscore.df$LME, levels = c("GC-like", "Mesenchymal", "Inflammatory", "Depleted"))
    print(table(LMEscore.df$LME))
    print(round(prop.table(table(LMEscore.df$LME)),2))
    
    print(table(LMEscore.df$LME, LMEscore.df$DLBCL90_Group_updated))
    
    LMEscore.df <- LMEscore.df %>% arrange(LME)
    anno <- LMEscore.df[,c("LME",exst)]
    mat <- LMEscore.df[,lme.targets]
    
    ph <- pheatmap::pheatmap(mat = t(mat), 
                             annotation_col = anno, 
                             show_colnames = F,
                             cluster_cols = F, cluster_rows = F,
                             scale = "row"
    )
    print(ph)
    
    cat("\n---- Draw LME heatmap (COO=ALL) ----\n")
    LMEscore.df$BCC_CD79b_Pattern <- factor(LMEscore.df$BCC_CD79b_Pattern, levels = c("M","C","NEG","UD"))
    LMEscore.df <- LMEscore.df %>% arrange(LME) %>% arrange(-BCC_CD79b_Hscore) 
    top_anno <- LMEscore.df[,c("DLBCL90_Group_updated","LME")]
    top_anno$LME <- factor(top_anno$LME, levels = c("GC-like","Mesenchymal","Inflammatory","Depleted"))
    top_anno$DLBCL90_Group_updated <- factor(top_anno$DLBCL90_Group_updated, levels = c("DZsig-pos","GCB","ABC","UNC"))
    top_anno <- 
      HeatmapAnnotation(df=top_anno, 
                        col = list(LME = c("GC-like"="#82b7ff","Mesenchymal"="#e5b200","Inflammatory"="#d39fff","Depleted"="#ffa242"),
                                   DLBCL90_Group_updated = c("DZsig-pos"="#078282","GCB"="#339e66","ABC"="#95DBe5","UNC"="gray90")
                                   ),
                                   simple_anno_size = unit(0.4, "cm")
                        )
    
    hscore <- LMEscore.df[,c("BCC_CD79b_Hscore", "BCC_CD79b_Hscore_class")]
    res <- summary(hscore$BCC_CD79b_Hscore)
    bottom_anno <- HeatmapAnnotation(CD79B_Hscore=anno_barplot(x=hscore$BCC_CD79b_Hscore,
                                                               gp=gpar(fill=hscore$BCC_CD79b_Hscore_class),
                                                               baseline=res[[3]]
                                                               ),
                                     height = unit(1, "cm")
                                     )
    
    mat <- LMEscore.df[,lme.targets]
    colnames(mat) <- c("LEC", "VEC", "CAF", "FRC", "ECM", "ECM_remodeling", "Granulocyte_traffic", 
                       "IS_cytokines", "FDC", "Macrophages", "Activated M1", "Tcells_traffic", 
                       "MHCI", "TFH", "Treg", "TIL", "IS_checkpoints", "NK", "Bcells_traffic", 
                       "Bcells", "Cell_proliferation")
    
    col_fun = colorRamp2(c(-1,0,1), c("#00bfff", "#ffffeb","#ff4000"))
    col_fun(seq(30))
    
    ha <- Heatmap(t(mat),
                  top_annotation = top_anno,
                  bottom_annotation = bottom_anno,
                  cluster_columns = F,
                  show_column_dend = F,
                  cluster_rows = F,
                  show_row_dend = F,
                  show_column_names = F,
                  col = col_fun,
                  use_raster = F, 
                  heatmap_height = unit(13.5, "cm"),
                  heatmap_width = unit(16.5, "cm")
                  )
    draw(ha)
    
    save_pdf(.heatmap = ha, "./saved/figs/BCC/LME_Heatmap(COO=ALL_Ordered=LME).pdf", width = 10, height = 8)
    
    print(table(LMEscore.df$LME, LMEscore.df$BCC_CD79b_Hscore_class) %>% chisq.test())
    print(table(LMEscore.df$LME, LMEscore.df$BCC_CD79b_Hscore_class) %>% chisq.test() %>% .[["stdres"]])
    
    save_pdf(.heatmap = ha, "./saved/figs/BCC/LME_Heatmap(COO=ABC_Ordered=LME).pdf", width = 10, height = 8)
    
    print(table(LMEscore.df[LMEscore.df$DLBCL90_Group_updated=="ABC",]$LME, 
                LMEscore.df[LMEscore.df$DLBCL90_Group_updated=="ABC",]$BCC_CD79b_Hscore_class) %>% chisq.test())
    print(table(LMEscore.df[LMEscore.df$DLBCL90_Group_updated=="ABC",]$LME, 
                LMEscore.df[LMEscore.df$DLBCL90_Group_updated=="ABC",]$BCC_CD79b_Hscore_class) %>% chisq.test() %>% .[["stdres"]])
    
    write.csv(LMEscore.df[,c("Patient_ID","LME",lme.targets)], "./saved/bulkRNAseq/BCC_LME.csv")
    }

scGCcluster <-
  function(){
    cat("\n---- scGC-cluster from Holmes AB et al. JEM 2020 ----\n")
    cat("\n---- load H-score data with clinical information (with UNCLASS) ----\n")
    hscore <- data.frame(data.table::fread("./saved/Scores/BCC_defined_cohort_withUNC.csv"))
    hscore <- hscore[hscore$available_clinicaldata=="Available"&hscore$available_IHCdata=="Available", ]
    ids <- hscore$Patient_ID
    
    cat("\n---- load GEP data (only protein coding genes) ----\n")
    cnt.nrm <- data.frame(data.table::fread("./saved/bulkRNAseq/BCC_CodingGENE_NrmCount(edgeR).csv"), row.names = 1)
    cnt.nrm <- cnt.nrm[rownames(cnt.nrm) %in% ids, ]
    
    holmes <- read.xlsx("./src/BCC/JEM2020/jem_20200483_tables5.xlsx", rowNames = F, startRow = 3, sheet = "Table_S5")
    exst <- c("Sample.ID","sc-COO.Class","sc-COO.Score","sc-COO.Group")
    holmes <- holmes[holmes$Dataset=="BCCA-DLBCL", exst]
    colnames(holmes) <- c("Patient_ID","JEMsupp_scGC_Class","JEMsupp_scGC_Score","JEMsupp_scGC_Group")

    hscore <- left_join(hscore, holmes, by="Patient_ID")
    
    hscore$JEMsupp_scGC_Class_Simple <- 
      ifelse(hscore$JEMsupp_scGC_Class %in% c("DZ a","DZ b","DZ c"), "DZ-like",
             ifelse(hscore$JEMsupp_scGC_Class %in% c("LZ a","LZ b"), "LZ-like",
                    ifelse(hscore$JEMsupp_scGC_Class %in% c("PBL a","PBL b"), "PBL-like",
                           ifelse(hscore$JEMsupp_scGC_Class == "PreM", "PreM-like", "INT-like"))))
    hscore$JEMsupp_scGC_Class_Simple <- factor(hscore$JEMsupp_scGC_Class_Simple,
                                               levels = c("DZ-like","INT-like","LZ-like","PBL-like","PreM-like"))
    
    
    cat("\n---- In BCC cohort, scGC cluster data is available in n =", 
        nrow(hscore[!is.na(hscore$JEMsupp_scGC_Class_Simple), ]), "samples (with COO=UNCLASS) ----\n")
    print(table(hscore$JEMsupp_scGC_Class_Simple))
    print(round(prop.table(table(hscore$JEMsupp_scGC_Class_Simple))*100,2))
    
    write.csv(hscore, "./saved/bulkRNAseq/BCC_scGCcluster.csv")
  }
































