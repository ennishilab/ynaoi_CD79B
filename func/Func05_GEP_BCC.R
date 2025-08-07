# --------------------- #
library(data.table)
library(dplyr)
library(edgeR)
library(DESeq2)
library(tibble)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GSVA)
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
    
    cnt.raw <- read.csv("./saved/bulkRNAseq/BCC_CodingGENE_RowCount.csv", row.names = 1)
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
    cnt.raw$TotalExpr_log <- log(cnt.raw$TotalExpr)
    cnt.raw$TotalExpr_log_scaled <- scale(cnt.raw$TotalExpr_log)
    print(table(cnt.raw$rand))
    
    pre <- ggboxplot(cnt.raw, 
                     x = "rand", 
                     y = "TotalExpr_log_scaled", 
                     outlier.shape = F,
                     ylim = c(-2,2)
    ) + 
      ylab("log total expression")+
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
    cnt.nrm$TotalExpr_log <- log(cnt.nrm$TotalExpr)
    cnt.nrm$TotalExpr_log_scaled <- scale(cnt.nrm$TotalExpr_log)
    
    post <- ggboxplot(cnt.nrm, 
                      x = "rand", 
                      y = "TotalExpr_log_scaled", 
                      outlier.shape = F,
                      ylim = c(-2,2)
    ) + 
      ylab("log total expression")+
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
    cnt.nrm$TotalExpr_log <- log(cnt.nrm$TotalExpr)
    cnt.nrm$TotalExpr_log_scaled <- scale(cnt.nrm$TotalExpr_log)
    
    post <- ggboxplot(cnt.nrm, 
                       x = "rand", 
                       y = "TotalExpr_log_scaled", 
                       outlier.shape = F,
                       ylim = c(-2,2)
    ) + 
      ylab("log total expression")+
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
    hscore <- read.csv("./saved/Scores/BCC_defined_cohort_withUNC.csv")
    hscore <- hscore[hscore$available_clinicaldata=="Available"&hscore$available_IHCdata=="Available", ]
    ids <- hscore$Patient_ID
    
    cat("\n---- load GEP data (only protein coding genes) ----\n")
    cnt.nrm <- read.csv("./saved/bulkRNAseq/BCC_CodingGENE_NrmCount(edgeR).csv", row.names = 1)
    cnt.nrm <- cnt.nrm[rownames(cnt.nrm) %in% ids, ]
    
    cat("\n---- In this study cohort, GEP data is available in n =", nrow(cnt.nrm), "samples ----\n")
    cat("\n---- In n =", nrow(cnt.nrm), "samples, we're going to calculate LME scores ----\n")
    
    cat("\n---- Basically, LME score is GSVA score ----\n")
      
    
    lme <- read.csv("./src/BCC/CancerDiscovery2021/CANCER_DISCOVERY_2021_LME_signatures.csv")
    LMEsig <- list()
    LMEsig <- apply(lme, MARGIN = 2, FUN = function(x){LMEsig <- c(LMEsig, x)})
    LMEsig <- lapply(LMEsig, FUN = function(x){x <- x[ x !=  ""] %>% unlist() })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
