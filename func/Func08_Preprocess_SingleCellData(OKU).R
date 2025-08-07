# --------------------- #
library(Seurat)
library(sctransform)
library(rlist)
library(pipeR)
# --------------------- #

MergeSingleCellData <- 
  function(path, mode){
    cat("\n=== Merge Seurat Objects ===\n")
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    
    files <- list.files(path)
    for(file in files){
      cat("\n", file, "")
    }
    
    if(mode=="Healthy_OKU"){
      targets <- files[grep("citeRLN_OKU", files, invert = F)]
    }
    
    if(mode=="Healthy_Public"){
      targets <- files[grep("citeRLN_OKU|scFL|scDLBCL|scCHL", files, invert = T)]
      print(targets)
    }
    
    merge.list <- list()
    for(target in targets){
      cat("\n==> loading", target, "==>")
      sc.obj = readRDS(paste0(path, target))
      cat("\n originally, ", nrow(sc.obj), "genes,", ncol(sc.obj), "cells  \n")
      merge.list[[target]] <- sc.obj
    }
    
    all.genes <- c()
    for(target in targets){
      all.genes <- append(all.genes, rownames(merge.list[[target]]))
    }
    all.genes <- data.frame(all.genes)
    colnames(all.genes)[1] <- "genes"
    dupli.genes <- all.genes %>%
      group_by(genes) %>%
      filter(n()==length(targets))
    cat("\n n =", nrow(dupli.genes), "genes are duplicated across dataset ----\n")
    
    merged <- merge.list[[1]]
    for(i in 2:length(merge.list)){
      merged <- merge(x=merged, y=merge.list[[i]])
    }
    
    gc()
    cat("\n---- Merged data:", nrow(merged), "genes,", ncol(merged), "cells ----\n")
    print(table(merged@meta.data$orig.ident))
    DefaultAssay(merged) <- "RNA"
    merged <- JoinLayers(merged, assay="RNA")
    merged@assays$RNA <- split(x = merged@assays$RNA, f = merged@meta.data$orig.ident)
    
    saveRDS(merged, paste0("./saved/scRNAseq/merged/",mode,"_prenrm.rds"))

    cat("\n---- Filtering ----\n")
    cat("\n---- subset = nFeature_RNA > 200 & nCount_RNA < 25000 & percent.mt < 5 ----\n")
    merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
    merged <- subset(merged, subset = nFeature_RNA > 200 & nCount_RNA < 25000 & percent.mt < 5)
    cat("\n---- Filtered data:", nrow(merged), "genes,", ncol(merged), "cells ----\n")

    cat("\n---- Normalizing (SCTransform) ----\n")
    merged <- SCTransform(merged, assay = "RNA", verbose = T)
    cat("\n---- Normalized data:", nrow(merged), "genes,", ncol(merged), "cells ----\n")

    cat("\n---- PCA ----\n")
    merged <- RunPCA(merged, assay = "SCT", verbose = T)

    cat("\n---- UMAP ----\n")
    merged <- RunUMAP(merged, reduction = "pca", dims = 1:50, verbose = T)

    cat("\n---- Processing Done ----\n")
    cat("\n---- Processed data:", nrow(merged), "genes,", ncol(merged), "cells ----\n")
    print(table(merged@meta.data$orig.ident))

    cat("\n---- Saving ----\n")
    saveRDS(merged, paste0("./saved/scRNAseq/merged/",mode,"_postnrm.rds"))
    cat("\n---- Saving Done ----\n")
  }
