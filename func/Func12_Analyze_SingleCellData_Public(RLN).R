# --------------------- #
library(Seurat)
library(ggplot2)
library(harmony)
library(sctransform)
library(dplyr)
library(escape)
library(UCell)
library(NMF)
library(pheatmap)
library(dittoSeq)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(tidyverse)
# --------------------- #

Analyze_Public_RLNs <-
  function(){
    .myfunc.env = new.env()
    sys.source( "./func/HolmesJEM_GCmarkers.R", envir = .myfunc.env )
    attach( .myfunc.env )
    Public_GCmarkers()
    
    markers.GC.UP <- readRDS("./src/scRNAseq/JEM2020/article/markers_UP.rds")
    
    markers.base <- c("PTPRC", #lymphocytes
                      "XBP1","IRF4", #Plasma-cell
                      "CD19","CD79A","MS4A1", #B-cell
                      "CD3D","CD3E", "CD4", "CD8A",#T-cell
                      "CD68","CD163","CD14", #MoMac
                      "COL1A1", #FB
                      "ACTA2", #SMC
                      "PECAM1", #EC
                      "SPARCL1","VWF", #BEC
                      "LYVE1", "PROX1", #LEC
                      "VCAM1", #FDC/FRC
                      #"ICAM1",
                      #"PDGFRA","PDGFRB",
                      "ITGA1", "CXCL12","CXCL14", #FRC
                      "CXCL13" #FDC
    )
    
    markers.base.list <- list(
      Non_GC_Bcells = c("CD19","IGHD"), 
      GC_Bcells = c("CD19","CD38"), 
      Plasmablast_plasmacytes = c("PRDM1","IRF4", "XBP1"), 
      MemoryBcells_PrecursorMemoryBcells = c("CCR6"),
      CD4_Tcells = c("CD3D","CD4"),
      CD8_Tcells = c("CD3D","CD8A"),
      Macrophages = c("CD68","CD163"),
      NKcells = c("NCAM1"), 
      FB <- c("COL1A1"),   
      EC <- c("PECAM1")
    )
    
    cat("\n=== Load SCTransform Normalized data ===\n")
    sc.obj <- readRDS("./saved/scRNAseq/merged/Healthy_Public_postnrm.rds")
    sc.obj <- FindNeighbors(sc.obj, dims = 1:50, assay = "SCT")
    sc.obj <- FindClusters(sc.obj, resolution = 0.02)
    
    dot <- DotPlot(sc.obj,
                   features = markers.base,
                   group.by = "SCT_snn_res.0.02"
    ) &
      scale_color_gradient2(high = "firebrick1", mid = "snow2", low = "dodgerblue2", midpoint = 0) &
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    dot
    
    Idents(sc.obj) <- "SCT_snn_res.0.02"
    new.cluster.ids <- c("T", "B", "B", "B", "B", "B", "B", "doublet_B/T", "PC", "MoMac")
    names(new.cluster.ids) <- levels(sc.obj)
    sc.obj <- RenameIdents(sc.obj, new.cluster.ids)
    sc.obj@meta.data$CellTypes_ver0 <- sc.obj@active.ident
    
    Idents(sc.obj) <- "SCT_snn_res.0.02"
    #dput(levels(sc.obj)) 
    new.cluster.ids <- c("T", "B1", "B2", "B3", "B4", "B5", "B6", "doublet_B/T", "PC", "MoMac")
    names(new.cluster.ids) <- levels(sc.obj)
    sc.obj <- RenameIdents(sc.obj, new.cluster.ids)
    sc.obj@meta.data$CellTypes_ver1 <- sc.obj@active.ident
    
    DimPlot(sc.obj, group.by = "CellTypes_ver0", pt.size = 1, raster = T) +
      DimPlot(sc.obj, group.by = "orig.ident", pt.size = 1, raster = T)
    
    # Extract B-cells & Plasmacytes ----
    sub.obj <- subset(sc.obj, subset = CellTypes_ver0 %in% c("B","PC"))
    cat("\n---- B&PC from public scRNA dataset has n =", nrow(sub.obj), "genes, 
        and n =", ncol(sub.obj), "cells ----\n")
    sub.obj <- SCTransform(sub.obj, assay = "RNA")
    sub.obj <- RunPCA(sub.obj, assay = "SCT")
    sub.obj <- RunUMAP(sub.obj, dims = 1:50)
    cat("\n---- After SCTransform, there are n =", nrow(sub.obj), "genes, 
        and n =", ncol(sub.obj), "cells ----\n")
    # sub.obj <- FindNeighbors(sub.obj, dims = 1:50, assay = "SCT")
    # sub.obj <- FindClusters(sub.obj, resolution = c(0.1,0.2,0.3))
    sub.obj <- JoinLayers(sub.obj, assay = "RNA")
    
    # Calculate UCell Score ----
    sub.obj <- AddModuleScore_UCell(sub.obj, features = markers.GC.UP)
    markers.GC.UP.UCell <- paste0(names(markers.GC.UP), "_UCell")
    f <- FeaturePlot(sub.obj, 
                     reduction = "umap", 
                     features = markers.GC.UP.UCell, 
                     ncol = 3, 
                     pt.size = 3,
                     order = T,
                     raster = T) & 
      scale_color_gradientn(colours = c("darkblue","skyblue", "lightgreen", "red", "darkred", "black"))
    print(f)
    
    saveRDS(sub.obj, "./saved/scRNAseq/processed/Public_RLNs_UCell.rds")
    
    
    # NMF clustering ----
    sub.obj <- readRDS("./saved/scRNAseq/processed/Public_RLNs_UCell.rds")
    df <- sub.obj@meta.data[, markers.GC.UP.UCell]
    colnames(df) <- gsub(pattern = "_UCell", replacement = "", x = colnames(df))
    res <- nmf(df, rank = 8, nrun=100, seed = 123456, .options="p")
    basismap(res, subsetRow=T)
    coefmap(res)
    consensusmap(res)
    
    res_w <- basis(res)
    res_h <- coef(res)
    
    ph <- pheatmap(data.frame(res_h))
    ph
    
    colnames(res_w) <- c("PreM/INTd","PBLa/PBLb","INTa",
                         "LZa/LZb","INTc",
                         "DZb","DZa/DZc/INTb","INTe")
    cluster <- apply(res_w, 1, FUN = function(x){
      colnames(res_w)[which(x == max(x))]
    })
    
    sub.obj$GC_cluster <- cluster
    sub.obj$GC_cluster <- factor(sub.obj$GC_cluster, 
                                 levels = c("DZb","INTa","DZa/DZc/INTb",
                                            "INTc","INTe","LZa/LZb",
                                            "PreM/INTd","PBLa/PBLb"
                                            ))
    Idents(sub.obj) <- "GC_cluster"
    
    cat("\n---- After SCTransform Normalization, Public_RLNs has",
        nrow(sub.obj), "genes, and", ncol(sub.obj), "cells ----\n")
    print(table(sub.obj$GC_cluster))
    print(round(prop.table(table(sub.obj$GC_cluster))*100,1))
    
    saveRDS(sub.obj, "./saved/scRNAseq/processed/Public_RLNs_UCell_NMF_GCclust.rds")
    }   









