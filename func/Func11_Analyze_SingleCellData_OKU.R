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

Analyze_citeRLN_OKU <-
  function(){
    
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
    sc.obj <- readRDS("./saved/scRNAseq/merged/Healthy_OKU_postnrm.rds")
    
    ref.genes <- rownames(sc.obj)
    write.table(ref.genes, "./saved/scRNAseq/processed/Healthy_OKU_genes.txt")
    
    sc.obj <- FindNeighbors(sc.obj, dims = 1:50, assay = "SCT")
    sc.obj <- FindClusters(sc.obj, resolution = c(0.02,0.04,0.06,0.1))
    
    dot <- DotPlot(sc.obj,
                   features = markers.base,
                   group.by = "SCT_snn_res.0.02"
    ) &
      scale_color_gradient2(high = "firebrick1", mid = "snow2", low = "dodgerblue2", midpoint = 0) &
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    dot
    
    Idents(sc.obj) <- "SCT_snn_res.0.02"
    levels(sc.obj) # [1] "0" "1" "2" "3" "4" "5" "6" "7"
    new.cluster.ids <- c("T", "B", "B", "PC/PBL", "MoMac", "T")
    names(new.cluster.ids) <- levels(sc.obj)
    sc.obj <- RenameIdents(sc.obj, new.cluster.ids)
    sc.obj@meta.data$CellTypes_ver0 <- sc.obj@active.ident
    
    Idents(sc.obj) <- "SCT_snn_res.0.02"
    levels(sc.obj) # [1] "0" "1" "2" "3" "4" "5" "6" "7"
    new.cluster.ids <- c("T1", "B1", "B2", "PC/PBL", "MoMac", "T2")
    names(new.cluster.ids) <- levels(sc.obj)
    sc.obj <- RenameIdents(sc.obj, new.cluster.ids)
    sc.obj@meta.data$CellTypes_ver1 <- sc.obj@active.ident
    
    DimPlot(sc.obj, group.by = "CellTypes_ver0", pt.size = 1, raster = T)
    
    # Extract GC B-cells & Plasmacytes & Memory Bcells ----
    sub.obj <- subset(sc.obj, subset = CellTypes_ver0 %in% c("B","PC/PBL"))
    sub.obj <- SCTransform(sub.obj, assay = "RNA")
    sub.obj <- RunPCA(sub.obj, assay = "SCT")
    sub.obj <- RunUMAP(sub.obj, dims = 1:50)
    # sub.obj <- FindNeighbors(sub.obj, dims = 1:50, assay = "SCT")
    # sub.obj <- FindClusters(sub.obj, resolution = c(0.1,0.2,0.3))
    sub.obj <- JoinLayers(sub.obj, assay = "ADT")
    
    # Calculate UCell Score ----
    gc()
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
    
    ggsave(filename = "./saved/figs/scRNAseq/scRLN_OKU_GC-UCell_Featureplots.pdf", 
           plot = f, width = 16, height = 20)
    
    saveRDS(sub.obj, "./saved/scRNAseq/processed/scRLN_OKU_UCell.rds")
    
    # NMF clustering ----
    gc()
    df <- sub.obj@meta.data[, markers.GC.UP.UCell]
    colnames(df) <- gsub(pattern = "_UCell", replacement = "", x = colnames(df))
    res <- nmf(df, rank = 7, nrun=100, seed = 123456, .options="p")
    basismap(res, subsetRow=T)
    coefmap(res)
    consensusmap(res)
    
    res_w <- basis(res)
    res_h <- coef(res)
    
    ph <- pheatmap(data.frame(res_h))
    ph
    
    colnames(res_w) <- c("PBLa/PBLb","INTc/INTe","DZa/DZc/INTb",
                         "LZa/LZb","INTa","DZb","PreM/INTd")
    cluster <- apply(res_w, 1, FUN = function(x){
      colnames(res_w)[which(x == max(x))]
    })
    
    sub.obj$GC_cluster <- cluster
    sub.obj$GC_cluster <- factor(sub.obj$GC_cluster, levels = sort(levels(factor(sub.obj$GC_cluster))))
    Idents(sub.obj) <- "GC_cluster"
    
    cat("\n---- After SCTransform Normalization, citeRLN_OKU has",
        nrow(sub.obj), "genes, and", ncol(sub.obj), "cells ----\n")
    print(table(sub.obj$GC_cluster))
    print(round(prop.table(table(sub.obj$GC_cluster))*100,1))
    
    saveRDS(sub.obj, "./saved/scRNAseq/processed/scRLN_OKU_UCell_NMF_GCclust.rds")
  }  


# CD79 expression ----
## part1 ----
sub.obj <- readRDS("./saved/scRNAseq/processed/scRLN_OKU_UCell_NMF_GCclust.rds")
dim(sub.obj) #[1] 16453  2447
DefaultAssay(sub.obj) <- "SCT"

# Supervised UMAP ----
cnt.sct.norm <- FetchData(object = sub.obj,
                          vars = rownames(sub.obj),
                          assay = "SCT",
                          layer = "data"
)
write.csv(cnt.sct.norm, "./saved/scRNAseq/processed/scRLN_OKU_sct-norm-counts.csv")

meta <- sub.obj@meta.data
write.csv(meta, "./saved/scRNAseq/processed/scRLN_OKU_metadata.csv")

## Supervised UMAP on Python
supervised_umap <- read.csv("./saved/scRNAseq/processed/scRLN_OKU_sUMAP_coords.csv")
#dim(supervised_umap)
colnames(supervised_umap)[1] <- "id"

umap_res <- sub.obj@reductions[["umap"]]@cell.embeddings
umap_res <- data.frame(umap_res)
umap_res$id <- rownames(umap_res)

nrow(umap_res) == nrow(supervised_umap) # [1] TRUE

sumap <- left_join(umap_res, supervised_umap, by = "id")
rownames(sumap) <- sumap$id
sumap <- sumap[,c("sUMAP_1", "sUMAP_2")]

sub.obj@reductions[["supervised.umap"]] <- sub.obj@reductions[["umap"]]
sub.obj@reductions[["supervised.umap"]]@cell.embeddings <- as.matrix(sumap)
sub.obj@reductions[["supervised.umap"]]@key <- "sUMAP_"

orders <- c("INTa","DZb","DZa/DZc/INTb","PreM/INTd","INTc/INTe","LZa/LZb","PBLa/PBLb")
sub.obj$GC_cluster <- factor(sub.obj$GC_cluster, levels = orders)

color_RefscGC <-
  c("#2ca02c", # INTa
    "#ff7f0e", # DZb
    "#d62728", # DZa/DZc/INTb,
    "#1f77b4", # PreM/INTd
    "#e377c2", # INTc/INTe
    "#9467bd", # LZa/LZb
    "#17becf"  # PBLa/PBLb
  )

d <- DimPlot(sub.obj, label = TRUE, label.size = 7,
             reduction = "supervised.umap",
             group.by = "GC_cluster", pt.size = 2, raster = TRUE,
             cols = color_RefscGC
)
f.rna <- FeaturePlot(sub.obj, features = "CD79B",
                     reduction = "supervised.umap", pt.size = 2, raster = TRUE) &
  scale_color_gradientn(colours = c("darkblue","skyblue", "red", "darkred"))
f.adt <- FeaturePlot(sub.obj, features = "Hu.CD79b",
                     reduction = "supervised.umap", pt.size = 2, raster = TRUE) &
  scale_color_gradientn(colours = c("darkblue","skyblue", "red", "darkred"))

ggsave("./saved/figs/scRNAseq/scRLN_OKU_GCclust_sUMAPs.pdf",
       plot = (d|f.rna|f.adt), width = 25, height = 7)

f <- FeaturePlot(sub.obj,
                 reduction = "supervised.umap",
                 features = markers.GC.UP.UCell,
                 ncol = 3,
                 pt.size = 3,
                 order = T,
                 raster = T) &
  scale_color_gradientn(colours = c("darkblue","skyblue", "lightgreen", "red", "darkred", "black"))
print(f)

ggsave(filename = "./saved/figs/scRNAseq/scRLN_OKU_GC-UCell_Featureplots_sUMAP.pdf",
       plot = f, width = 16, height = 20)

saveRDS(sub.obj, "./saved/scRNAseq/processed/scRLN_OKU_UCell_NMF_sUMAP.rds")

# fetch count data
bind <- data.frame(matrix(ncol = 2, nrow = 0))# 空のデータフレームを準備
colnames(bind) <- c("CD79B","Hu.CD79b")
foo <- FetchData(sub.obj, vars = "CD79B", assay="SCT", slot = "data")
bind[1:nrow(foo),"CD79B"] <- foo[,"CD79B"]
rownames(bind) <- rownames(foo)
foo <- FetchData(sub.obj, vars = "Hu.CD79b", assay="ADT", slot = "data")
bind[1:nrow(foo),"Hu.CD79b"] <- foo[,"adt_Hu.CD79b"]
rownames(bind) <- rownames(foo)

meta.sub <- sub.obj@meta.data[,c("orig.ident","GC_cluster")]
meta.sub$orig.ident <- NULL
meta.sub$cell_ids <- rownames(meta.sub)
bind$cell_ids <- rownames(bind)
bind <- left_join(bind, meta.sub, by="cell_ids")
rownames(bind) <- bind$cell_ids
bind$cell_ids <- NULL

medians_wide <- bind %>%
  dplyr::group_by(GC_cluster) %>%
  summarise(`CD79B(mRNA)`=median(CD79B), `CD79B(Protein)`=median(Hu.CD79b))

medians_long <- medians_wide %>%
  pivot_longer(cols = c("CD79B(mRNA)",
                        "CD79B(Protein)"),
               names_to = "CountType",
               values_to = "medianNrmCounts"
  )

medians_long <- medians_long %>%
  group_by(CountType) %>%
  mutate(medianNrmCounts_scaled = scale(medianNrmCounts))

medians_long$CountType <- factor(medians_long$CountType, levels = c("CD79B(Protein)","CD79B(mRNA)"))

p <- ggplot(medians_long,
            aes(x = reorder(GC_cluster, -medianNrmCounts),
                y = CountType,
                size = medianNrmCounts
            )) +
  geom_point(aes(colour = medianNrmCounts_scaled)) +
  scale_color_gradient2(high = "red", mid = "lightblue", low = "blue", midpoint = -0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./saved/figs/scRNAseq/scRLN_OKU_GCclust_CD79B_DotPlot.pdf",
       plot = p, width = 6, height = 4)


# CD79B protein and mRNA ----
cor_sp <- cor.test(bind$CD79B, bind$Hu.CD79b, method = "spearman", exact = FALSE)
print(cor_sp)

cors <- bind %>%
  summarize(cor = round(cor.test(CD79B, Hu.CD79b, method = "spearman",exact = F)$estimate,3), p = round(cor.test(CD79B, Hu.CD79b, method = "spearman",exact = F)$p.value,3)) %>%
  as.data.frame() %>% mutate(p = ifelse(p < 0.05, "p < 0.05", paste("p =", p )))

g <- ggplot(bind,aes(x=CD79B, y=Hu.CD79b))+
  geom_point(size=0.25)+
  geom_smooth(method="lm",formula='y~x',se=F, colour = "red")+
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 1, size = 16,
                                   hjust = 1, face = "bold"),
        axis.text.y = element_text(hjust = 1, size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),
        plot.title = element_text(size=22),
        plot.subtitle = element_text(size=20),
        legend.position = "none") +
  geom_text(data=cors,mapping = aes(x = -Inf, y = -Inf, label = paste("r =",cor)),hjust = -3,vjust = -3) +
  geom_text(data=cors,mapping = aes(x = -Inf, y = -Inf, label = p), hjust = -3,vjust = -1)
g

ggsave("./saved/figs/scRNAseq/scRLN_OKU_CD79B_CorrPlot.pdf",
       plot = g, width = 4, height = 4)

# Boxplot
stat_rna <- reshape2::melt(bind[,c("GC_cluster","CD79B")])
stat_rna$value_scaled <- scale(stat_rna$value)
res.kruskal.rna <- stat_rna %>% kruskal_test(value ~ GC_cluster)
stat.test.rna <-  stat_rna  %>%
  dunn_test(value ~ GC_cluster) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.rna <- stat.test.rna %>% mutate(y.position = 2)
cat("---- RNA ----")
print(stat.test.rna)

bxp.rna <- ggboxplot(
  data = stat_rna,
  x = "GC_cluster",
  y = "value_scaled",
  fill = "GC_cluster",
  outlier.shape = NA, ylim = c(-2,2.5)
)

g.rna <- bxp.rna +
  # stat_pvalue_manual(stat.test.rna, hide.ns = T)+
  labs(subtitle = get_test_label(res.kruskal.rna, detailed = TRUE),
       caption = get_pwc_label(stat.test.rna))+
  ylab("CD79B(mRNA)")+
  xlab("scGC-cluster")+
  theme(
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 4.2, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 5, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5, vjust = -3),
    plot.subtitle = element_text(size = 6, face = "bold", hjust = 1),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent', colour = NA),
    legend.box.background = element_rect(fill='transparent', colour = NA),
    text = element_text(family = "sans")
  ) +
  scale_fill_manual(values = color_RefscGC)

stat_adt <- reshape2::melt(bind[,c("GC_cluster","Hu.CD79b")])
stat_adt$value_scaled <- scale(stat_adt$value)
res.kruskal.adt <- stat_adt %>% kruskal_test(value ~ GC_cluster)
stat.test.adt <-  stat_adt  %>%
  dunn_test(value ~ GC_cluster) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.adt <- stat.test.adt %>% mutate(y.position = 2)
cat("---- adt ----")
print(stat.test.adt)

bxp.adt <- ggboxplot(
  data = stat_adt,
  x = "GC_cluster",
  y = "value_scaled",
  fill = "GC_cluster",
  outlier.shape = NA, ylim = c(-2,2.5)
)

g.adt <- bxp.adt +
  # stat_pvalue_manual(stat.test.adt, hide.ns = T)+
  labs(subtitle = get_test_label(res.kruskal.adt, detailed = TRUE),
       caption = get_pwc_label(stat.test.adt))+
  ylab("CD79B(Protein)")+
  xlab("scGC-cluster")+
  theme(
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 4.2, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 5, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5, vjust = -3),
    plot.subtitle = element_text(size = 6, face = "bold", hjust = 1),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent', colour = NA),
    legend.box.background = element_rect(fill='transparent', colour = NA),
    text = element_text(family = "sans")
  ) +
  scale_fill_manual(values = color_RefscGC)

print(g.rna|g.adt)

ggsave("./saved/figs/scRNAseq/scRLN_OKU_CD79B_BoxPlots.pdf",
       plot = (g.rna|g.adt), width = 8, height = 4)








