# --------------------- #
library(Seurat)
library(sctransform)
library(ggplot2)
library(dplyr)
library(dsb)
library(data.table)
library(tidyr)
library(tidyverse)
library(magrittr)
# --------------------- #

# CITEseq
CreateSeuratObj_citeRLN_OKU <- 
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n ====== Start Processing OKU CITEseq data ====== \n")
    
    cat("\n=== Loading S010 Data ===\n")
    # read raw data using the Seurat function "Read10X" 
    raw = Seurat::Read10X("./src/scRNAseq/CITEseq/src/S010_raw/outs/multi/count/raw_feature_bc_matrix/")
    cells = Seurat::Read10X("./src/scRNAseq/CITEseq/src/S010/count/sample_filtered_feature_bc_matrix/")
    # define a vector of cell-containing barcodes and remove them from unfiltered data 
    stained_cells = colnames(cells$`Gene Expression`)
    background = setdiff(colnames(raw$`Gene Expression`), stained_cells)
    # split the data into separate matrices per assay 
    prot = raw$`Antibody Capture`
    rna = raw$`Gene Expression`
    # create metadata of droplet QC stats used in standard scRNAseq processing
    rna_size = log10(Matrix::colSums(rna))
    prot_size = log10(Matrix::colSums(prot))
    ngene = Matrix::colSums(rna > 0)
    mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE)
    propmt = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
    md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
    md$bc = rownames(md)
    md$droplet_class = ifelse(test = md$bc %in% stained_cells, yes = 'cell', no = 'background')
    # filter barcodes to only include those with data for both assays 
    md = md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )
    # Quality control on cell-containing and background droplets
    ggplot(md, aes(x = log10(ngene), y = prot_size )) +
      theme_bw() + 
      geom_bin2d(bins = 300) + 
      scale_fill_viridis_c(option = "C") + 
      facet_wrap(~droplet_class) 
    
    cellmd = md %>% dplyr::filter(droplet_class == 'cell')
    plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
    p1 = ggplot(cellmd, aes(x = rna_size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
    p2 = ggplot(cellmd, aes(x = propmt)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
    p3 = ggplot(cellmd, aes(x = log10(ngene), y = rna_size, fill = propmt )) + plot_aes
    p4 = ggplot(cellmd, aes(x = ngene, y = prot_size, fill = propmt )) + plot_aes
    p1+p2+p3+p4
    
    # calculate statistical thresholds for droplet filtering. 
    rna_size_min = median(cellmd$rna_size) - (3*mad(cellmd$rna_size))
    rna_size_max = median(cellmd$rna_size) + (3*mad(cellmd$rna_size))
    prot_size_min = median(cellmd$prot_size) - (3*mad(cellmd$prot_size))
    prot_size_max = median(cellmd$prot_size) + (3*mad(cellmd$prot_size))
    
    # filter rows based on droplet qualty control metrics
    positive_cells = cellmd[
      cellmd$prot_size > prot_size_min & 
        cellmd$prot_size < prot_size_max & 
        cellmd$propmt < 0.14 &  
        cellmd$rna_size > rna_size_min & 
        cellmd$rna_size < rna_size_max, ]$bc
    cells_mtx_rawprot = as.matrix(prot[ , positive_cells])
    
    # Sanity check: are the number of QCd cells in line with the expected recovery from the experiment?
    length(positive_cells) # 2920
    
    # define a vector of background droplet barcodes based on protein library size and mRNA content
    background_drops = md[md$prot_size > 1.5 & md$prot_size < 3 & md$ngene < 100, ]$bc
    negative_mtx_rawprot = as.matrix(prot[ , background_drops])
    
    #normalize protein data for the cell containing droplets with the dsb method. 
    dsb_norm_prot = dsb::DSBNormalizeProtein(
      cell_protein_matrix = cells_mtx_rawprot, 
      empty_drop_matrix = negative_mtx_rawprot, 
      denoise.counts = TRUE, 
      use.isotype.control = TRUE, 
      isotype.control.name.vec = rownames(cells_mtx_rawprot)[29:31] 
    )
    dsb_norm_prot = apply(dsb_norm_prot, 2, function(x){ ifelse(test = x < -10, yes = 0, no = x)})
    
    # filter raw protein, RNA and metadata to only include cell-containing droplets 
    cells_rna = rna[ ,positive_cells]
    md2 = md[positive_cells, ]
    
    # create Seurat object !note: min.cells is a gene filter, not a cell filter
    s = Seurat::CreateSeuratObject(counts = cells_rna, 
                                   meta.data = md2, 
                                   assay = "RNA", 
                                   project = "citeRLN_OKU_S010")
    
    # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
    s@assays[["ADT"]] = CreateAssay5Object(data = dsb_norm_prot, key = "adt_")
    s -> S010_dsb
    cat("\n---- OKU_citeRLN_S010 has", 
        nrow(S010_dsb), "genes,", 
        ncol(S010_dsb), "cells ----\n")
    
    S010_dsb@meta.data[["orig.ident"]] <- rep("citeRLN_OKU_S010")
    saveRDS(S010_dsb, "./saved/scRNAseq/src/citeRLN_OKU_S010.rds")
    
    cat("\n=== Loading S073 Data ===\n")
    # read raw data using the Seurat function "Read10X" 
    raw = Seurat::Read10X("./src/scRNAseq/CITEseq/src/S073_raw/outs/multi/count/raw_feature_bc_matrix/")
    cells = Seurat::Read10X("./src/scRNAseq/CITEseq/src/S073/count/sample_filtered_feature_bc_matrix/")
    # define a vector of cell-containing barcodes and remove them from unfiltered data 
    stained_cells = colnames(cells$`Gene Expression`)
    background = setdiff(colnames(raw$`Gene Expression`), stained_cells)
    # split the data into separate matrices per assay 
    prot = raw$`Antibody Capture`
    rna = raw$`Gene Expression`
    # create metadata of droplet QC stats used in standard scRNAseq processing
    rna_size = log10(Matrix::colSums(rna))
    prot_size = log10(Matrix::colSums(prot))
    ngene = Matrix::colSums(rna > 0)
    mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE)
    propmt = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
    md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
    md$bc = rownames(md)
    md$droplet_class = ifelse(test = md$bc %in% stained_cells, yes = 'cell', no = 'background')
    # filter barcodes to only include those with data for both assays 
    md = md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )
    # Quality control on cell-containing and background droplets
    ggplot(md, aes(x = log10(ngene), y = prot_size )) +
      theme_bw() + 
      geom_bin2d(bins = 300) + 
      scale_fill_viridis_c(option = "C") + 
      facet_wrap(~droplet_class) 
    
    cellmd = md %>% dplyr::filter(droplet_class == 'cell')
    plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
    p1 = ggplot(cellmd, aes(x = rna_size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
    p2 = ggplot(cellmd, aes(x = propmt)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
    p3 = ggplot(cellmd, aes(x = log10(ngene), y = rna_size, fill = propmt )) + plot_aes
    p4 = ggplot(cellmd, aes(x = ngene, y = prot_size, fill = propmt )) + plot_aes
    p1+p2+p3+p4
    
    # calculate statistical thresholds for droplet filtering. 
    rna_size_min = median(cellmd$rna_size) - (3*mad(cellmd$rna_size))
    rna_size_max = median(cellmd$rna_size) + (3*mad(cellmd$rna_size))
    prot_size_min = median(cellmd$prot_size) - (3*mad(cellmd$prot_size))
    prot_size_max = median(cellmd$prot_size) + (3*mad(cellmd$prot_size))
    
    # filter rows based on droplet qualty control metrics
    positive_cells = cellmd[
      cellmd$prot_size > prot_size_min & 
        cellmd$prot_size < prot_size_max & 
        cellmd$propmt < 0.14 &  
        cellmd$rna_size > rna_size_min & 
        cellmd$rna_size < rna_size_max, ]$bc
    cells_mtx_rawprot = as.matrix(prot[ , positive_cells])
    
    # Sanity check: are the number of QCd cells in line with the expected recovery from the experiment?
    length(positive_cells) # 2474
    
    # define a vector of background droplet barcodes based on protein library size and mRNA content
    background_drops = md[md$prot_size > 1.5 & md$prot_size < 3 & md$ngene < 100, ]$bc
    negative_mtx_rawprot = as.matrix(prot[ , background_drops])
    
    #normalize protein data for the cell containing droplets with the dsb method. 
    dsb_norm_prot = dsb::DSBNormalizeProtein(
      cell_protein_matrix = cells_mtx_rawprot, 
      empty_drop_matrix = negative_mtx_rawprot, 
      denoise.counts = TRUE, 
      use.isotype.control = TRUE, 
      isotype.control.name.vec = rownames(cells_mtx_rawprot)[29:31] 
    )
    dsb_norm_prot = apply(dsb_norm_prot, 2, function(x){ ifelse(test = x < -10, yes = 0, no = x)})
    
    # filter raw protein, RNA and metadata to only include cell-containing droplets 
    cells_rna = rna[ ,positive_cells]
    md2 = md[positive_cells, ]
    
    # create Seurat object !note: min.cells is a gene filter, not a cell filter
    s = Seurat::CreateSeuratObject(counts = cells_rna, 
                                   meta.data = md2, 
                                   assay = "RNA", 
                                   project = "citeRLN_OKU_S073")
    
    # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
    s@assays[["ADT"]] = CreateAssay5Object(data = dsb_norm_prot, key = "adt_")
    s -> S073_dsb
    cat("\n---- OKU_citeRLN_S073 has", 
        nrow(S073_dsb), "genes,", 
        ncol(S073_dsb), "cells ----\n")
    
    S073_dsb@meta.data[["orig.ident"]] <- rep("citeRLN_OKU_S073")
    saveRDS(S073_dsb, "./saved/scRNAseq/src/citeRLN_OKU_S073.rds")
    
    cat("\n ====== Finished Processing OKU CITEseq data ====== \n")
    }

# scRNAseq
CreateSeuratObj_scRLN_JEM2020 <- 
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== citeLungLN: Start Processing Holmes AB et al. JEM 2020 ----\n")
    
    filepath <- list.files(path = "./src/scRNAseq/JEM2020/src/", 
                           pattern = "umi", 
                           full.names = T)
    GClist <- lapply(filepath, function(x){data.frame(fread(x, header = TRUE, sep = "\t"), row.names = 1)})
    names(GClist) <- paste0("scRLN_JEM2020_", stringr::str_split(string = filepath, pattern = "_", simplify = T)[,2])
    GClist <-
      lapply(GClist, function(x){
        x <- x[!duplicated(stringr::str_split(string = rownames(x), pattern = ";", simplify = T)[,2]),]
      })
    
    for(i in 1:length(GClist)){
      rownames(GClist[[i]]) <- stringr::str_split(string = rownames(GClist[[i]]), pattern = ";", simplify = T)[,2]
    }
    
    GClist <-
      mapply(
        function(x,y){CreateSeuratObject(counts = x, project = y)},
        x = GClist,
        y= names(GClist),
        SIMPLIFY = FALSE
      )
    
    GClist <-
      mapply(
        function(x,y){RenameCells(x, add.cell.id=y)},
        x = GClist,  y= names(GClist),
        SIMPLIFY = FALSE
      )
    
    seurat.object <- GClist[[1]]
    for(i in 2:length(GClist)){
      seurat.object <- merge(x=seurat.object, y=GClist[[i]])
      seurat.object@meta.data$orig.ident <- rep("scRLN_JEM2020")
    }
    cat("\n---- scRLN JEM2020 has", 
        nrow(seurat.object), "genes,", 
        ncol(seurat.object), "cells ----\n")
    
    print(table(seurat.object@meta.data$orig.ident))
    
    saveRDS(seurat.object, "./saved/scRNAseq/src/scRLN_JEM2020.rds")
  }

CreateSeuratObj_scRLN_SciImmunol2021 <-
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== scRLN: Start Processing King et al. Sci Immunol 2021 ----\n")
    seurat.object <- readRDS("./src/scRNAseq/SciImmunol2021/src/E-MTAB-9005/SeuratGeneExpressionObjects/SEURAT_OBJECTS/CompleteIntegrated_scRNA_SeuratObject.rds")
    seurat.object <- UpdateSeuratObject(object = seurat.object)
    DefaultAssay(seurat.object) <- "RNA"
    
    meta <- seurat.object@meta.data
    exst <- c("Donor", "CellType", "Subset",
              "S.Score", "G2M.Score", "DoubletFinderPrediction", 
              "ISOTYPE")
    meta <- meta[,exst]
    colnames(meta) <- c("Donor","CellType_Simple","CellType_Detail",
                        "S.Score", "G2M.Score", "DoubletFinderPrediction", 
                        "IG_ISOTYPE")
    
    cnt.raw <- FetchData(object = seurat.object, 
                         vars = rownames(seurat.object), 
                         assay = "RNA",
                         layer = "counts")
    gc() 
    sc.obj <- CreateSeuratObject(counts = t(cnt.raw), project = "scRLN_SciImmunol2021")
    sc.obj@meta.data$orig.ident <- rep("scRLN_SciImmunol2021")
    sc.obj <- AddMetaData(object = sc.obj, metadata = meta)
    sc.obj@assays$RNA <- split(x = sc.obj@assays$RNA, f = sc.obj@meta.data$Donor)
    
    cat("\n---- scRLN_SciImmunol2021 has", 
        nrow(sc.obj), "genes,", 
        ncol(sc.obj), "cells ----\n")
    print(table(sc.obj@meta.data$orig.ident))
    
    saveRDS(sc.obj, "./saved/scRNAseq/src/scRLN_SciImmunol2021.rds")
  }

CreateSeuratObj_scMixed_CancerCell2021 <-
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== scRLN: Start Processing Steen et al. Cancer Cell 2021 ----\n")
    cnt <- data.frame(fread(file = "./src/scRNAseq/CancerCell2021/src/GSE182434_raw_count_matrix.txt" , header = TRUE, sep = "\t"), row.names=1)
    meta <- data.frame(fread(file = "./src/scRNAseq/CancerCell2021/src/GSE182434_cell_annotation.txt" , header = TRUE, sep = "\t"), row.names=1)
    print(table(meta$Tissue))
    
    seurat.object <- CreateSeuratObject(counts = cnt)
    seurat.object <- AddMetaData(object = seurat.object, metadata = meta)
    
    tissues <- unique(meta$Tissue)
    for(tissue in tissues){
      sub <- subset(seurat.object, subset = Tissue == tissue)
      sub@meta.data$orig.ident <- rep(paste0("sc",tissue,"_CancerCell2021"))
      sub@project.name <- paste0("sc",tissue,"_CancerCell2021")
      cat("project =",sub@project.name)
      print(table(sub$orig.ident))
      cat("\n---- Steen et al. Cancer Cell 2021", tissue, "has", 
          nrow(sub), "genes,", 
          ncol(sub), "cells ----\n")
      saveRDS(sub, paste0("./saved/scRNAseq/src/sc",tissue,"_CancerCell2021.rds"))
    }
  }

CreateSeuratObj_scMixed_CancerDiscov2019 <-
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== scMixed: Start Processing Aoki et al. Cancer Discov 2019 ----\n")
    sc.obj <- readRDS("./src/scRNAseq/CancerDiscov2019/src/SingleCellExperiment_object")
    count <- counts(sc.obj)
    seurat.object <- CreateSeuratObject(counts = count)
    
    cols <- data.frame(colnames(seurat.object))
    colnames(cols) <- "ColName"
    rownames(cols) <-cols$ColName
    cols <- cols %>% separate(ColName, c("ColName_omit", "orig.ident"),sep = "_")
    cols$ColName_omit <- NULL
    for (i in 1:nrow(cols)){
      cols$Type[i] <- substr(cols$orig.ident[i], 1, 3)
    }
    seurat.object <- AddMetaData(seurat.object, metadata = cols)
    print(table(seurat.object@meta.data$Type))
    
    types <- unique(seurat.object@meta.data$Type)
    for(type in types){
      sub <- subset(seurat.object, subset = Type == type)
      sub@project.name <- paste0("sc",type,"_CancerDiscov2019")
      sub$orig.ident <- rep(paste0("sc",type,"_CancerDiscov2019"))
      cat("\n---- Aoki et al. Cancer Discov 2019", type, "has", 
          nrow(sub), "genes,", 
          ncol(sub), "cells ----\n")
      print(table(sub$orig.ident))
      saveRDS(sub, paste0("./saved/scRNAseq/src/sc",type,"_CancerDiscov2019.rds"))
    }
  }

CreateSeuratObj_scRLN_SciImmunol2022 <- 
  function(){
    gc()
    cat("\n=== scRLN: Start analyzing Siu JHY et al. SciImmunol 2022 ----\n")
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    
    targets <- c("A_MLN","B_MLN","C_MLN")
    merge.list <- list()
    for(target in targets){
      cat("\n---- loading", target, "----\n")
      counts.rna = Read10X(paste0("./src/scRNAseq/SciImmunol2022/src/", target))
      sc.obj <- CreateSeuratObject(counts = counts.rna, project = "scRLN_SciImmunol2022")
      cat("\n----", target, "has", nrow(sc.obj), "genes,", ncol(sc.obj), "cells ----\n")
      
      merge.list[[target]] <- sc.obj
      
      # vdj.b <- data.frame(data.table::fread(paste0("./src/scRNAseq/SciImmunol2022/src/vdj_b/",target,".csv"), sep=","))
      # vdj.b <- vdj.b[!duplicated(vdj.b$barcode),]
      # rownames(vdj.b) <- vdj.b$barcode
      # vdj.b <- vdj.b[,c("raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene")]
      # sc.obj <- AddMetaData(sc.obj, metadata = vdj.b)
    }
    
  seurat.object <- merge.list[[1]]
  for(i in 2:length(merge.list)){
    seurat.object <- merge(x=seurat.object, y=merge.list[[i]])
  }
  seurat.object@meta.data$orig.ident <- rep("scRLN_SciImmunol2022")
  cat("\n---- scRLN_SciImmunol2022 has", 
      nrow(seurat.object), "genes,", 
      ncol(seurat.object), "cells ----\n")
  saveRDS(seurat.object, "./saved/scRNAseq/src/scRLN_SciImmunol2022.rds")
  }

CreateSeuratObj_scMixed_NatCellBiol2020 <- 
  function(){
    gc()
    date <- Sys.Date()
    cat("\n Date:",format(date, "%d/%m/%Y"),"\n")
    cat("\n=== scMixed: Start analyzing Roider et al. NatCellBiol 2020 ----\n")
    
    targets <- c("rLN1","rLN2","rLN3")
    merge.list <- list()
    for(target in targets){
      cat("\n---- loading", target, "----\n")
      counts.rna = Read10X(paste0("./src/scRNAseq/NatCellBiol2020/src/", target))
      sc.obj <- CreateSeuratObject(counts = counts.rna, project = "scRLN_NatCellBiol2020")
      merge.list[[target]] <- sc.obj
    }
    
    seurat.object <- merge.list[[1]]
    for(i in 2:length(merge.list)){
      seurat.object <- merge(x=seurat.object, y=merge.list[[i]])
    }
    seurat.object@meta.data$orig.ident <- rep("scRLN_NatCellBiol2020")
    cat("\n---- scRLN_NatCellBiol2020 has", 
        nrow(seurat.object), "genes,", 
        ncol(seurat.object), "cells ----\n")
    saveRDS(seurat.object, "./saved/scRNAseq/src/scRLN_NatCellBiol2020.rds")
    
    
    targets <- c("DLBCL1","DLBCL2","DLBCL3")
    for(target in targets){
      cat("\n---- loading", target, "----\n")
      counts.rna = Read10X(paste0("./src/scRNAseq/NatCellBiol2020/src/", target))
      sc.obj <- CreateSeuratObject(counts = counts.rna, project = paste0("scDLBCL_NatCellBiol2020_",target))
      cat("\n----", target, "has", nrow(sc.obj), "genes,", ncol(sc.obj), "cells ----\n")
       print(table(sc.obj$orig.ident))
      saveRDS(sc.obj, paste0("./saved/scRNAseq/src/scDLBCL_NatCellBiol2020_",target,".rds"))
    }
    
    targets <- c("FL1","FL2","FL3","FL4","tFL1","tFL2")
    for(target in targets){
      cat("\n---- loading", target, "----\n")
      counts.rna = Read10X(paste0("./src/scRNAseq/NatCellBiol2020/src/", target))
      sc.obj <- CreateSeuratObject(counts = counts.rna, project = paste0("scFL_NatCellBiol2020_",target))
      cat("\n----", target, "has", nrow(sc.obj), "genes,", ncol(sc.obj), "cells ----\n")
      print(table(sc.obj$orig.ident))
      saveRDS(sc.obj, paste0("./saved/scRNAseq/src/scFL_NatCellBiol2020_",target,".rds"))
    }
  }




























