# Hello, world!
#' @param rds1 The path to the RDS file.
#' @param exported_file_name The name of the exported CSV file.
#' @param plot1_name The name of the AD-specific genes heatmap plot.
#' @param plot2_name The name of the PV-specific genes heatmap plot.
#' @param plot3_name The name of the sample-level means plot without individual cells.
#' @param plot4_name The name of the sample-level means plot with individual cells.
#' @param original_file_name The name of the original file.
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
.onAttach <- function(libname, pkgname) {
  bioconductor_packages <- c("MAST", "multtest", "scDblFinder", "ComplexHeatmap", "batchelor")
  missing_packages <- bioconductor_packages[!bioconductor_packages %in% rownames(installed.packages())]

  if (length(missing_packages) > 0) {
    message("To use all functionalities of ", pkgname, ", install the following Bioconductor packages: ",
            paste(missing_packages, collapse = ", "),
            "\nYou can install them using the following command:\n",
            "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); ",
            "BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))")
  }
}
#' @import Seurat
#' @import metap
#' @import MAST
#' @import scDblFinder
#' @import ggplot2
#' @import cowplot
#' @import patchwork
#' @import stringr
#' @import RColorBrewer
#' @import ComplexHeatmap
#' @import circlize
#' @import tidyverse
#' @import concaveman
#' @import ggforce
#' @import monocle3
#' @import dplyr
#' @import Rmisc
#' @import ggrepel
#' @import raster
#' @import vegan
#' @import plyr
#' @import SummarizedExperiment
#' @import grid
#' @export
process_requests <- function(data_file_path,rds1, exported_file_name,
                             plot1_name = "plot1.png",
                             plot2_name = "plot2.png",
                             plot3_name = "plot3.png",
                             plot4_name = "plot4.png",
                             original_file_name = "original_data.csv") {
  print("Hello, world!")
  # data_file_path <- system.file("extdata", "skin_integrated.rds", package = "RashxPackage")
  skin.integrated <- readRDS(data_file_path)
  print(colnames(skin.integrated))

  print(original_file_name)

  tryCatch({ #Start error catch block
    skinX.query <- readRDS(rds1)

  }, error=function(e){print("ERROR: Submission dataset read error")}) #End error catch block and print out error message

  tryCatch({ #Start error catch block

    #5. Cell type classification using an integrated reference
    #(i.e. Mapping and annotating external query datasets based on our reference datasets to define Trm1 cells).
    #find anchors and predict.id
    skinX.anchors <-
      FindTransferAnchors(
        reference = skin.integrated,
        query = skinX.query,
        dims = 1:30,
        reference.reduction = 'pca'
      )

  }, error=function(e){print("ERROR: Metadata transfer anchors not found.")}) #End error catch block and print out error message


  tryCatch({ #Start error catch block

    predictionsX <-
      TransferData(anchorset = skinX.anchors,
                   refdata = skin.integrated$ID,
                   dims = 1:30)
    skinX.query <- AddMetaData(skinX.query, metadata = predictionsX)

  }, error=function(e){print("ERROR: Cell attribute transfer unsuccessful after integration anchor identification.")}) #End error catch block and print out error message



  skinX.query <- NormalizeData(skinX.query)

  #6. DEGs heatmap of external query dataset for predicted.id = 2 (or Trm1) cells
  # subset cells with predicted.id = 2
  skinX.query2 <-
    subset(x = skinX.query, subset = predicted.id == "2")

  #assign condition
  skinX.query2$dis <- "RashX"

  skinX.query2$donor <- "RashX"

  # subset reference data cells with ID=2
  skin.integrated2 <- subset(x = skin.integrated, subset = ID == "2")

  tryCatch({ #Start error catch block


    #Get integration features
    features <- SelectIntegrationFeatures(object.list = list(skin.integrated2,skinX.query2))

    #Integration anchors
    #Speed up FindInetegrationAnchors: https://github.com/satijalab/seurat/discussions/3999
    #plan("multisession", workers = 4) # Enable parallelization from future package
    options(future.globals.maxSize = 20000 * 1024^2)
    cell.anchors <- FindIntegrationAnchors(object.list = list(skin.integrated2,skinX.query2),
                                           anchor.features = features,
                                           #reference = c(1,2),
                                           reduction = "rpca",
                                           dims = 1:30
    )

    # merge reference and query datasets
    Unknown <- IntegrateData(anchorset = cell.anchors,
                             k.weight = 10) # https://github.com/satijalab/seurat/issues/3930

  }, error=function(e){print("ERROR: Unable to merge submitted query data with reference data, likely no Trm cells identified.")}) #End error catch block and print out error message


  # run DEGs for AD- or Pso-specific genes
  genes <-
    c(
      "TWIST1",
      "LGALS1",
      "IL32",
      "CAPG",
      "ITM2C",
      "MFHAS1",
      "ANXA1",
      "SOS1",
      "CSGALNACT1",
      "LMO4",
      "IFITM2",
      "S100A10",
      "MT-ND5",
      "CYSLTR1",
      "PLA2G16",
      "SYNE2",
      "THADA",
      "NEAT1",
      "IL17RB",
      "RPL36A",
      "ARHGAP21",
      "NBAS",
      "ACTG1",
      "PRKX",
      "TGFBR3",
      "TIMP1",
      "TNFSF10",
      "AHNAK",
      "MT-ND2",
      "ISG15",
      "RPL17",
      "LONRF2",
      "CD99",
      "TSHZ2",
      "MMP25",
      "IFITM1",
      "MT-ND1",
      "BIRC3",
      "FAM102A",
      "LPCAT2",
      "NRIP3",
      "CRIP1",
      "CLU",
      "PLP2",
      "ZFP36",
      "ZFP36L2",
      "TUBA1B",
      "GATA3",
      "SLC5A3",
      "SFXN1",
      "FANK1",
      "TAGLN2",
      "CXCL13",
      "MTRNR2L12",
      "CD7",
      "MGAT4A",
      "FTH1",
      "LAYN",
      "IL17F",
      "KLRB1",
      "GNLY",
      "CPM",
      "CTSH",
      "GBP5",
      "SOX4",
      "CLEC2B",
      "GZMB",
      "CD2",
      "CEBPD",
      "ODF2L",
      "LAG3",
      "LRRN3",
      "ARHGEF12",
      "PTPN13",
      "TNFAIP3",
      "TRPS1",
      "SNX9",
      "METRNL",
      "BTG1",
      "JUN",
      "SPOCK2",
      "GABARAPL1",
      "PMEPA1",
      "HIST1H1E",
      "RBPJ",
      "LINC01871",
      "MAP3K4",
      "H1FX",
      "UBC",
      "GALNT1",
      "PNRC1",
      "GABPB1-AS1",
      "RPS26",
      "MUC20-OT1",
      "CHN1",
      "NAP1L4",
      "PTMS",
      "F2R",
      "CTLA4",
      "DAPK2",
      "RAP1B",
      "CCR6",
      "B3GALT2",
      "YPEL2",
      "FYN",
      "PPDPF",
      "SLA2",
      "CBLB",
      "ADGRG1",
      "SARAF"
    )

  #updated ad_genes and pv_genes 6/24/23
  ad_genes = c(
    "TWIST1",
    "LGALS1",
    "IL32",
    "CAPG",
    "ITM2C",
    "MFHAS1",
    "ANXA1",
    "SOS1",
    "CSGALNACT1",
    "LMO4",
    "IFITM2",
    "S100A10",
    "SYNE2",
    "THADA",
    "NEAT1",
    "IL17RB",
    "ARHGAP21",
    "NBAS",
    "ACTG1",
    "TGFBR3",
    "TNFSF10",
    "AHNAK",
    "ISG15",
    "RPL17",
    "LONRF2",
    "TSHZ2",
    "MMP25",
    "IFITM1",
    "BIRC3",
    "FAM102A",
    "LPCAT2",
    "NRIP3",
    "CRIP1",
    "CLU",
    "ZFP36",
    "ZFP36L2",
    "TUBA1B",
    "GATA3",
    "SLC5A3",
    "SFXN1",
    "FANK1",
    "TAGLN2",
    "RPS17",
    "SMKR1",
    "TC2N",
    "MYL12A",
    "LINC02195",
    "CRNDE",
    "CHDH",
    "HLA-DRB1",
    "SLC4A7",
    "MIB1",
    "SEC61G",
    "RPS29",
    "TRERF1",
    "S100A4"
  )

  pv_genes = c(
    "CXCL13",
    "CD7",
    "MGAT4A",
    "FTH1",
    "LAYN",
    "IL17F",
    "KLRB1",
    "GNLY",
    "CPM",
    "CTSH",
    "GBP5",
    "SOX4",
    "CLEC2B",
    "GZMB",
    "CD2",
    "CEBPD",
    "ODF2L",
    "LAG3",
    "LRRN3",
    "ARHGEF12",
    "PTPN13",
    "TNFAIP3",
    "TRPS1",
    "SNX9",
    "METRNL",
    "BTG1",
    "JUN",
    "SPOCK2",
    "GABARAPL1",
    "PMEPA1",
    "RBPJ",
    "LINC01871",
    "MAP3K4",
    "UBC",
    "GALNT1",
    "PNRC1",
    "GABPB1-AS1",
    "RPS26",
    "MUC20-OT1",
    "CHN1",
    "NAP1L4",
    "PTMS",
    "F2R",
    "CTLA4",
    "DAPK2",
    "RAP1B",
    "CCR6",
    "B3GALT2",
    "YPEL2",
    "FYN",
    "PPDPF",
    "SLA2",
    "CBLB",
    "ADGRG1",
    "SARAF",
    "DUSP1"
  )


  # RashX DEGs
  set1 = list()
  tryCatch({
    gc()
    set1[[paste0("RashX")]] <-
      FindMarkers(
        Unknown,
        ident.1 = "RashX",
        ident.2 = c("150", "154", "155", "169", "195", "204", "207"),
        group.by = "donor",
        verbose = TRUE,
        assay = "RNA",
        slot = "data",
        test.use = "wilcox",
        min.cells.feature = 0,
        min.cells.group = 0,
        logfc.threshold = 0,
        min.pct = 0,
        min.diff.pct = 0,
        features = genes
      )
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })

  tryCatch({ #Start error catch block

    # ref DEGs
    marker <-
      FindMarkers(
        skin.integrated2,
        ident.1 = "AD1",
        ident.2 = "Pso",
        group.by = "dis",
        verbose = FALSE,
        assay = "RNA",
        slot = "data",
        test.use = "wilcox",
        min.cells.feature = 0,
        min.cells.group = 0,
        logfc.threshold = 0,
        min.pct = 0,
        min.diff.pct = 0,
        features = genes
      )


    # Add a column to split AD and PV

    marker$group <- ifelse(rownames(marker) %in% ad_genes, "AD", "PV")
    set1[["ADvsPso"]] <- marker

    library(plyr)

    for (i in 1:length(set1)) {
      colnames(set1[[i]]) <-
        paste0(names(set1)[i], "_", colnames(set1[[i]]))
      set1[[i]]$ROWNAMES <- rownames(set1[[i]])
    }

    data <- join_all(set1, by = "ROWNAMES", type = "full")
    rownames(data) <- data$ROWNAMES
    data$ROWNAMES <- NULL

    RashX1 <- subset(data, ADvsPso_group == "AD")
    RashX1 <- RashX1[, which(str_detect(colnames(RashX1), "_log2FC"))]
    RashX2 <- subset(data, ADvsPso_group == "PV")
    RashX2 <- RashX2[, which(str_detect(colnames(RashX2), "_log2FC"))]
    RashX2$ADvsPso_avg_log2FC = RashX2$ADvsPso_avg_log2FC * (-1)

    # make heatmap figures
    #AD-specific genes heatmap, user can downlaod it
    png(
      plot1_name,
      width = 3,
      height = 7,
      units = "in",
      res = 1200
    )
    colnames(RashX1) <- c("RashX", "ADvsPV")
    ht1 <- Heatmap(
      as.matrix(RashX1),
      name = "RashX",
      #### annotate legend
      col = colorRamp2(
        c(-2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2),
        c(
          "#2D004B",
          "#542788",
          "#8073AC",
          "#B2ABD2",
          "#D8DAEB",
          "#F7F7F7",
          "#FEE0B6",
          "#FDB863",
          "#E08214",
          "#B35806",
          "#7F3B08"
        )
      ),
      #### set the color scales
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      #cluster_columns = F,
      clustering_distance_columns = "canberra",
      #canberra, euclidean
      column_title = "AD-specific genes",
      show_column_names = T,
      show_row_names = T,
      row_dend_side = "left",
      column_dend_side = "top",
      column_title_gp = gpar(fontsize = 12)
    )


    draw(ht1, auto_adjust = FALSE)
    dev.off()
    #Pso-specific genes heatmap, user can downlaod it
    png(
      plot2_name,
      width = 3,
      height = 7,
      units = "in",
      res = 1200
    )
    colnames(RashX2) <- c("RashX", "PVvsAD")
    ht2 <- Heatmap(
      as.matrix(RashX2),
      name = "RashX",
      #### annotate legend
      col = colorRamp2(
        c(-2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2),
        c(
          "#7F3B08",
          "#B35806",
          "#E08214",
          "#FDB863",
          "#FEE0B6",
          "#F7F7F7",
          "#D8DAEB",
          "#B2ABD2",
          "#8073AC",
          "#542788",
          "#2D004B"
        )
      ),
      #### set the color scales
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      #cluster_columns = F,
      clustering_distance_columns = "canberra",
      column_title = "PV-specific genes",
      show_column_names = T,
      show_row_names = T,
      row_dend_side = "left",
      column_dend_side = "top",
      column_title_gp = gpar(fontsize = 12)
    )

    draw(ht2, auto_adjust = FALSE)
    dev.off()

  }, error=function(e){print("ERROR: DEG analysis between AD and PV failed.")}) #End error catch block and print out error message


  tryCatch({ #Start error catch block

    Idents(Unknown) <- Unknown$donor
    Unknown <-
      RenameIdents(
        Unknown,
        `150` = "NML1",
        `154` = "NML2",
        `155` = "NML3",
        `169` = "NML4",
        `195` = "NML5",
        `204` = "NML6",
        `207` = "NML7",
        `170` = "AD1",
        `198` = "AD2",
        `230` = "AD3",
        `231` = "AD4",
        `232` = "AD5",
        `233` = "AD6",
        `236` = "AD7",
        `165` = "Pso1",
        `173` = "Pso2",
        `194` = "Pso3",
        `199` = "Pso4",
        `211` = "Pso5",
        `222` = "Pso6",
        `234` = "Pso7",
        `235` = "Pso8",
        `RashX` = "RashX"
      )

    Unknown[["STATUS"]] <- Idents(object = Unknown)

    #Seurat object --> CDS
    exp_mat <-
      Unknown@assays[["RNA"]]@data #pull out NORMALIZED counts from Seurat object
    cell_metadat <-
      Unknown@meta.data #pull out cell meta-data from Seurat object
    gene_annot = data.frame(Unknown@assays[["RNA"]]@counts@Dimnames[[1]]) #pull out gene names from Seurat object
    names(gene_annot) = "gene_short_name"
    row.names(gene_annot) = gene_annot$gene_short_name #row.names of gene_metadata must be equal to row.names of expression_data

    data.frame(names(cell_metadat))

    human_human_big_cds <- new_cell_data_set(exp_mat,
                                             cell_metadata = cell_metadat,
                                             gene_metadata = gene_annot)

    #create 1 large object from all the cds's created above

    #Subset to patients of interest
    pts = c("RashX", "AD1", "Pso")
    human_human_big_cds = human_human_big_cds[, colData(human_human_big_cds)$dis %in% pts] #subset to all normal patients and a diseased patient
    colData(human_human_big_cds)$dis_updated = colData(human_human_big_cds)$donor
    colData(human_human_big_cds)$dis_updated = gsub("170", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("198", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("230", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("231", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("232", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("233", "AD", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("236", "AD", colData(human_human_big_cds)$dis_updated)

    colData(human_human_big_cds)$dis_updated = gsub("165", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("173", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("194", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("199", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("211", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("222", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("234", "PV", colData(human_human_big_cds)$dis_updated)
    colData(human_human_big_cds)$dis_updated = gsub("235", "PV", colData(human_human_big_cds)$dis_updated)

    colData(human_human_big_cds)$dis_updated = gsub("RashX", "RashX", colData(human_human_big_cds)$dis_updated)


    #--------------------------------------------------------------------#
    #Get the genes of interest
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    ad_mat = t(human_human_big_cds@assays@data@listData[["counts"]][ad_genes,])
    ad_mat[1:10, 1:10]
    ad_int = data.frame(
      colData(human_human_big_cds)$dis_updated,
      colData(human_human_big_cds)$donor,
      rowSums(ad_mat),
      colData(human_human_big_cds)$STATUS
    )
    names(ad_int) = c("dis", "sample", "gene_sig", "status")
    names(ad_int) = c("dis", "sample", "gene_sig", "status")
    ad_dfc = summarySE(ad_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "status"))
    head(ad_dfc)
    names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")

    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = t(human_human_big_cds@assays@data@listData[["counts"]][pv_genes,])
    pv_mat[1:10, 1:10]
    pv_int = data.frame(
      colData(human_human_big_cds)$dis_updated,
      colData(human_human_big_cds)$donor,
      rowSums(pv_mat),
      colData(human_human_big_cds)$STATUS
    )
    names(pv_int) = c("dis", "sample", "gene_sig", "status")
    names(pv_int) = c("dis", "sample", "gene_sig", "status")
    pv_dfc = summarySE(pv_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "status"))
    head(pv_dfc)
    names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")

    #--------------------------------------------------------------------#
    re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures

    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#

    p1 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
      geom_point(alpha = 1, aes(color = dis), size = 2) +
      geom_text_repel(
        data = re_int,
        aes(x = ad_gene_sig, y = pv_gene_sig, label = status),
        box.padding = 0.4
      ) + #add labels
      geom_errorbarh(aes(
        xmax = ad_gene_sig + ad_se,
        xmin = ad_gene_sig - ad_se,
        color = dis
      )) +
      geom_errorbar(aes(
        ymax = pv_gene_sig + pv_se,
        ymin = pv_gene_sig - pv_se,
        color = dis
      )) +
      #stat_ellipse(aes(color=dis)) +
      theme_classic() +
      xlab("AD-specific genes") +
      ylab("PV-specific genes") +
      theme(legend.position = "right") +
      scale_color_manual(
        name = "",
        values = c(
          "AD" = "darkblue",
          "PV" = "darkred",
          "RashX" = "darkorange"
        )
      ) +
      scale_fill_manual(
        name = "",
        values = c(
          "AD" = "darkblue",
          "PV" = "darkred",
          "RashX" = "darkorange"
        )
      ) +
      geom_mark_hull(aes(fill = dis, label = dis), concavity = 5) +
      theme(
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        plot.title = element_text(size = 15, color = 'black')
      )

    ggsave(plot = p1,
           plot3_name,
           width = 6.5,
           height = 5)

    #----------------------------------------------------------------#
    ### Plot: Sample-level means with individual cells
    #----------------------------------------------------------------#

    #make df of individual cells
    int_dat = data.frame(ad_int, pv_int)

    #Calculate disease centroids
    ad_centroid = data.frame(mean(subset(re_int, dis == "AD")$ad_gene_sig), mean(subset(re_int, dis ==
                                                                                          "AD")$pv_gene_sig))
    pv_centroid = data.frame(mean(subset(re_int, dis == "PV")$ad_gene_sig), mean(subset(re_int, dis ==
                                                                                          "PV")$pv_gene_sig))
    names(ad_centroid) = c("xpt", "ypt")
    names(pv_centroid) = c("xpt", "ypt")

    p2 <- ggplot() +
      geom_point(
        data = int_dat,
        alpha = 0.25,
        aes(x = gene_sig, y = gene_sig.1, color = dis),
        size = 0.4
      ) + #individual cells
      geom_point(
        data = re_int,
        alpha = 1,
        aes(
          x = ad_gene_sig,
          y = pv_gene_sig,
          color = dis,
          shape = dis
        ),
        size = 4
      ) + #sample centroids
      #geom_errorbarh(data=re_int, aes(x=pc1,xmax = pc1 + pc1_se, xmin = pc1 - pc1_se, color=dis_updated)) +
      #geom_errorbar(data=re_int, aes(y=pc2,ymax = pc2 + pc2_se, ymin = pc2 - pc2_se, color=dis_updated)) +
      #stat_ellipse(data=re_int, aes(x=ad_gene_sig, y=pv_gene_sig, color=dis)) +
      stat_ellipse(data = int_dat, aes(x = gene_sig, y = gene_sig.1, color =
                                         dis)) +
      theme_classic() +
      ggtitle("Optimized gene program, individual cells") +
      xlab("AD program") +
      ylab("PV program") +
      theme(legend.position = "right") +
      scale_color_manual(
        name = "",
        values = c(
          "AD" = "darkblue",
          "PV" = "darkred",
          "RashX" = "darkorange"
        )
      ) +
      scale_fill_manual(
        name = "",
        values = c(
          "AD" = "darkblue",
          "PV" = "darkred",
          "RashX" = "darkorange"
        )
      ) +
      geom_point(
        data = ad_centroid,
        aes(x = xpt, y = ypt),
        color = "darkblue",
        shape = 10,
        size = 7,
        stroke = 2
      ) +
      geom_point(
        data = pv_centroid,
        aes(x = xpt, y = ypt),
        color = "darkred",
        shape = 10,
        size = 7,
        stroke = 2
      ) +
      theme(
        axis.text = element_text(size = 20, color = 'black'),
        axis.title = element_text(size = 20, color = 'black'),
        plot.title = element_text(size = 20, color = 'black')
      )


    ggsave(plot = p2,
           plot4_name,
           width = 6.5,
           height = 5)

  }, error=function(e){print("ERROR: Hyperdimensionality plot failed.")})


  #----------------------------------------------------------------#
  # Statistic analysis
  #----------------------------------------------------------------#

  tryCatch({ #Start error catch block


    all_coords_mat = as.matrix(data.frame(re_int$ad_gene_sig, re_int$pv_gene_sig))
    row.names(all_coords_mat) = re_int$sstatus

    dist_mat2 = vegdist(
      all_coords_mat,
      method = "canberra",
      diag = FALSE,
      upper = TRUE,
      p = 2
    )
    dist_mat2 <- as.matrix(dist_mat2)
    rownames(dist_mat2) = re_int$status
    colnames(dist_mat2) = re_int$status

    ind_samps = subset(re_int, dis == "RashX")$status
    ind_samps <- as.character(ind_samps)

    statistic =
      #Loop through the indeterminate samples and get averages
      all_res = do.call(rbind, lapply(1:length(ind_samps), function(i) {
        #i=1
        dist_df = data.frame(re_int$dis, re_int$status, dist_mat2[ind_samps[i],]) #get distance of a sample from the matrix and combine with metadata
        names(dist_df) = c("dis", "status", "dist")

        ad_dist_df = subset(dist_df, dis == "AD") #remove other indeterminate samples
        pv_dist_df = subset(dist_df, dis == "PV") #remove other indeterminate samples

        #calculate means to see which sided test to use
        ad_mean = mean(ad_dist_df$dist)
        pv_mean = mean(pv_dist_df$dist)

        #now test based on mean direction, find which is least
        if (ad_mean > pv_mean) {
          wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative = "greater")
          wilrest = data.frame("PV",
                               ad_mean,
                               pv_mean,
                               ad_mean - pv_mean,
                               wiltest$statistic,
                               wiltest$p.value)
        } else {
          wiltest = wilcox.test(ad_dist_df$dist, pv_dist_df$dist, alternative = "less")
          wilrest = data.frame("AD",
                               ad_mean,
                               pv_mean,
                               ad_mean - pv_mean,
                               wiltest$statistic,
                               wiltest$p.value)
        }

        #comine with data and return
        res = data.frame(ind_samps[i], wilrest)
        names(res) = c("sample",
                       "proximity",
                       "AD_dist_mean",
                       "PV_dist_mean",
                       "AD_PV",
                       "W",
                       "p")
        res
      }))


    #### visualize the results directly
    statistic

    ###save the statistic results as csv format, user can downlaod it
    write.csv(statistic, exported_file_name)

    return("Processed")

  }, error=function(e){print("ERROR: Hyperdimensional proximity testing failed.")}) #End error catch block and print out error message



}
