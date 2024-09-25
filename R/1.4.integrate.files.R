#' Batch Integration of scRNA-Seq and Multiome Datasets
#'
#' Integrates a list of samples processed as Seurat objects.
#'
#' @param l_so List of processed Seurat objects to be integrated.
#' @param l_par List of processing parameters passed to function.
#' For multiome data integration, must have a gene annotation and
#' reference genome file.
#' @param proc_mode Processing mode to be used
#' (either "scRNA-Seq" or "multiome").
#' @param parl Logical indicating whether processing should be run
#' in parallel. Set to FALSE if running sequentially.
#' @param core_perc Proportion (as a numeric value) of available cores
#' to use if running in parallel. Set to 1 if running sequentially.
#' @return A Seurat object containing integrated data for all samples
#' present in a scRNA-Seq or multiome experiment.
#' @examples
#'
#' # d.integrated <- sc.integrate.data(
#' #   list_data,
#' #   l_params,
#' #   "multiome",
#' #   TRUE,
#' #   0.5
#' # )
#'
#' @export
sc_integrate_data <- function(
  l_so,
  l_par,
  proc_mode,
  parl,
  core_perc
) {

  if(proc_mode == "multiome") { # nolint
    if(unlist(packageVersion("Seurat"))[1] == 5) { # nolint
      if(file.exists(""))
      ld_int <- l_so
      # Extract Seurat objects for integration
      ld_int <- setNames(
        lapply(
          seq.int(1, length(ld_int), 1),
          function(x) {
            d1 <- ld_int[[x]][[1]]
            Seurat::DefaultAssay(d1) <- "SCT"
            d1 <- Seurat::RunPCA(d1)
            Seurat::DefaultAssay(d1) <- "peaks"
            return(d1)
          }
        ),
        names(ld_int)
      )
      # Unify ATAC peaks across data sets
      ufy_peaks <- Signac::UnifyPeaks(
        object.list = ld_int,
        mode = "reduce"
      )
      ## Filter low quality peaks
      ufy_width <- BSgenome::width(ufy_peaks)
      range(ufy_width)
      quantile(ufy_width)
      median(ufy_width)
      # Create new peak assays using unified list
      ld_int <- setNames(
        parallel::mclapply(
          mc.cores = ceiling(
            parallel::detectCores() *
              core_perc
          ),
          seq.int(1, length(ld_int), 1),
          function(x) {
            d1 <- ld_int[[x]]
            Seurat::DefaultAssay(d1) <- "peaks"
            d1_cnts <- Signac::FeatureMatrix(
              fragments = Signac::Fragments(d1),
              features = ufy_peaks,
              cells = colnames(d1)
            )
            d1[["ufy.peaks"]] <- Signac::CreateChromatinAssay(
              counts = d1_cnts,
              fragments = Signac::Fragments(d1),
              annotation = l_par[["ref.gtf"]]
            )
            Seurat::DefaultAssay(d1) <- "ufy.peaks"
            return(d1)
          }
        ),
        names(ld_int)
      )
      saveRDS(ld_int, "processed/data.processed.unified.rds")

      # Compute WNN for each sample, then integrate layers
      ## Run Harmony to correct for batch effects introduced by Code
      ### Remove previous objects prior to integration steps
      remove(list = ls())
      gc(reset = TRUE)

      # Weighted nearest neighbors (Integrate RNA and unified ATAC peak assays)
      ld_int <- readRDS("processed/data.processed.unified.rds")
      ld_int <- setNames(
        parallel::mclapply(
          mc.cores = ceiling(
            parallel::detectCores() *
              0.1
          ),
          seq.int(1, length(ld_int), 1),
          function(x) {
            d1 <- ld_int[[x]]
            d1 <- Seurat::FindMultiModalNeighbors(
              d1,
              reduction.list = list("pca", "lsi"),
              dims.list = list(1:50, 2:50)
            )
            Seurat::DefaultAssay(d1) <- "ufy.peaks"
            return(d1)
          }
        ),
        names(ld_int)
      )
      gc(reset = TRUE)
      library(Seurat)
      d_int <- merge(
        x = ld_int[[1]],
        y = ld_int[2:length(ld_int)],
        add.cell.ids = names(ld_int)
      )
      saveRDS(d_int, "processed/data.processed.merge.rds")
      remove(list = ls())
      gc(reset = TRUE)
      d_int <- readRDS("processed/data.processed.merge.rds")
        



        # Normalization
        ## GEX
        Seurat::DefaultAssay(d_int) <- "RNA"
        d_int <- Seurat::FindVariableFeatures(
          d_int,
          selection.method = "vst",
          nfeatures = 3000
        )
        d_int <- Seurat::NormalizeData(d_int)
        d_int <- Seurat::ScaleData(d_int)
        d_int <- Seurat::RunPCA(d_int)
        
        ## ATAC Peaks
        Seurat::DefaultAssay(d_int) <- "ufy.peaks"
        d_int <- Signac::FindTopFeatures(d_int, min.cutoff = 5)
        d_int <- Signac::RunTFIDF(d_int)
        d_int <- Signac::RunSVD(d_int)
        Signac::DepthCor(d_int)
        
        ## Integrate Layers
        library(Seurat)
        library(SeuratWrappers)
        Seurat::DefaultAssay(d_int) <- "RNA"
        
        d_int <- Seurat::IntegrateLayers(
          object = d_int,
          method = CCAIntegration,
          orig.reduction = "pca",
          new.reduction = "int.SCT.CCA",
          verbose = TRUE
        )
        d_int <- SeuratObject::JoinLayers(d_int)
        
        ## Harmony batch effect correction for GEX and ATAC
        d_int <- harmony::RunHarmony(
          d_int,
          assay.use = "ufy.peaks",
          group.by.vars = "Code",
          reduction.use = "lsi",
          reduction.save = "ufy.peaks.corrected",
          project.dim = FALSE
        )
        
        d_int <- harmony::RunHarmony(
          d_int,
          assay.use = "RNA",
          group.by.vars = "Code",
          reduction.use = "int.SCT.CCA",
          reduction.save = "RNA.corrected",
          project.dim = FALSE
        )
        
        ## WNN
        d_int <- Seurat::FindMultiModalNeighbors(
          d_int,
          reduction.list = list("RNA.corrected", "ufy.peaks.corrected"),
          dims.list = list(1:50, 2:50)
        )
        
        d_int <- Seurat::RunUMAP(
          d_int,
          nn.name = "weighted.nn",
          reduction.name = "wnn.umap",
          reduction.key = "wnnUMAP_"
        )
        
        d_int <- Seurat::FindClusters(
          d_int,
          graph.name = "wsnn",
          algorithm = 3,
          verbose = TRUE,
          resolution = 0.5
        )
        head(d_int@meta.data)
        d_int
        
        ## Visualize Clusters
        ggplot2::ggsave(
          "analysis/plot.umap.panel.WNN.png",
          sc_umap_panel(d_int,
            c("Group", "Code", "seurat_clusters"),
            "wnn.umap"
          ),
          height = 8,
          width = 24,
          dpi = 600
        )
        
        ## Integration Performance
        d1qc_pre <- Seurat::VlnPlot(
          object = d_int,
          features = c(
            "nCount_RNA",
            "nCount_ATAC",
            "TSS.enrichment",
            "nucleosome_signal"
          ),
          layer = "counts",
          ncol = 4,
          pt.size = 0.2
        )
        
        ggplot2::ggsave(
          "analysis/plot.integration.qc.png",
          d1qc_pre,
          width = 24,
          height = 8,
          dpi = 800
        )
        
        
        ## Save
        saveRDS(d_int, "analysis/data.integrated.rds")
    }
  }

  list_d <- l_so
  
  if(Sys.info()[["sysname"]] != "Windows" &
     parl == TRUE) {
    fun.int <- function(list.d.subset,future.size){
      ## Find integration anchors
      future::plan("multisession",workers = ceiling(parallel::detectCores()*core.perc))
      options(future.globals.maxSize = future.size * 1024^2)
      d.anchor <- Seurat::FindIntegrationAnchors(
        object.list = list.d.subset,
        anchor.features = 4000,
        dims = 1:50
        )
      ## Integrate data
      future::plan("sequential")
      options(future.globals.maxSize = 500 * 1024^2)
      d.merged <- Seurat::IntegrateData(
        anchorset = d.anchor,
        dims = 1:50
      )
      return(d.merged)
      }

    
    
    if(length(list.d) <= 12 & length(list.d) > 9) {
      print(paste(length(list.d),"samples present: dividing integration into 4 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:9],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
        )
      # Batch 4
      d.merged4 <- fun.int(list.d[10:length(list.d)],5000)
      saveRDS(
        d.merged4,
        "processed/data.integrated.batch.4.rds"
        )
      ## Combine batch 1&2 then 3&4
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
        )
      d.merged34 <- fun.int(list(d.merged3,d.merged4),10000)
      saveRDS(
        d.merged34,
        "processed/data.integrated.batch.3and4.rds"
        )
      remove(d.merged1,d.merged2,d.merged3,d.merged4)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged34),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged12,d.merged34)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.4.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
      file.remove("processed/data.integrated.batch.3and4.rds")
      }

    if(length(list.d) <= 9 & length(list.d) > 6) {
      print(paste(length(list.d),"samples present: dividing integration into 3 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:length(list.d)],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
        )
      ## Combine batch 1&2 then 3
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
        )
      remove(d.merged1,d.merged2)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged3),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged12,d.merged3)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
      }

    if(length(list.d) <= 6 & length(list.d) > 3) {
      print(paste(length(list.d),"samples present: dividing integration into 2 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      ## Combine batch 1&2
      d.integrated <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged1,d.merged2)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      }

    if(length(list.d) <= 3) {
      print(paste(length(list.d),"samples present: integrating all samples as 1 batch...",sep = " "))
      # Batch 1
      d.integrated <- fun.int(list.d[1:length(list.d)],5000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      }
    }

  if(Sys.info()[["sysname"]] != "Windows" &
     parl == FALSE) {
    fun.int <- function(list.d.subset,future.size){
      ## Find integration anchors
      d.anchor <- Seurat::FindIntegrationAnchors(
        object.list = list.d.subset,
        anchor.features = 4000,
        dims = 1:50
        )
      ## Integrate data
      d.merged <- Seurat::IntegrateData(
        anchorset = d.anchor,
        dims = 1:50
        )
      return(d.merged)
      }

    if(length(list.d) <= 12 & length(list.d) > 9) {
      print(paste(length(list.d),"samples present: dividing integration into 4 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:9],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
        )
      # Batch 4
      d.merged4 <- fun.int(list.d[10:length(list.d)],5000)
      saveRDS(
        d.merged4,
        "processed/data.integrated.batch.4.rds"
        )
      ## Combine batch 1&2 then 3&4
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
        )
      d.merged34 <- fun.int(list(d.merged3,d.merged4),10000)
      saveRDS(
        d.merged34,
        "processed/data.integrated.batch.3and4.rds"
        )
      remove(d.merged1,d.merged2,d.merged3,d.merged4)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged34),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged12,d.merged34)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.4.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
      file.remove("processed/data.integrated.batch.3and4.rds")
    }

    if(length(list.d) <= 9 & length(list.d) > 6) {
      print(paste(length(list.d),"samples present: dividing integration into 3 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:length(list.d)],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
        )
      ## Combine batch 1&2 then 3
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
        )
      remove(d.merged1,d.merged2)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged3),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged12,d.merged3)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
    }

    if(length(list.d) <= 6 & length(list.d) > 3) {
      print(paste(length(list.d),"samples present: dividing integration into 2 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
        )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
        )
      ## Combine batch 1&2
      d.integrated <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      remove(d.merged1,d.merged2)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      }

    if(length(list.d) <= 3) {
      print(paste(length(list.d),"samples present: integrating all samples as 1 batch...",sep = " "))
      # Batch 1
      d.integrated <- fun.int(list.d[1:length(list.d)],5000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
        )
      }
    }

  if(Sys.info()[["sysname"]] == "Windows") {
    fun.int <- function(list.d.subset,future.size){
      ## Find integration anchors
      d.anchor <- Seurat::FindIntegrationAnchors(
        object.list = list.d.subset,
        anchor.features = 4000,
        dims = 1:50
      )
      ## Integrate data
      d.merged <- Seurat::IntegrateData(
        anchorset = d.anchor,
        dims = 1:50
      )
      return(d.merged)
    }

    if(length(list.d) <= 12 & length(list.d) > 9) {
      print(paste(length(list.d),"samples present: dividing integration into 4 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
      )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
      )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:9],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
      )
      # Batch 4
      d.merged4 <- fun.int(list.d[10:length(list.d)],5000)
      saveRDS(
        d.merged4,
        "processed/data.integrated.batch.4.rds"
      )
      ## Combine batch 1&2 then 3&4
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
      )
      d.merged34 <- fun.int(list(d.merged3,d.merged4),10000)
      saveRDS(
        d.merged34,
        "processed/data.integrated.batch.3and4.rds"
      )
      remove(d.merged1,d.merged2,d.merged3,d.merged4)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged34),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
      )
      remove(d.merged12,d.merged34)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.4.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
      file.remove("processed/data.integrated.batch.3and4.rds")
    }

    if(length(list.d) <= 9 & length(list.d) > 6) {
      print(paste(length(list.d),"samples present: dividing integration into 3 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
      )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
      )
      # Batch 3
      d.merged3 <- fun.int(list.d[7:length(list.d)],5000)
      saveRDS(
        d.merged3,
        "processed/data.integrated.batch.3.rds"
      )
      ## Combine batch 1&2 then 3
      d.merged12 <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.merged12,
        "processed/data.integrated.batch.1and2.rds"
      )
      remove(d.merged1,d.merged2)
      ## Combine into final batch and save
      d.integrated <- fun.int(list(d.merged12,d.merged3),15000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
      )
      remove(d.merged12,d.merged3)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
      file.remove("processed/data.integrated.batch.3.rds")
      file.remove("processed/data.integrated.batch.1and2.rds")
    }

    if(length(list.d) <= 6 & length(list.d) > 3) {
      print(paste(length(list.d),"samples present: dividing integration into 2 batches...",sep = " "))
      # Batch 1
      d.merged1 <- fun.int(list.d[1:3],5000)
      saveRDS(
        d.merged1,
        "processed/data.integrated.batch.1.rds"
      )
      # Batch 2
      d.merged2 <- fun.int(list.d[4:6],5000)
      saveRDS(
        d.merged2,
        "processed/data.integrated.batch.2.rds"
      )
      ## Combine batch 1&2
      d.integrated <- fun.int(list(d.merged1,d.merged2),10000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
      )
      remove(d.merged1,d.merged2)
      file.remove("processed/data.integrated.batch.1.rds")
      file.remove("processed/data.integrated.batch.2.rds")
    }

    if(length(list.d) <= 3) {
      print(paste(length(list.d),"samples present: integrating all samples as 1 batch...",sep = " "))
      # Batch 1
      d.integrated <- fun.int(list.d[1:length(list.d)],5000)
      saveRDS(
        d.integrated,
        "processed/data.integrated.rds"
      )
    }
  }


  remove(list.d)
  return(d.integrated)
  }


#' scRNA-Seq Integration Quality Control
#'
#' Displays the feature number, average counts, and percentage of reads derived from the mitochondrial genome given an integrated Seurat object.
#'
#' @param so Integrated Seurat object.
#' @param cl.var Clustering variable provided as a character string (generally use "seurat_clusters").
#' @return A panel of violin plots providing QC measures for an integrated Seurat object.
#' @examples
#'
#' # d.integrated <- sc.integrate.data(d.integrated)
#'
#' @export
sc.integration.qc <- function(
    so,
    cl.var
    ) {
  d <- so
  df <- d@meta.data[,c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    cl.var)]
  p.qc <- ggpubr::ggarrange(
    plotlist = lapply(
      names(
        dplyr::select(
          df,
          -.data[[cl.var]]
          )
        ),
      function(x) {
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            x = .data[[cl.var]],
            y = .data[[x]],
            fill = .data[[cl.var]]
            )
          ) +
          ggplot2::scale_fill_manual(
            name = "Cluster",
            values = col.univ()
            ) +
          # Add violin plot and dotplot
          ggplot2::geom_violin(
            trim = F
            ) +
          ggplot2::geom_jitter(
            ggplot2::aes(
              alpha = 0.2
              ),
            shape = 16,
            size = 0.1,
            position = ggplot2::position_jitter(
              width = 0.4
              ),
            show.legend = F
            ) +
          # Add Theme
          sc.theme1() +
          ggplot2::labs(y = x) +
          ggplot2::theme(
            plot.margin = ggplot2::unit(c(0.1,0.1,0.1,0.1),"cm")
            )
        }
      ),
    common.legend = T,
    legend = "none",
    nrow = 1,
    ncol = 3
    )
  return(p.qc)
  }




