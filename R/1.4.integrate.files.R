#' Batch Integration of scRNA-Seq Datasets
#'
#' Integrates a list of samples processed as Seurat objects.
#'
#' @param l.so List of processed Seurat objects to be integrated.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Proportion (as a numeric value) of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A Seurat object containing integrated data for all samples present in a scRNA-Seq experiment.
#' @examples
#'
#' # d.integrated <- sc.integrate.data(
#' #    list.data,
#' #    TRUE,
#' #    0.5
#' #    )
#'
#' @export
sc.integrate.data <- function(
    l.so,
    parl,
    core.perc
    ){
  list.d <- l.so
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




