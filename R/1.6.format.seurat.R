#' Create Seurat Object from QC Dataframe
#'
#' Converts a filtered phenocycler data frame from pc_qc()
#' to a Seurat object.
#'
#' @param df Input data frame.
#' @param samp_id Column containing sample IDs for selected data.
#' @param mcc Cores to use for processing input data.
#' @return List of Seurat objects containing filtered phenocycler data.
#' @examples
#'
#' # ld_seurat <- pc_create_seurat(
#' #   df = pheno_qc[["data_ind"]]
#' # )
#'
#' @export
pc_create_seurat <- function(
  df,
  samp_id = "Code",
  mcc = 2
) {
  # Load data
  d1 <- df
  sid <- samp_id
  # Split data by individual sample
  d1 <- setNames(
    parallel::mclapply(
      mc.cores = mcc,
      seq.int(1, length(unique(d1[[sid]])), 1),
      function(x) {
        d2 <- d1[d1[[sid]] == unique(d1[[sid]])[[x]], ]
        return(d2) # nolint
      }
    ),
    unique(d1[[sid]])
  )
  lapply(d1, function(x) head(x))
  gc(reset = TRUE)
  # Cast data frames
  d2 <- setNames(
    parallel::mclapply(
      mc.cores = mcc,
      seq.int(1, length(d1), 1),
      function(x) {
        d1 <- d1[[x]]
        d2 <- reshape2::dcast(
          d1,
          ID + Slide + Code + Group + X + Y ~ Channel,
          value.var = "value"
        )
        return(d2) # nolint
      }
    ),
    names(d1)
  )
  lapply(d2, function(x) head(x))
  remove(d1)
  gc(reset = TRUE)
  # Convert each data frame into Seurat objects
  d3 <- setNames(
    parallel::mclapply(
      mc.cores = mcc + 2,
      seq.int(1, length(d2), 1),
      function(x) {
        d1 <- Seurat::CreateSeuratObject(
          counts = t(as.matrix(d2[[x]][, 7:ncol(d2[[x]])])),
          meta.data = d2[[x]][, 1:6],
          assay = "PC"
        )
        return(d1) # nolint
      }
    ),
    names(d2)
  )
  d3
  parallel::mclapply(
    mc.cores = mcc + 2,
    seq.int(1, length(d3), 1),
    function(x) {
      saveRDS(
        d3[[x]],
        paste("data/d_seurat_", names(d3)[[x]], ".rds", sep = "")
      )
    }
  )
  remove(d2)
  gc(reset = TRUE)
  # Convert Seurat objects to individual AnnData files
  parallel::mclapply(
    mc.cores = 2 + 6,
    seq.int(1, length(d3), 1),
    function(x) {
      d2 <- d3[[x]]
      d1 <- Seurat::CreateAssayObject(counts = d2[["PC"]]$counts)
      d2[["PCv3"]] <- d1
      Seurat::DefaultAssay(d2) <- "PCv3"
      d2[["PC"]] <- NULL
      d2 <- SeuratObject::RenameAssays(object = d2, PCv3 = "PC")
      SeuratDisk::SaveH5Seurat(
        d2,
        filename = paste("data/pc_", names(d3)[[x]], ".h5Seurat", sep = "")
      )
      SeuratDisk::Convert(
        paste("data/pc_", names(d3)[[x]], ".h5Seurat", sep = ""),
        dest = "h5ad"
      )
      file.remove(paste("data/pc_", names(d3)[[x]], ".h5Seurat", sep = ""))
      return(d2) # nolint
    }
  )
  return(print(paste(".h5ad file created for ", names(d3), sep = "")))
}
