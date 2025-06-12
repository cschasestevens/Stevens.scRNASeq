#' Create Seurat Object from QC Dataframe
#'
#' Converts a filtered phenocycler data frame from pc_qc()
#' to a Seurat object.
#'
#' @param df Input data frame.
#' @param dir1 Directory name for saving Seurat objects.
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
  dir1 = "data/",
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
  parallel::mclapply(
    mc.cores = mcc + 2,
    seq.int(1, length(d3), 1),
    function(x) {
      saveRDS(
        d3[[x]],
        paste(dir1, "d_seurat_", names(d3)[[x]], ".rds", sep = "")
      )
    }
  )
  remove(d2)
  gc(reset = TRUE)
  return(print(paste(".rds file created for ", names(d3), sep = "")))
}

#' Create AnnData from Seurat
#'
#' Converts a Seurat Object to a h5ad file.
#'
#' @param type1 Data type; provide a character string indicating
#' if the data have been normalized.
#' @param so_list Input Seurat object list.
#' @param asy Assay to use for file conversion.
#' @param mcc Cores to use for file conversion.
#' @return A h5ad file with the same name as the input object.
#' @examples
#'
#' # ld_seurat <- pc_create_ad(
#' #   so = readRDS("data/seuratobject.rds")
#' # )
#'
#' @export
pc_create_ad <- function(
  type1 = "norm",
  so_list,
  asy = "PC",
  mcc = 8
) {
  ifelse(
    dir.exists(paste("data/", sep = "")) == FALSE,
    dir.create(paste("data/", sep = "")),
    print(
      paste("Saving QC cellmap in data/", sep = "")
    )
  )
  # Load data
  d3 <- so_list
  # Convert Seurat objects to individual AnnData files
  parallel::mclapply(
    mc.cores = mcc,
    seq.int(1, length(d3), 1),
    function(x) {
      d2 <- d3[[x]]
      d1 <- Seurat::CreateAssayObject(counts = d2[[asy]]$counts)
      d2[["PCv3"]] <- d1
      Seurat::DefaultAssay(d2) <- "PCv3"
      d2[[asy]] <- NULL
      d2 <- SeuratObject::RenameAssays(object = d2, PCv3 = asy)
      SeuratDisk::SaveH5Seurat(
        d2,
        filename = paste(
          "data/d_", names(d3)[[x]], "_", type1, ".h5Seurat", sep = ""
        )
      )
      SeuratDisk::Convert(
        paste("data/d_", names(d3)[[x]], "_", type1, ".h5Seurat", sep = ""),
        dest = "h5ad"
      )
      file.remove(paste(
        "data/d_", names(d3)[[x]], "_", type1, ".h5Seurat", sep = ""
      ))
      return(d2) # nolint
    }
  )
}
