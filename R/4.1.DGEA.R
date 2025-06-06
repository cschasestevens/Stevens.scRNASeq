#' scRNA-Seq DGEA
#'
#' Performs DGEA per cell type of a Seurat Object.
#'
#' @param so An object of class Seurat. Must contain the columns
#' 'CellType' and 'CellGroup' in the metadata slot.
#' @param asy Character string providing the name of the assay
#' to use for differential analysis.
#' @param slot1 Character string selecting the slot from the
#' Seurat object to pull from the chosen assay. Either type "data"
#' if using expression data and "counts" if analyzing ATAC data.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param ct Character string vector containing the name(s) of up to
#' 2 cell type groups present in the CellGroup column.
#' The provided Seurat object is then subsetted to include
#' only clusters from the provided group names.
#' @param mast_comp Character string indicating the name of the MAST
#' group comparison for conducting DGEA. MAST names are comprised of the chosen
#' variable name and the leading factor level within that variable.
#' @param mast_name User-defined name of a DGEA comparison,
#' given as a character string.
#' @param form1 Formula to use for MAST generalized linear model.
#' @param parl Logical indicating whether processing should be run in parallel
#' (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if running
#' in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @param atac_type (optional) Character string providing the name of a cell
#' type for performing differential analysis. Only used for scATAC-Seq data.
#' @return A list of DGEA results per cell type for the chosen group comparison,
#' including genes missing fold changes and cell type DGEA results
#' with errors.
#' @examples
#'
#' # diff.output <- sc_diff(
#' #   # Seurat object
#' #   d,
#' #   # Assay
#' #   "ATAC",
#' #   # Slot
#' #   "counts",
#' #   # metadata column list
#' #   c("Code","Airway","CellType","nFeature_ATAC"),
#' #   # Cell type column name
#' #   "CellType",
#' #   # MAST comparison name
#' #   "AirwaySAE",
#' #   # MAST name (user-provided)
#' #   "SAE vs. LAE",
#' #   # Formula
#' #   as.formula(
#' #     paste(
#' #       "~","Airway","+",
#' #       "nFeature_ATAC",
#' #       sep = " "
#' #     )
#' #   ),
#' #   # run in parallel? (Set to FALSE if on Windows)
#' #   TRUE,
#' #   # core percentage to use
#' #   0.08,
#' #   # (optional) CellType to use if assay is "ATAC"
#' #   atac_type = x
#' # )
#'
#' @export
sc_diff <- function( # nolint
  so,
  asy = "sct",
  slot1 = "scale.data",
  md_list,
  ct = "CellType",
  mast_comp,
  mast_name,
  form1,
  parl = TRUE,
  core_perc = 0.2,
  atac_type = NULL
) {
  if(missing(atac_type) == TRUE) { # nolint
    # Seurat object
    d <- so
    # Metadata columns
    lc <- md_list
    # Cell type column
    c <- ct
  }
  if(missing(atac_type) == FALSE) { # nolint
    # Seurat object
    d <- so
    # Metadata columns
    lc <- md_list
    # Cell type column
    c <- ct
    ## Subset
    d <- d[, d@meta.data[[c]] == atac_type]
  }
  # MAST Comparison (combines column name and leading factor level for name)
  mc <- mast_comp
  # MAST Comparison name
  mn <- mast_name
  ## Input
  deg_mat <- as.matrix(
    SeuratObject::GetAssayData(d, slot = slot1, assay = asy)
  )
  if(missing(atac_type) == FALSE) { # nolint
    if(asy == "chromvar") { # nolint
      rownames(deg_mat) <- rownames(d@assays[[asy]])
    }
    if(asy == "ATAC") { # nolint
      rownames(deg_mat) <- paste(
        d@assays[["chromvar"]]@meta.features$nearestGene,
        seq.int(1, nrow(d@assays[[asy]]@meta.features), 1),
        sep = "."
      )
    }
  }
  deg_cols <- data.frame(
    d@meta.data[, lc]
  )
  ## Format input as DGEA object
  if(missing(atac_type) == FALSE) { # nolint
    dgea_sc <- MAST::FromMatrix(
      deg_mat,
      cData = deg_cols,
      check_sanity = FALSE
    )
  }
  if(missing(atac_type) == TRUE) { # nolint
    dgea_sc <- MAST::FromMatrix(
      deg_mat,
      cData = deg_cols
    )
  }
  dgea_celltype <- unique(
    as.character(SingleCellExperiment::colData(dgea_sc)[[c]])
  )
  list_dgea <- list("SCE" = dgea_sc, "CellType" = dgea_celltype)
  remove(d)
  if(Sys.info()[["sysname"]] != "Windows" && parl == TRUE) { #nolint
    list_dgea_res <- setNames(parallel::mclapply(
      mc.cores = ceiling(parallel::detectCores() * core_perc),
      seq.int(1, length(list_dgea[[2]]), 1),
      function(x) {
        tryCatch(
          {
            # Subset data
            s1 <- list_dgea[[1]][ , SingleCellExperiment::colData(list_dgea[[1]])[[c]] == list_dgea[[2]][[x]]] #nolint
            s1_sum <- rowSums(SummarizedExperiment::assay(s1) > 0)
            s1_sum[is.na(s1_sum)] <- 0
            s1 <- s1[s1_sum / ncol(s1) >= 0.05, ]
            ### create glm (generalized linear model for each variable)
            s1_fit <- MAST::zlm(
              formula = form1,
              s1,
              method = "glm",
              ebayes = FALSE,
              parallel = FALSE
            )
            ### Output DFs
            d_mast_sum_fun <- function(
              comp1,
              ct2,
              comp1_name
            ) {
              s1_res <- MAST::summary(
                s1_fit,
                doLRT = comp1,
                logFC = TRUE,
                parallel = FALSE
              )
              ### make dfs to display summary results by comp
              s1_dt <- reshape2::melt(
                dplyr::select(
                  dplyr::filter(
                    s1_res$datatable,
                    contrast == comp1 & # nolint
                      component != "S" # nolint
                  ),
                  -c("contrast")
                ),
                id.vars = c("primerid", "component")
              )
              s1_dt[["vars"]] <- paste(
                s1_dt$component,
                s1_dt$variable,
                sep = "."
              )
              s1_dt <- dplyr::select(
                dplyr::mutate(
                  reshape2::dcast(
                    dplyr::select(
                      dplyr::filter(
                        s1_dt,
                        vars != "logFC.Pr(>Chisq)" & # nolint
                          vars != "H.ci.hi" &
                          vars != "H.ci.lo" &
                          vars != "H.coef" &
                          vars != "H.z"
                      ),
                      -c(
                        "component",
                        "variable"
                      )
                    ),
                    primerid ~ vars
                  ),
                  "CellType" = ct2,
                  "Comparison" = comp1_name
                ),
                c(
                  "CellType", "Comparison", "primerid",
                  "logFC.coef", "H.Pr(>Chisq)", "C.Pr(>Chisq)",
                  "D.Pr(>Chisq)", everything() # nolint
                )
              )
              names(s1_dt) <- c(
                "CellType", "Comparison", "GENE",
                "logFC", "H.pval", "C.pval",
                "D.pval",
                names(
                  s1_dt[8:ncol(
                    s1_dt
                  )]
                )
              )
              return(s1_dt) # nolint
            }
            d1 <- d_mast_sum_fun(
              mc,
              list_dgea[[2]][[1]],
              mn
            )
            return(d1)
          },
          error = function(e) {
            print("Differential analysis unsuccessful for cell type...")
          }
        )
      }
    ), as.character(list_dgea[[2]]))
  }
  if(Sys.info()[["sysname"]] == "Windows" | parl == FALSE) { # nolint
    list_dgea_res <- setNames(
      lapply( # nolint
        seq.int(1, length(list_dgea[[2]]), 1),
        function(x) {
          tryCatch(
            {
              # Subset data
              s1 <- list_dgea[[1]][,
              SingleCellExperiment::colData( # nolint
                list_dgea[[1]]
              )[[c]] == list_dgea[[2]][[x]]]
              s1_sum <- rowSums(SummarizedExperiment::assay(s1) > 0) # nolint
              s1_sum[is.na(s1_sum)] <- 0
              s1 <- s1[s1_sum / ncol(s1) >= 0.05, ]
              ### create glm (generalized linear model for each variable)
              s1_fit <- MAST::zlm(
                formula = form1,
                s1,
                method = "glm",
                ebayes = FALSE,
                parallel = FALSE
              )
              ### Output DFs
              d_mast_sum_fun <- function(
                comp1,
                ct2,
                comp1_name
              ) {
                s1_res <- MAST::summary(
                  s1_fit,
                  doLRT = comp1,
                  logFC = TRUE,
                  parallel = FALSE
                )
                ### make dfs to display summary results by comp
                s1_dt <- reshape2::melt(
                  dplyr::select(
                    dplyr::filter(
                      s1_res$datatable,
                      contrast == comp1 & # nolint
                        component != "S" # nolint
                    ),
                    -c("contrast")
                  ),
                  id.vars = c("primerid", "component")
                )
                s1_dt[["vars"]] <- paste(
                  s1_dt$component,
                  s1_dt$variable,
                  sep = "."
                )
                s1_dt <- dplyr::select(
                  dplyr::mutate(
                    reshape2::dcast(
                      dplyr::select(
                        dplyr::filter(
                          s1_dt,
                          vars != "logFC.Pr(>Chisq)" & # nolint
                            vars != "H.ci.hi" &
                            vars != "H.ci.lo" &
                            vars != "H.coef" &
                            vars != "H.z"
                        ),
                        -c(
                          "component",
                          "variable"
                        )
                      ),
                      primerid ~ vars
                    ),
                    "CellType" = ct2,
                    "Comparison" = comp1_name
                  ),
                  c(
                    "CellType", "Comparison", "primerid",
                    "logFC.coef", "H.Pr(>Chisq)", "C.Pr(>Chisq)",
                    "D.Pr(>Chisq)", everything() # nolint
                  )
                )
                names(s1_dt) <- c(
                  "CellType", "Comparison", "GENE",
                  "logFC", "H.pval", "C.pval",
                  "D.pval",
                  names(
                    s1_dt[8:ncol(
                      s1_dt
                    )]
                  )
                )
                return(s1_dt) # nolint
              }
              d1 <- d_mast_sum_fun( # nolint
                mc,
                list_dgea[[2]][[x]],
                mn
              )
              return(d1) # nolint
            },
            error = function(e) {
              print("Differential analysis unsuccessful for cell type...")
            }
          )
        }
      ), as.character(list_dgea[[2]]))
  }
  # Combine results
  dgea_comb <- list_dgea_res[lengths(list_dgea_res) > 1]
  ## isolate DGEA results with errors
  dgea_error <- list_dgea_res[lengths(list_dgea_res) <= 1]
  ## return genes for each result with missing logFC
  dgea_miss <- lapply(
    dgea_comb,
    function(x) x[is.na(x[["logFC"]]), ]
  )
  dgea_res <- lapply(
    dgea_comb,
    function(x) x[!is.na(x[["logFC"]]), ]
  )
  dgea_sum <- list(
    "D.results" = dplyr::bind_rows(dgea_res),
    "D.missing" = dplyr::bind_rows(dgea_miss),
    "D.errors" = dplyr::bind_rows(dgea_error)
  )
  dgea_sum[["D.results"]][["H.qval"]] <- p.adjust(
    dgea_sum[["D.results"]][["H.pval"]],
    method = "BH"
  )
  dgea_sum[["D.results"]][["log2FC"]] <- log2(
    exp(dgea_sum[["D.results"]][["logFC"]])
  )
  dgea_sum[["D.results"]] <- dplyr::select(
    dgea_sum[["D.results"]],
    1:4,
    H.qval, # nolint
    log2FC, # nolint
    everything() # nolint
  )
  return(dgea_sum)
}
