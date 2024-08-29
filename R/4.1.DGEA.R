#' scRNA-Seq DGEA
#'
#' Performs DGEA per cell type of a Seurat Object.
#'
#' @param so An object of class Seurat. Must contain the columns
#' 'CellType' and 'CellGroup' in the metadata slot.
#' @param asy Character string providing the name of the assay
#' to use for differential analysis.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param ct Character string vector containing the name(s) of up to
#' 2 cell type groups present in the CellGroup column.
#' The provided Seurat object is then subsetted to include
#' only clusters from the provided group names.
#' @param MAST_comp Character string indicating the name of the MAST
#' group comparison for conducting DGEA. MAST names are comprised of the chosen
#' variable name and the leading factor level within that variable.
#' @param MAST_name User-defined name of a DGEA comparison,
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
#' # dgea.output <- sc.diff(
#' # # Seurat object
#' # d1,
#' # # Assay
#' # "RNA",
#' # # metadata column list
#' # c(list.p.cols,"CellType","nFeature_RNA"),
#' # # Cell type column name
#' # "CellType",
#' # # MAST comparison name
#' # "KnockoutKO",
#' # # MAST name (user-provided)
#' # "KO vs. NG ",
#' # # Formula
#' # as.formula(
#' #   paste(
#' #     "~","Knockout","+",
#' #     "nFeature_RNA",
#' #     sep = " "
#' #     )
#' #   ),
#' # # run in parallel? (Set to FALSE if on Windows)
#' # TRUE,
#' # # core percentage to use
#' # 0.5
#' # )
#'
#' @export
sc_diff <- function(
  so,
  asy,
  md_list,
  ct,
  mast_comp,
  mast_name,
  form1,
  parl,
  core_perc,
  atac_type = NULL
) {
  # Load an existing DGEA results object and skip DGEA if present
  if(file.exists("analysis/object.diff.result.rds")) { # nolint
    print("A Differential analysis results object already exists
    for this data set! Loading existing .rds object...")
    dgea_sum <- readRDS("analysis/object.diff.result.rds")
  }

  if(!file.exists("analysis/object.diff.result.rds")) { # nolint
    if(is.null(atac_type)) { # nolint
      # Seurat object
      d <- so
      # Metadata columns
      lc <- md_list
      # Cell type column
      c <- ct
    }

    if(!is.null(atac_type)) { # nolint
      library(Seurat)
      # Seurat object
      d <- so
      # Metadata columns
      lc <- md_list
      # Cell type column
      c <- ct

      ## Seurat
      d <- subset(d, subset = CellType == ct) # nolint
    }

    # MAST Comparison (combines column name and leading factor level for name)
    mc <- mast_comp
    # MAST Comparison name
    mn <- mast_name

    ## Input
    deg_mat <- as.matrix(SeuratObject::GetAssayData(d, "data", assay = asy))
    deg_cols <- data.frame(d@meta.data[, c(lc)])

    ## Format input as DGEA object
    dgea_sc <- MAST::FromMatrix(
      deg_mat,
      cData = deg_cols
    )
    dgea_celltype <- unique(SingleCellExperiment::colData(dgea_sc)[[c]])
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
                      contrast == comp1 &
                        component != "S"
                      ),
                    -c("contrast")
                    ),
                  id.vars = c("primerid","component")
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
                          vars != "logFC.Pr(>Chisq)" &
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
                    "CellType","Comparison","primerid",
                    "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                    "D.Pr(>Chisq)",everything()
                    )
                  )
                names(s1_dt) <- c(
                  "CellType","Comparison","GENE",
                  "logFC","H.pval","C.pval",
                  "D.pval",
                  names(
                    s1_dt[8:ncol(
                      s1_dt
                    )]
                  )
                )
                return(s1_dt)
                }

            d1 <- d_mast_sum_fun(
              mc,
              list_dgea[[2]][[x]],
              mn
              )

              return(d1)
              },
          error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
          )
        }
      ),as.character(list_dgea[[2]]))
    }

    if(Sys.info()[["sysname"]] != "Windows" &
       parl == FALSE){
      list_dgea_res <- setNames(lapply(
        seq.int(1,length(list_dgea[[2]]),1),
        function(x) {
          tryCatch(
            {
              # Subset data
              s1 <- list_dgea[[1]][,SingleCellExperiment::colData(list_dgea[[1]])[[c]] == list_dgea[[2]][[x]]]
              s1_sum <- rowSums(SummarizedExperiment::assay(s1)>0)
              s1 <- s1[s1_sum/ncol(s1) >= 0.05,]
              ### create glm (generalized linear model for each variable)
              s1_fit <- MAST::zlm(
                formula = form1,
                s1,
                method = "glm",
                ebayes = F,
                parallel = F
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
                  logFC = T,
                  parallel = F)
                ### make dfs to display summary results by comp
                s1_dt <- reshape2::melt(
                  dplyr::select(
                    dplyr::filter(
                      s1_res$datatable,
                      contrast == comp1 &
                        component != "S"
                    ),
                    -c("contrast")
                  ),
                  id.vars = c("primerid","component")
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
                          vars != "logFC.Pr(>Chisq)" &
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
                    "CellType","Comparison","primerid",
                    "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                    "D.Pr(>Chisq)",everything()
                  )
                )
                names(s1_dt) <- c(
                  "CellType","Comparison","GENE",
                  "logFC","H.pval","C.pval",
                  "D.pval",
                  names(
                    s1_dt[8:ncol(
                      s1_dt
                    )]
                  )
                )
                return(s1_dt)
              }

              d1 <- d_mast_sum_fun(
                mc,
                list_dgea[[2]][[x]],
                mn
              )

              return(d1)
            },
      error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
          )
        }
      ),as.character(list_dgea[[2]]))
    }

    if(Sys.info()[["sysname"]] == "Windows"){
      list_dgea_res <- setNames(
        lapply(
          seq.int(1,length(list_dgea[[2]]),1),
          function(x) {
            tryCatch(
              {
                # Subset data
                s1 <- list_dgea[[1]][,SingleCellExperiment::colData(list_dgea[[1]])[[c]] == list_dgea[[2]][[x]]]
                s1_sum <- rowSums(SummarizedExperiment::assay(s1)>0)
                s1 <- s1[s1_sum/ncol(s1) >= 0.05,]
                ### create glm (generalized linear model for each variable)
                s1_fit <- MAST::zlm(
                  formula = form1,
                  s1,
                  method = "glm",
                  ebayes = F,
                  parallel = F
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
                    logFC = T,
                    parallel = F)
                  ### make dfs to display summary results by comp
                  s1_dt <- reshape2::melt(
                    dplyr::select(
                      dplyr::filter(
                        s1_res$datatable,
                        contrast == comp1 &
                          component != "S"
                      ),
                      -c("contrast")
                    ),
                    id.vars = c("primerid","component")
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
                            vars != "logFC.Pr(>Chisq)" &
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
                      "CellType","Comparison","primerid",
                      "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                      "D.Pr(>Chisq)",everything()
                    )
                  )
                  names(s1_dt) <- c(
                    "CellType","Comparison","GENE",
                    "logFC","H.pval","C.pval",
                    "D.pval",
                    names(
                      s1_dt[8:ncol(
                        s1_dt
                      )]
                    )
                  )
                  return(s1_dt)
                }

                d1 <- d_mast_sum_fun(
                  mc,
                  list_dgea[[2]][[x]],
                  mn
                )

                return(d1)
              },
      error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
          )
        }
      ),as.character(list_dgea[[2]]))
    }

    # Combine results
    dgea_comb <- list_dgea_res[lengths(list_dgea_res) > 1]
    ## isolate DGEA results with errors
    dgea_error <- list_dgea_res[lengths(list_dgea_res) <= 1]
    ## return genes for each result with missing logFC
    dgea_miss <- lapply(
      dgea_comb,
      function(x) x[is.na(x[["logFC"]]),]
      )
    dgea_res <- lapply(
      dgea_comb,
      function(x) x[!is.na(x[["logFC"]]),]
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
      H.qval,
      log2FC,
      everything()
      )

    ## Export list elements and save as RDS
    write.table(
      dgea_sum[[1]],
      file = "analysis/table.diff.results.txt",
      sep = "\t",
      col.names = T,
      row.names = F
      )
    write.table(
      dgea_sum[[2]],
      file = "analysis/table.diff.missFC.txt",
      sep = "\t",
      col.names = T,
      row.names = F
      )
    write.table(
      dgea_sum[[3]],
      file = "analysis/table.diff.errors.txt",
      sep = "\t",
      col.names = T,
      row.names = F
      )
    saveRDS(dgea_sum,"analysis/object.diff.result.rds")
    }

    return(dgea_sum)
    }


