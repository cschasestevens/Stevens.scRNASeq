#' scRNA-Seq DGEA
#'
#' Performs DGEA per cell type of a Seurat Object.
#'
#' @param so An object of class Seurat. Must contain the columns 'CellType' and 'CellGroup' in the metadata slot.
#' @param asy Character string providing the name of the assay to use for differential analysis.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.
#' @param ct Character string vector containing the name(s) of up to 2 cell type groups present in the CellGroup column.
#' The provided Seurat object is then subsetted to include only clusters from the provided group names.
#' @param MAST.comp Character string indicating the name of the MAST group comparison for conducting DGEA. MAST names are comprised of the chosen
#' variable name and the leading factor level within that variable.
#' @param MAST.name User-defined name of a DGEA comparison, given as a character string.
#' @param form1 Formula to use for MAST generalized linear model.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A list of DGEA results per cell type for the chosen group comparison, including genes missing fold changes and cell type DGEA results
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
sc.diff <- function(
    so,
    asy,
    md.list,
    ct,
    MAST.comp,
    MAST.name,
    form1,
    parl,
    core.perc
    ) {
  # Load an existing DGEA results object and skip DGEA if present
  if(file.exists("analysis/object.diff.result.rds")){
    print("A Differential analysis results object already exists for this data set! Loading existing .rds object...")
    dgea.sum <- readRDS("analysis/object.diff.result.rds")
    }

  if(!file.exists("analysis/object.diff.result.rds")){
  # Seurat object
  d <- so
  # Metadata columns
  lc <- md.list
  # Cell type column
  c <- ct
  # MAST Comparison (uses a combination of column name and leading factor level for name)
  mc <- MAST.comp
  # MAST Comparison name
  mn <- MAST.name

  ## Input
  deg.mat <- as.matrix(
    SeuratObject::GetAssayData(d,'data',assay = asy))
  deg.cols <- data.frame(d@meta.data[,c(lc)])

  ## Format input as DGEA object
  dgea.sc <- MAST::FromMatrix(
    deg.mat,
    cData = deg.cols
    )
  dgea.celltype <- unique(
    SingleCellExperiment::colData(dgea.sc)[[c]])
  list.dgea <- list("SCE" = dgea.sc,"CellType" = dgea.celltype)
  remove(d)

  if(Sys.info()[["sysname"]] != "Windows" &
     parl == TRUE){
    list.dgea.res <- setNames(parallel::mclapply(
      mc.cores = ceiling(parallel::detectCores()*core.perc),
      seq.int(1,length(list.dgea[[2]]),1),
      function(x) {
        tryCatch(
          {
          # Subset data
          s1 <- list.dgea[[1]][,SingleCellExperiment::colData(list.dgea[[1]])[[c]] == list.dgea[[2]][[x]]]
          s1.sum <- rowSums(SummarizedExperiment::assay(s1)>0)
          s1 <- s1[s1.sum/ncol(s1) >= 0.05,]
          ### create glm (generalized linear model for each variable)
          s1.fit <- MAST::zlm(
            formula = form1,
            s1,
            method = "glm",
            ebayes = F,
            parallel = F
            )
          ### Output DFs
          d.mast.sum.fun <- function(
            comp1,
            ct2,
            comp1.name
            ) {
            s1.res <- MAST::summary(
              s1.fit,
              doLRT = comp1,
              logFC = T,
              parallel = F)
            ### make dfs to display summary results by comp
            s1.dt <- reshape2::melt(
              dplyr::select(
                dplyr::filter(
                  s1.res$datatable,
                  contrast == comp1 &
                    component != "S"
                  ),
                -c("contrast")
                ),
              id.vars = c("primerid","component")
              )
            s1.dt[["vars"]] <- paste(
              s1.dt$component,
              s1.dt$variable,
              sep = "."
              )
            s1.dt <- dplyr::select(
              dplyr::mutate(
                reshape2::dcast(
                  dplyr::select(
                    dplyr::filter(
                      s1.dt,
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
                "Comparison" = comp1.name
                ),
              c(
                "CellType","Comparison","primerid",
                "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                "D.Pr(>Chisq)",everything()
                )
              )
            names(s1.dt) <- c(
              "CellType","Comparison","GENE",
              "logFC","H.pval","C.pval",
              "D.pval",
              names(
                s1.dt[8:ncol(
                  s1.dt
                )]
              )
            )
            return(s1.dt)
            }

          d1 <- d.mast.sum.fun(
            mc,
            list.dgea[[2]][[x]],
            mn
            )

            return(d1)
            },
        error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
        )
      }
    ),as.character(list.dgea[[2]]))
  }

  if(Sys.info()[["sysname"]] != "Windows" &
     parl == FALSE){
    list.dgea.res <- setNames(lapply(
      seq.int(1,length(list.dgea[[2]]),1),
      function(x) {
        tryCatch(
          {
            # Subset data
            s1 <- list.dgea[[1]][,SingleCellExperiment::colData(list.dgea[[1]])[[c]] == list.dgea[[2]][[x]]]
            s1.sum <- rowSums(SummarizedExperiment::assay(s1)>0)
            s1 <- s1[s1.sum/ncol(s1) >= 0.05,]
            ### create glm (generalized linear model for each variable)
            s1.fit <- MAST::zlm(
              formula = form1,
              s1,
              method = "glm",
              ebayes = F,
              parallel = F
              )
            ### Output DFs
            d.mast.sum.fun <- function(
            comp1,
            ct2,
            comp1.name
            ) {
              s1.res <- MAST::summary(
                s1.fit,
                doLRT = comp1,
                logFC = T,
                parallel = F)
              ### make dfs to display summary results by comp
              s1.dt <- reshape2::melt(
                dplyr::select(
                  dplyr::filter(
                    s1.res$datatable,
                    contrast == comp1 &
                      component != "S"
                  ),
                  -c("contrast")
                ),
                id.vars = c("primerid","component")
              )
              s1.dt[["vars"]] <- paste(
                s1.dt$component,
                s1.dt$variable,
                sep = "."
              )
              s1.dt <- dplyr::select(
                dplyr::mutate(
                  reshape2::dcast(
                    dplyr::select(
                      dplyr::filter(
                        s1.dt,
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
                  "Comparison" = comp1.name
                ),
                c(
                  "CellType","Comparison","primerid",
                  "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                  "D.Pr(>Chisq)",everything()
                )
              )
              names(s1.dt) <- c(
                "CellType","Comparison","GENE",
                "logFC","H.pval","C.pval",
                "D.pval",
                names(
                  s1.dt[8:ncol(
                    s1.dt
                  )]
                )
              )
              return(s1.dt)
            }

            d1 <- d.mast.sum.fun(
              mc,
              list.dgea[[2]][[x]],
              mn
            )

            return(d1)
          },
    error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
        )
      }
    ),as.character(list.dgea[[2]]))
  }

  if(Sys.info()[["sysname"]] == "Windows"){
    list.dgea.res <- setNames(
      lapply(
        seq.int(1,length(list.dgea[[2]]),1),
        function(x) {
          tryCatch(
            {
              # Subset data
              s1 <- list.dgea[[1]][,SingleCellExperiment::colData(list.dgea[[1]])[[c]] == list.dgea[[2]][[x]]]
              s1.sum <- rowSums(SummarizedExperiment::assay(s1)>0)
              s1 <- s1[s1.sum/ncol(s1) >= 0.05,]
              ### create glm (generalized linear model for each variable)
              s1.fit <- MAST::zlm(
                formula = form1,
                s1,
                method = "glm",
                ebayes = F,
                parallel = F
                )
              ### Output DFs
              d.mast.sum.fun <- function(
              comp1,
              ct2,
              comp1.name
              ) {
                s1.res <- MAST::summary(
                  s1.fit,
                  doLRT = comp1,
                  logFC = T,
                  parallel = F)
                ### make dfs to display summary results by comp
                s1.dt <- reshape2::melt(
                  dplyr::select(
                    dplyr::filter(
                      s1.res$datatable,
                      contrast == comp1 &
                        component != "S"
                    ),
                    -c("contrast")
                  ),
                  id.vars = c("primerid","component")
                )
                s1.dt[["vars"]] <- paste(
                  s1.dt$component,
                  s1.dt$variable,
                  sep = "."
                )
                s1.dt <- dplyr::select(
                  dplyr::mutate(
                    reshape2::dcast(
                      dplyr::select(
                        dplyr::filter(
                          s1.dt,
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
                    "Comparison" = comp1.name
                  ),
                  c(
                    "CellType","Comparison","primerid",
                    "logFC.coef","H.Pr(>Chisq)","C.Pr(>Chisq)",
                    "D.Pr(>Chisq)",everything()
                  )
                )
                names(s1.dt) <- c(
                  "CellType","Comparison","GENE",
                  "logFC","H.pval","C.pval",
                  "D.pval",
                  names(
                    s1.dt[8:ncol(
                      s1.dt
                    )]
                  )
                )
                return(s1.dt)
              }

              d1 <- d.mast.sum.fun(
                mc,
                list.dgea[[2]][[x]],
                mn
              )

              return(d1)
            },
    error = function(e) {print("Differential analysis unsuccessful for selected cell type...")}
        )
      }
    ),as.character(list.dgea[[2]]))
  }

  # Combine results
  dgea.comb <- list.dgea.res[lengths(list.dgea.res) > 1]
  ## isolate DGEA results with errors
  dgea.error <- list.dgea.res[lengths(list.dgea.res) <= 1]
  ## return genes for each result with missing logFC
  dgea.miss <- lapply(
    dgea.comb,
    function(x) x[is.na(x[["logFC"]]),]
    )
  dgea.res <- lapply(
    dgea.comb,
    function(x) x[!is.na(x[["logFC"]]),]
    )
  dgea.sum <- list(
    "D.results" = dplyr::bind_rows(dgea.res),
    "D.missing" = dplyr::bind_rows(dgea.miss),
    "D.errors" = dplyr::bind_rows(dgea.error)
    )

  dgea.sum[["D.results"]][["H.qval"]] <- p.adjust(
    dgea.sum[["D.results"]][["H.pval"]],
    method = "BH"
    )

  dgea.sum[["D.results"]][["log2FC"]] <- log2(
    exp(dgea.sum[["D.results"]][["logFC"]])
    )

  dgea.sum[["D.results"]] <- dplyr::select(
    dgea.sum[["D.results"]],
    1:4,
    H.qval,
    log2FC,
    everything()
    )

  ## Export list elements and save as RDS
  write.table(
    dgea.sum[[1]],
    file = "analysis/table.diff.results.txt",
    sep = "\t",
    col.names = T,
    row.names = F
    )
  write.table(
    dgea.sum[[2]],
    file = "analysis/table.diff.missFC.txt",
    sep = "\t",
    col.names = T,
    row.names = F
    )
  write.table(
    dgea.sum[[3]],
    file = "analysis/table.diff.errors.txt",
    sep = "\t",
    col.names = T,
    row.names = F
    )
  saveRDS(dgea.sum,"analysis/object.diff.result.rds")
  }

  return(dgea.sum)
  }


