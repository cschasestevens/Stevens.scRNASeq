#' Phenocycler QC
#'
#' Runs quality control of segmented phenocycler data for all samples
#'
#' @param qc_type Specific QC steps to perform (
#' select from "all" or "sum_only")
#' @param ld A list of phenocycler segmentation results containing summary data
#' and segmented cells.
#' @param col_nums1 Vector of column numbers to include from segmentation
#' summary data.
#' @param md_var Vector for renaming summary data column names.
#' @param md_num Number of metadata columns in exported cell segmentation
#' results. Only used if qc_type is "all".
#' @param samp_id Sample ID column name.
#' @param g_col (optional) Add a column providing grouping information.
#' @param ch_list List of channels exported with segmentation results. Note
#' that the channel list must match the order of channels included within
#' the results.
#' @param ch_nuc Nuclei channel name.
#' @param col_nums2 If qc_type is "all", indicate the column number to
#' use as the sample ID as well as the X and Y coordinates
#' for formatting the segmentation results table.
#' @param loess_norm Normalize signal intensities by performing LOESS
#' prior to z-score scaling.
#' @return 
#' @examples
#'
#' # pheno_data <- pc_qc(
#' #   ld = ld1,
#' #   col_nums1 = c(4, 10, 12, 13),
#' #   col_names1 = c("Code", "Area.um2", "Cells.raw", "Description"),
#' #   g_col = ifelse(grepl("KK", qc1[["Code"]]), "CF", "Norm"),
#' #   ch_list = s3,
#' #   col_nums2 = c(6, 8, 9)
#' # )
#'
#' @export
pc_qc <- function(
  qc_type = "all",
  ld,
  col_nums1,
  md_var,
  md_num = 9,
  samp_id = "Code",
  g_col = NULL,
  ch_list = NULL,
  col_nums2 = NULL,
  ch_nuc = "DAPI",
  loess_norm = FALSE
) {
  # Load data
  d1 <- ld
  cn1 <- col_nums1
  cn2 <- md_var
  sid <- samp_id
  # Perform QC
  if(qc_type == "all") { # nolint
    # Summary statistics
    ## Count section dimensions and detected cells per section
    qc1 <- setNames(dplyr::bind_rows(lapply(
      seq.int(1, length(d1), 1),
      function(x) {
        data.frame(
          "Slide" = x,
          dplyr::select(
            d1[[x]][[1]],
            sid, cn1
          )
        )
      }
    )), c("Slide", sid, cn2))
    if(is.null(g_col) == FALSE) { # nolint
      qc1 <- dplyr::select(
        dplyr::mutate(
          qc1,
          "Group" = g_col
        ),
        c("Slide", sid, "Group"), everything() # nolint
      )
    }
    # Intensity-based filtering (by nuclei and total signal)
    ## Format segmentation results (Cellpose results contain summary stats!)
    ### Formatting for exported Cellpose measurements
    qc2 <- dplyr::bind_rows(lapply(
      seq.int(1, length(d1), 1),
      function(x) {
        setNames(
          data.frame(
            "Slide" = rep(x, nrow(d1[[x]][[2]])),
            d1[[x]][[2]]
          ),
          c(
            "Slide",
            names(d1[[x]][[2]][1:md_num]),
            unlist(
              lapply(
                seq.int(1, nrow(ch_list), 1),
                function(x) {
                  st1 <- c("mean", "median", "min", "max", "sd")
                  st2 <- unlist(
                    lapply(
                      seq.int(1, length(st1), 1),
                      function(y) {
                        paste(
                          ch_list[x, 1],
                          st1[[y]],
                          sep = "_"
                        )
                      }
                    )
                  )
                  return(st2)
                }
              )
            )
          )
        )
      }
    ))
    qc2 <- dplyr::select(
      dplyr::mutate(
        qc2, "ID" = seq.int(1, nrow(qc2), 1)
      ), "ID", everything() # nolint
    )
    ## remove all columns except mean intensity per cell
    qc2 <- dplyr::select(
      qc2,
      names(qc2[, 1:(md_num + 2)]),
      names(qc2[, grepl("mean", names(qc2))])
    )
    qc2 <- setNames(
      reshape2::melt(
        qc2,
        id.vars = 1:(md_num + 2),
        variable.name = "column",
        value.name = "value"
      )[, c(
        1, 2,
        (col_nums2[[1]] + 2),
        (col_nums2[[2]] + 2),
        (col_nums2[[3]] + 2),
        (md_num + 3),
        (md_num + 4)
      )],
      c("ID", "Slide", sid, "X", "Y", "column", "value")
    )
    if(is.null(g_col) == FALSE) { # nolint
      qc2 <- dplyr::select(
        dplyr::mutate(
          qc2,
          "Group" = g_col,
          "Channel" = gsub("\\_.*", "", qc2[["column"]])
        ),
        c("ID", "Slide", sid, "Group", "X", "Y", "Channel", "value")
      )
    }
    if(is.null(g_col) == TRUE) { # nolint
      qc2 <- dplyr::select(
        dplyr::mutate(
          qc2,
          "Channel" = gsub("\\_.*", "", qc2[["column"]])
        ),
        c("ID", "Slide", sid, "X", "Y", "Channel", "value")
      )
    }
    # Calculate sum intensity per cell
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
      qc_sum <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = 2,
          seq.int(1, length(unique(qc2[[sid]])), 1),
          function(x) {
            md1 <- unique(qc2[[sid]])
            d1 <- qc2[qc2[[sid]] == md1[[x]], ]
            d1 <- setNames(
              aggregate(
                d1[["value"]],
                list(d1[["ID"]]),
                function(y) sum(y)
              ),
              c("ID", "sum.int.raw")
            )
            d1 <- dplyr::select(
              dplyr::mutate(
                d1,
                "Code" = rep(md1[[x]], nrow(d1))
              ),
              "Code",
              everything() # nolint
            )
            return(d1)
          }
        )
      )
    }
    if(Sys.info()[["sysname"]] == "Windows") { # nolint
      qc_sum <- dplyr::bind_rows(
        lapply(
          seq.int(1, length(unique(qc2[[sid]])), 1),
          function(x) {
            md1 <- unique(qc2[[sid]])
            d1 <- qc2[qc2[[sid]] == md1[[x]], ]
            d1 <- setNames(
              aggregate(
                d1[["value"]],
                list(d1[["ID"]]),
                function(y) sum(y)
              ),
              c("ID", "sum.int.raw")
            )
            d1 <- dplyr::select(
              dplyr::mutate(
                d1,
                "Code" = rep(md1[[x]], nrow(d1))
              ),
              "Code",
              everything() # nolint
            )
            return(d1)
          }
        )
      )
    }
    ### Flag cells based on Total intensity
    filt1 <- qc_sum[
      qc_sum[["sum.int.raw"]] < quantile(qc_sum[["sum.int.raw"]], 0.01) |
        qc_sum[["sum.int.raw"]] > quantile(qc_sum[["sum.int.raw"]], 0.99),
      "ID"
    ]
    ### Flag cells based on DAPI intensity
    filt2 <- qc2[qc2[["Channel"]] == ch_nuc, ]
    filt2 <- filt2[
      filt2[["value"]] < quantile(filt2[["value"]], 0.01) |
        filt2[["value"]] > quantile(filt2[["value"]], 0.99),
      "ID"
    ]
    ### Remove flagged cells from individual and aggregated data
    qc2_filt <- qc2[
      qc2[["ID"]] %in% filt1 == FALSE &
        qc2[["ID"]] %in% filt2 == FALSE,
    ]
    qc_sum_filt <- qc_sum[
      qc_sum[["ID"]] %in% filt1 == FALSE &
        qc_sum[["ID"]] %in% filt2 == FALSE,
    ]
    fprp <- round((1 - (nrow(qc2_filt) / nrow(qc2))), digits = 2) * 100
    print(paste(fprp, "% ", "of cells removed...", sep = ""))
    qc_sum_cnt <- setNames(dplyr::count(
      qc_sum_filt,
      .data[["Code"]] # nolint
    ), c(sid, "Cells.filt"))
    qc1 <- dplyr::left_join(
      qc1,
      qc_sum_cnt,
      by = sid
    )
    ### Merge coordinates with sum intensity data frame
    qc_sum_filt <- dplyr::left_join(
      qc_sum_filt,
      qc2_filt[!duplicated(qc2_filt[["ID"]]), c("ID", "X", "Y")],
      by = "ID"
    )
    # Normalizations
    ## Z-score only
    #### formula (x - mean)/sd # nolint
    qc2_filt[["value.z"]] <- (qc2_filt[["value"]] - mean(qc2_filt[["value"]])) /
      sd(qc2_filt[["value"]])
    ## LOESS + Z-score (Experimental) # nolint
    if(loess_norm == TRUE) { # nolint
      
    }
    #### START HERE 5/15/25 LOESS norm, plot all slides, antibodies, and format output for Python scVI dim reduction and UMAP

  }



  
  return(d)
}
