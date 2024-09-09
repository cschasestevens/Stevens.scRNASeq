#' Define Processing Parameters
#'
#' Creates a data frame of processing parameters
#' used in scRNA-Seq data processing by Seurat.
#'
#' @param d.path Character string indicating the path
#' to a set of CellRanger files for data processing.
#' @param study.md A data frame containing sample
#' metadata variables specific to an individual study.
#' @return Data frame containing a list of parameters
#' to use for scRNA-Seq data processing by Seurat.
#' @examples
#'
#' # # Dataset input parameters
#' # list_params <- create_proc_params(
#' #   "data/",
#' #   data.frame(
#' #     # ID column name
#' #     # (splits by underscore, should be listed
#' #     # first in folder name as in 's01_1_KO_a')
#' #     "Code" = unlist(
#' #        lapply(
#' #          strsplit(
#' #            basename(list.files("data/")),
#' #            "_",
#' #            fixed = T
#' #          ),
#' #          "[",
#' #          1
#' #        )
#' #      ),
#' #     # 1st metadata column (include in CellRanger folder name)
#' #     "Knockout" = as.factor(
#' #                    ifelse(
#' #                      grepl("NG",basename(list.files("data/"))),
#' #                      "ctrl",
#' #                      "KO"
#' #                    )
#' #                  ),
#' #     # 2nd metadata column
#' #     "Region" = as.factor(
#' #                  ifelse(
#' #                    grepl("LAE",basename(list.files("data/"))),
#' #                    "1",
#' #                    "2"
#' #                  )
#' #                ),
#' #     # 3rd metadata column (add/remove columns as needed)
#' #     "Time" = as.factor(
#' #                ifelse(
#' #                  grepl("D28",basename(list.files("data/"))),
#' #                  "a",
#' #                  "b"
#' #                )
#' #              )
#' #   )
#' # )
#'
#' @export
create_proc_params <- function(d_path, study_md) {

  list_params <- data.frame(
    # Universal columns
    data.frame(
      # Sample Number
      Sample.No = seq(1:length( # nolint
        basename(
          list.files(d_path)
        )
      )
      ),
      # Individual file names (uses CellRanger folder name by default)
      File.ID = basename(list.files(d_path)),
      # Data file paths (Location of CellRanger files: 'Data/' by default)
      Path = paste(d_path, list.files(d_path), sep = ""),
      # Path to feature files
      Path.feat = paste(
        d_path,
        list.files(d_path),
        "/filtered_feature_bc_matrix/features.tsv.gz",
        sep = ""
      )
    ),

    # Dataset-specific columns
    study_md
  )

  return(list_params)
}