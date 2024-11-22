#' CellChat Based Cell-cell Communication
#'
#' Conducts CellChat for a selected Seurat object, either provided as
#' a whole data set or split by a treatment variable for comparison.
#'
#' @param so A Seurat object. Must contain a gene expression assay.
#' @param asy1 RNA assay to use.
#' @param g_name Group name for comparing cell-cell communication between
#' treatments, provided as a character string. If only one group is present,
#' include an additional variable in the selected Seurat data set indicating
#' that all cells belong to the same group.
#' @param s_name Variable name containing individual sample names.
#' @param ct_col Cell type column name to use in the selected Seurat object.
#' @return A data frame including either pathway or individual
#' ligand-receptor interactions.
#' @examples
#'
#' # ccmerge <- sc_cc_run(
#' #   so = d,
#' #   asy1 = "RNA",
#' #   g_name = "Group",
#' #   s_name = "Code",
#' #   ct_col = "CellType"
#' # )
#'
#' @export
sc_cc_run <- function(
  so,
  asy1,
  g_name,
  s_name,
  ct_col
) {
  # Load and extract data (subset before creating CellChat if desired)
  dcc <- so
  Seurat::DefaultAssay(dcc) <- asy1
  ## Subset by condition
  dcc_l <- setNames(
    parallel::mclapply(
      mc.cores = 2,
      seq.int(1, length(unique(dcc@meta.data[[g_name]])), 1),
      function(x) {
        dcc <- dcc[, dcc[[g_name]] == unique(dcc@meta.data[[g_name]])[[x]]]
        dcc@meta.data[["CellType"]] <- factor(
          as.character(dcc@meta.data[["CellType"]]),
          levels = gtools::mixedsort(
            unique(as.character(dcc@meta.data[["CellType"]]))
          )
        )
        dcc_d <- dcc[[asy1]]$data
        dcc_m <- dcc@meta.data
        dcc_m[["samples"]] <- as.factor(dcc_m[[s_name]])
        ct_col1 <- ct_col
        gc(reset = TRUE)
        # Create CellChat object
        cc1 <- CellChat::createCellChat(
          object = dcc_d,
          meta = dcc_m,
          group.by = ct_col1
        )
        # Set database
        CellChatDB <- CellChatDB.human # nolint
        CellChatDB.use <- CellChatDB # nolint
        cc1@DB <- CellChatDB.use
        # Pre-process to determine over-expressed ligands and receptors
        cc1 <- CellChat::subsetData(cc1)
        future::plan("multisession", workers = parallel::detectCores() * 0.25)
        cc1 <- CellChat::identifyOverExpressedGenes(cc1)
        cc1 <- CellChat::identifyOverExpressedInteractions(cc1)
        cc1 <- CellChat::smoothData(cc1, adj = PPI.human) # nolint
        # Compute communication probability
        options(future.globals.maxSize = 3000 * 1024^2)
        cc1 <- CellChat::computeCommunProb(cc1, type = "triMean")
        cc1 <- CellChat::computeCommunProbPathway(cc1)
        # Calculate aggregated cc communication network
        cc1 <- CellChat::aggregateNet(cc1)
        options(future.globals.maxSize = 500 * 1024^2)
        return(cc1)
      }
    ),
    unique(dcc@meta.data[[g_name]])
  )
  saveRDS(dcc_l, "analysis/data.cellchat.list.rds")
  return(dcc_l)
}

#' CellChat Chord Diagram
#'
#' Creates a chord diagram from a provided CellChat results data frame.
#'
#' @param ccdf A CellChat data frame generat.
#' @param title1 Plot title, provided as a character string.
#' @param cc_type Type of plot to use (either "total" of "comp"). Selecting
#' "total" plots the total number of interactions between cell types, whereas
#' "comp" compares the fold difference between two treatment groups
#' @param pw_sel Character string of a pathway to use for filtering the
#' provided CellChat object.
#' @param g_name Group name to split data by if cc_type is "comp."
#' @param pw Plot width (in cm).
#' @param ph Plot height (in cm).
#' @param fs1 Font size (between 0.1 and 1).
#' @param fd1 Gap between cell type labels and plot (generally between 1 and 3).
#' @return A chord diagram displaying either pathway or individual
#' ligand-receptor interactions.
#' @examples
#'
#' # pchrd <- sc_cc_chrd(
#' #   ccdf = ccm[[1]],
#' #   title1 = "CellChat: Total L-R Interactions",
#' #   cc_type = "total",
#' #   pw = 12,
#' #   ph = 12,
#' #   fs1 = 0.6,
#' #   fd1 = 2.8
#' # )
#'
#' @export
sc_cc_chrd <- function(
  ccdf,
  title1,
  title2,
  cc_type,
  g_name = NULL,
  pw_sel = NULL,
  pw,
  ph,
  fs1,
  fd1,
  lgx = 0.1,
  lgy = 0.1,
  g_ind,
  list_ct
) {
  # Input CC data frame
  chrd1 <- ccdf
  # Set parameters
  ## circlize params
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.degree = 5,
    track.margin = c(-0.1, 0.1),
    points.overflow.warning = FALSE
  )
  ## base params
  par(mar = rep(0.5, 4))

  # Plots
  ## All communications
  if(is.null(pw_sel)) { # nolint
    ## Total interactions
    if(cc_type == "total") { # nolint
      col1 <- setNames(
        col_univ()[1:length(unique( # nolint
        c(
          as.character(ccm[[1]][["target"]]),
          as.character(ccm[[1]][["source"]])
        )
        ))],
        gtools::mixedsort(unique(
          c(
            as.character(ccm[[1]][["target"]]),
            as.character(ccm[[1]][["source"]])
          )
        ))
      )
col1
      chrd1 <- chrd1[
        chrd1[["Group"]] == "NonCF" &
          chrd1[["target"]] %in% list_ct,
      ]
      col_fun <- col1[col1 %in% unique(
        c(
          as.character(chrd1[["target"]]),
          as.character(chrd1[["source"]])
        )
        )]

      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = chrd1,
        grid.col = col_fun, # nolint
        transparency = 0.25,
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      dev.off()
    }
    # Comparison between two datasets
    if(cc_type == "comp") { # nolint
      col_fun <- circlize::colorRamp2(
        c(min(chrd1[[3]]), 1, max(chrd1[[3]])),
        c("dodgerblue3", "white", "firebrick2"),
        transparency = 0.25
      )
      grid_col <- setNames(
        col_univ()[1:length(unique(as.character(chrd1[[2]])))], # nolint
        unique(as.character(chrd1[[2]]))
      )
      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = chrd1,
        grid.col = grid_col, # nolint
        col = col_fun,
        transparency = 0.25,
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      dev.off()
    }
  }
    if(!is.null(pw_sel)) { # nolint
    chrd1 <- chrd1[chrd1[["pathway_name"]] == pw_sel, ]
    ## Total interactions
    if(cc_type == "total") { # nolint
      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = chrd1,
        grid.col = col_univ()[1:length(unique(chrd1[["target"]]))], # nolint
        transparency = 0.25,
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      dev.off()
    }
    # Comparison between two datasets
    if(cc_type == "comp") { # nolint
      # calculate differences in interactions between groups
      ## add fold change column
      cc_comp <- chrd1[chrd1[["pathway_name"]] == pw_sel, ]
      cc_comp2 <- dplyr::count(
        cc_comp,
        .data[[g_name]], # nolint
        source, # nolint
        target # nolint
      )
      cc_comp3 <- reshape2::dcast(
        cc_comp2,
        source + target ~ cc_comp2[[g_name]], value.var = "n"
      )
      cc_comp3[is.na(cc_comp3)] <- 0
      cc_comp3[["ratio"]] <- cc_comp3[[3]] / (cc_comp3[[3]] + cc_comp3[[4]])
      cc_comp <- dplyr::left_join(
        cc_comp,
        cc_comp3[, c("source", "target", "ratio")],
        by = c("source", "target")
      )
      write.table(
        cc_comp,
        "analysis/table.cellchat.LR.IL1.compare.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
      cc_plot <- cc_comp[, c("source", "target", "Group")]
      # color schemes
      col_fun <- circlize::colorRamp2(
        c(min(cc_comp[["ratio"]]), 0.5, max(cc_comp[["ratio"]])),
        c("dodgerblue3", "white", "firebrick2"),
        transparency = 0.25
      )
      grid_col <- setNames(
        col_univ()[ # nolint
          1:( # nolint
            length(unique(c(
              as.character(cc_comp[["target"]]),
              as.character(cc_comp[["source"]])
            )))
          )
        ],
        gtools::mixedsort(
          unique(c(
            as.character(cc_comp[["target"]]),
            as.character(cc_comp[["source"]])
          ))
        )
      )
      png(
        paste(
          "analysis/",
          title1,
          ".png",
          sep = ""
        ),
        width = pw,
        height = ph,
        res = 1200,
        units = "cm"
      )
      circlize::chordDiagram(
        x = cc_plot,
        grid.col = grid_col, # nolint
        col = ifelse(
          cc_plot[[g_name]] == unique(cc_plot[[g_name]])[[1]],
          col_univ()[[1]], # nolint
          col_univ()[[2]]
        ),
        directional = 1,
        direction.type = c("arrows", "diffHeight"),
        diffHeight  = -0.04,
        annotationTrack = "grid",
        annotationTrackHeight = c(0.05, 0.1),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.largest.ontop = TRUE
      )
      circlize::circos.trackPlotRegion(
        track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          sector_index <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(
            x = mean(xlim),
            y = fd1,
            labels = sector_index,
            facing = "bending",
            cex = fs1
          )
        }
      )
      lgd <- ComplexHeatmap::Legend(
        at = unique(cc_plot[[g_name]]),
        type = "grid",
        legend_gp = grid::gpar(fill = col_univ()), # nolint
        title = "Group"
      )
      ComplexHeatmap::draw(
        lgd,
        x = grid::unit(1, "npc") - grid::unit(lgx, "mm"), # nolint
        y = grid::unit(lgy, "mm"),
        just = c("right", "bottom")
      )
      circlize::circos.clear()
      text(-0, 1.02, title2, cex = 1)
      dev.off()
    }
  }
}
