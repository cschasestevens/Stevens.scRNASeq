#' scATAC-Seq Coverage Plot
#'
#' Generates a coverage plot from a Signac ChromatinAssay.
#' Requires scATAC-Seq peak information and a reference GRanges object
#' for plotting gene positions.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @param g_name Gene to plot, provided as a character string.
#' @param bp_window Numeric value indicating the number of base pairs to
#' extend the plotting window on each end of the selected gene's location.
#' Useful for visualizing peaks corresponding to neighboring genes.
#' @param md_list Character vector of up to 2 metadata variables for
#' stratifying peak data.
#' @return A coverage plot including the specified gene track and all other
#' genes present within the specified window.
#' @examples
#'
#' # p_cov <- sc_coverage_plot(
#' #   readRDS("analysis/data.annotated.withTFs.rds"),rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' # ),"TMEM45A", 20000, c("CellType", "Airway"))
#'
#' @export
sc_coverage_plot <- function(
  so,
  dref,
  g_name,
  bp_window,
  md_list
) {
  d <- so
  ref_gene <- dref

  # Format reference gene annotation file
  ref_gene <- dplyr::as_tibble(ref_gene)
  # Extract detected genes from reference gene list
  head(d@assays$RNA$counts)
  g_data <- ref_gene[
    ref_gene[["gene_name"]] %in% rownames(d@assays$RNA$counts),
  ]
  head(g_data)

  # Extract individual gene location from Seurat object and map
  # against reference
  g_loc <- g_data[
    g_data$gene_name == g_name &
      g_data$type == "gene",
    c("start", "end", "width", "seqnames")
  ]
  p <- d@assays$ATAC@meta.features
  p[["ID"]] <- seq.int(1, nrow(p), 1)
  g_range <- p[
    p$start >= g_loc$start - bp_window &
      p$end <= g_loc$end + bp_window,
  ]
  g_range <- g_range[g_range$seqnames == as.character(g_loc$seqnames), ]
  g_range[["p.ID"]] <- paste(
    g_range$ID,
    ".",
    g_range$nearestGene,
    sep = ""
  )
  g_range <- g_range[
    as.character(g_range[["seqnames"]]) == g_loc[["seqnames"]],
  ]
  nrow(g_range)
  head(g_range)

  # Extract raw peak counts from ATAC assay
  if(length(md_list) == 1) { # nolint
    p2 <- setNames(data.frame(
      "col1" = d@meta.data[[md_list[[1]]]],
      setNames(
        as.data.frame(t(as.matrix(d@assays$ATAC$counts[g_range$ID, ]))),
        c(g_range$p.ID)
      )
    ), c(md_list[[1]], g_range[["p.ID"]]))
  }
  if(length(md_list) == 2) { # nolint
    p2 <- setNames(data.frame(
      "col1" = d@meta.data[[md_list[[1]]]],
      "col2" = d@meta.data[[md_list[[2]]]],
      setNames(
        as.data.frame(t(as.matrix(d@assays$ATAC$counts[g_range$ID, ]))),
        c(g_range$p.ID)
      )
    ), c(md_list[[1]], md_list[[2]], g_range[["p.ID"]]))
  }
  head(p2)
  nrow(p2)

  ## Calculate sum of reads per fragment per cell type
  ## Normalize reads using the following: (f.raw/n.cells)*mean(n.reads)
  if(length(md_list) == 1) { # nolint
    p2 <- dplyr::group_by(
      p2,
      .data[[md_list[[1]]]] # nolint
    )
    # raw signal
    d_raw <- data.frame(
      "col1" = levels(p2[[md_list[[1]]]]),
      as.data.frame(
        lapply(
          p2[3:ncol(p2)],
          function(y) {
            aggregate(
              y,
              list(p2[[md_list[[1]]]]),
              function(x) sum(x)
            )[[2]]
          }
        )
      )
    )
    head(d_raw)

    ## group scaling factor
    sc_norm_atac <- function(
      # raw data frame
      df_raw,
      # raw data frame (sum frequency per group)
      df_sum,
      # raw signal
      f_raw,
      # number of metadata columns
      md
    ) {
      (
        # raw signal
        f_raw /
          (
            # total cells per group
            aggregate(
              df_raw[[(md + 1)]],
              list(df_raw[[md_list[[1]]]]),
              function(x) length(x)
            )[[2]]
          )
      ) *
        (
          # average sequencing depth per group
          rowMeans(
            df_sum[2:ncol(df_sum)]
          )
        )
    }

    # normalized signal
    d_norm <- setNames(
      data.frame(
        "col1" = d_raw[["col1"]],
        as.data.frame(
          lapply(
            seq.int(2, ncol(d_raw), 1),
            function(y) {
              sc_norm_atac(
                p2,
                d_raw,
                d_raw[[y]],
                2
              )
            }
          )
        )
      ),
      c(md_list[[1]], names(p2[3:ncol(p2)]))
    )
    head(d_norm)

    # format for plotting
    d_norm <- setNames(
      reshape2::melt(
        d_norm,
        id.vars = 1
      ),
      c(md_list[[1]], "p.ID", "freq.norm")
    )
    head(d_norm)

    d_norm <- dplyr::left_join(
      d_norm,
      g_range[, c("p.ID", "start", "end")],
      by = "p.ID"
    )
    head(d_norm)

    ## 1.Split by CellType;
    ## 2.Merge each with bp interval
    ## (in steps of 500 bp [this is the fragment length]);
    ## 3.Smooth each track
    ## 4.Bind rows and plot
    d_norm <- dplyr::bind_rows(
      setNames(
        lapply(
          unique(d_norm[[md_list[[1]]]]),
          function(x) {
            d <- d_norm[d_norm[[md_list[[1]]]] == x, ]
            p2_int <- data.frame(
              "start" = seq.int(
                min(d[["start"]]) - 5000,
                max(d[["start"]] + 5000),
                by = 501
              )
            )
            d <- dplyr::full_join(
              p2_int,
              d,
              by = "start"
            )
            d <- d[order(d[["start"]]), ]
            d[is.na(d[["freq.norm"]]), "freq.norm"] <- 0
            d[is.na(d[[md_list[[1]]]]), md_list[[1]]] <- x
            ## Perform loess smoothing of tracks
            d[["freq.loess"]] <- lowess(
              x = d[["start"]],
              y = d[["freq.norm"]],
              f = 0.05,
              iter = 5,
              delta = 0
            )[[2]]
            d[d[["freq.loess"]] < 0, "freq.loess"] <- 0
            return(d)
          }
        ),
        unique(d_norm[[md_list[[1]]]])
      )
    )
    head(d_norm)
    nrow(d_norm)

    d_norm[[md_list[[1]]]] <- factor(
      d_norm[[md_list[[1]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[1]]]])
        )
      )
    )
  }

  if(length(md_list) == 2) { # nolint
    p2 <- dplyr::group_by(
      p2,
      .data[[md_list[[1]]]], # nolint
      .data[[md_list[[2]]]]
    )
    # raw signal
    d_raw <- data.frame(
      "col1" = levels(p2[[md_list[[1]]]]),
      "col2" = levels(p2[[md_list[[2]]]]),
      as.data.frame(
        lapply(
          p2[3:ncol(p2)],
          function(y) {
            aggregate(
              y,
              list(
                p2[[md_list[[1]]]],
                p2[md_list[[2]]]
              ),
              function(x) sum(x)
            )[[2]]
          }
        )
      )
    )
    head(d_raw)

    ## group scaling factor
    sc_norm_atac <- function(
      # raw data frame
      df_raw,
      # raw data frame (sum frequency per group)
      df_sum,
      # raw signal
      f_raw,
      # number of metadata columns
      md
    ) {
      (
        # raw signal
        f_raw /
          (
            # total cells per group
            aggregate(
              df_raw[[(md + 1)]],
              list(
                df_raw[[md_list[[1]]]],
                df_raw[[md_list[[2]]]]
              ),
              function(x) length(x)
            )[[2]]
          )
      ) *
        (
          # average sequencing depth per group
          rowMeans(
            df_sum[3:ncol(df_sum)]
          )
        )
    }

    # normalized signal
    d_norm <- setNames(
      data.frame(
        "col1" = d_raw[["col1"]],
        "col2" = d_raw[["col2"]],
        as.data.frame(
          lapply(
            seq.int(3, ncol(d_raw), 1),
            function(y) {
              sc_norm_atac(
                p2,
                d_raw,
                d_raw[[y]],
                3
              )
            }
          )
        )
      ),
      c(md_list[[1]], md_list[[2]], names(p2[3:ncol(p2)]))
    )
    head(d_norm)

    # format for plotting
    d_norm <- setNames(
      reshape2::melt(
        d_norm,
        id.vars = 1
      ),
      c(md_list[[1]], md_list[[2]], "p.ID", "freq.norm")
    )
    head(d_norm)

    d_norm <- dplyr::left_join(
      d_norm,
      g_range[, c("p.ID", "start", "end")],
      by = "p.ID"
    )
    head(d_norm)

    ## 1.Split by CellType;
    ## 2.Merge each with bp interval
    ## (in steps of 500 bp [this is the fragment length]);
    ## 3.Smooth each track
    ## 4.Bind rows and plot
    d_norm <- dplyr::bind_rows(
      setNames(
        lapply(
          seq.int(1, unique(d_norm[, md_list]), 1),
          function(x) {
            d1 <- unique(d_norm[, md_list])
            d <- d_norm[
              d_norm[[md_list[[1]]]] == d1[x, md_list[[1]]] &
                d_norm[[md_list[[2]]]] == d1[x, md_list[[2]]],
            ]
            p2_int <- data.frame(
              "start" = seq.int(
                min(d[["start"]]) - 5000,
                max(d[["start"]] + 5000),
                by = 501
              )
            )
            d <- dplyr::full_join(
              p2_int,
              d,
              by = "start"
            )
            d <- d[order(d[["start"]]), ]
            d[is.na(d[["freq.norm"]]), "freq.norm"] <- 0
            d[is.na(d[[md_list[[1]]]]), md_list[[1]]] <- d1[x, md_list[[1]]]
            d[is.na(d[[md_list[[2]]]]), md_list[[2]]] <- d1[x, md_list[[2]]]

            ## Perform loess smoothing of tracks
            d[["freq.loess"]] <- lowess(
              x = d[["start"]],
              y = d[["freq.norm"]],
              f = 0.05,
              iter = 5,
              delta = 0
            )[[2]]
            d[d[["freq.loess"]] < 0, "freq.loess"] <- 0
            return(d)
          }
        ),
        paste(
          unique(d_norm[, md_list])[[1]],
          unique(d_norm[, md_list])[[2]],
          sep = "."
        )
      )
    )
    head(d_norm)
    nrow(d_norm)

    d_norm[[md_list[[1]]]] <- factor(
      d_norm[[md_list[[1]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[1]]]])
        )
      )
    )

    d_norm[[md_list[[2]]]] <- factor(
      d_norm[[md_list[[2]]]],
      levels = c(
        gtools::mixedsort(
          unique(d_norm[[md_list[[2]]]])
        )
      )
    )

  }

  # Plot
  p_cov <- ggplot2::ggplot(
    data = d_norm
  ) +
    ggplot2::geom_area(
      ggplot2::aes(
        x = .data[["start"]], # nolint
        y = .data[["freq.loess"]],
        fill = .data[[md_list[[1]]]]
      ),
      alpha = 0.7,
      show.legend = FALSE
    ) +
    sc_theme1() +
    ggplot2::scale_fill_manual(
      values = col_univ()
    ) +
    ggplot2::scale_x_continuous(
      name = "Position (bp)",
      limits = c(
        (min(d_norm[["start"]]) - 5000),
        max(d_norm[["start"]] + 5000)
      )
    ) +
    ggplot2::facet_wrap(
      ~ .data[[md_list[[1]]]],
      ncol = 1,
      strip.position = "top"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      strip.text.x.top = ggplot2::element_text(
        angle = 0,
        margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        size = 10
      ),
      strip.background.x = ggplot2::element_rect(fill = "grey95"),
      panel.spacing.y = grid::unit(0, "cm"),
      axis.text.y = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(color = "black"),
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.title.x.bottom = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = paste(
        "LOESS-smoothed Normalized Signal",
        paste(
          "(Range ",
          round(min(d_norm[["freq.loess"]]), digits = 2),
          " - ",
          round(max(d_norm[["freq.loess"]]), digits = 2), ")", sep = ""
        ),
        sep = " "
      )
    ) +
    ggplot2::ggtitle(
      paste(
        g_loc[["seqnames"]], ":",
        min(d_norm[["start"]] - 5000),
        "-",
        max(d_norm[["start"]] + 5000),
        " (",
        g_name,
        ")",
        sep = ""
      )
    )


  p_pos <- g_data[
    g_data$start >= g_loc$start - bp_window &
      g_data$end <= g_loc$end + bp_window,
  ]
  p_pos
  p_pos <- p_pos[p_pos[["seqnames"]] == as.character(g_loc[["seqnames"]]), ]
  head(p_pos)
  p_pos <- p_pos[p_pos[["type"]] == "exon", ]
  head(p_pos)
  p_pos <- p_pos[p_pos[["seqnames"]] == g_loc[["seqnames"]], ]
  nrow(p_pos)
  names(p_pos)
  unique(p_pos[["gene_name"]])
  p_peak <- p[p[["nearestGene"]] %in% unique(p_pos[["gene_name"]]), ]
  p_peak <- p_peak[!is.na(p_peak[["seqnames"]]), ]
  head(p_peak)
  unique(p_pos[, c("gene_name", "gene_id")])
  p_g_rng <- dplyr::bind_rows(
    setNames(
      lapply(
        seq.int(
          1,
          nrow(unique(p_pos[, c("gene_name", "gene_id")])),
          1
        ),
        function(x) {
          r <- unique(p_pos[, c("gene_name", "gene_id")])
          d <- p_pos[
            p_pos[["gene_name"]] == r[x, ][["gene_name"]] &
              p_pos[["gene_id"]] == r[x, ][["gene_id"]],
          ]
          d <- data.frame(
            "gene_name" = unique(d[["gene_name"]]),
            "gene_id" = unique(d[["gene_id"]]),
            "start" = min(d[["start"]]),
            "end" = max(d[["end"]])
          )
          return(d)
        }
      ),
      unique(p_pos[, c("gene_name", "gene_id")][["gene_id"]])
    )
  )
  # Plot
  p_seq <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = p_peak,
      ggplot2::aes(
        x = .data[["start"]], # nolint
        xend = .data[["end"]],
        y = 0
      ),
      linewidth = 12,
      alpha = 0.6,
      color = col_univ()[[20]]
    ) +
    ggplot2::geom_segment(
      data = p_pos[
        p_pos[["type"]] == "exon",
      ],
      ggplot2::aes(
        x = start,
        xend = end,
        color = g_name,
        y = 0
      ),
      linewidth = 8,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        xend = end,
        y = 0
      ),
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = p_g_rng,
      ggplot2::aes(
        x = start,
        y = -0.1,
        label = g_name,
        color = g_name
      ),
      bg.color = "white",
      show.legend = FALSE,
      size = 5
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.2, 0.2)
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(
        size = 14,
        angle = 90,
        color = "white"
      ),
      axis.line.x.bottom = ggplot2::element_line(color = "black"),
      axis.text.x = ggplot2::element_text(
        color = "grey40",
        size = 14,
        face = "bold"
      ),
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
      axis.ticks.length.x = grid::unit(0.2, "cm"),
      axis.ticks.x = ggplot2::element_line(color = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Position (bp)", y = "Seq") +
    ggplot2::scale_x_continuous(
      limits = c(
        (min(d_norm[["start"]]) - 5000),
        max(d_norm[["start"]] + 5000)
      )
    )

  p_cov_comb <- ggpubr::ggarrange(
    p_cov,
    p_seq,
    ncol = 1,
    nrow = 2,
    heights = c(0.9, 0.1)
  ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")
    )

  return(p_cov_comb)
}

#' Add JASPAR TF Motifs
#'
#' Creates a ChromatinAssay from a scATAC-Seq counts assay
#' and annotates transcription factor motifs based on
#' the JASPAR database.
#'
#' @param so An object of class Seurat. Must contain an ATAC assay.
#' @param dref Path to a .gtf file containing reference gene annotations.
#' @return An annotated ChromatinAssay object.
#' @examples
#'
#' # d_chrmasy <- sc_atac_motifs(
#' #   d,
#' #   rtracklayer::import(
#' #    "ref/gencode.v45.primary_assembly.annotation.gtf"
#' #    )
#' # )
#'
#' @export
sc_atac_motifs <- function(
  so,
  dref
) {
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Hsapiens.UCSC.hg38)
  ## Extract counts from Seurat ATAC assay
  d <- so
  Seurat::DefaultAssay(d) <- "ATAC"
  dc <- SeuratObject::GetAssayData(d, slot = "counts")
  rownames(dc) <-  paste(
    d@assays$ATAC@meta.features[["seqnames"]],
    paste(
      d@assays$ATAC@meta.features[["start"]],
      d@assays$ATAC@meta.features[["end"]],
      sep = "-"
    ),
    sep = ":"
  )
  ## Format gene locations and create ChromatinAssay
  ref1 <- dref
  gene_coords <- ref1[ref1$type == "gene"]
  gene_coords$gene_biotype <- gene_coords$gene_type
  GenomeInfoDb::seqlevelsStyle(gene_coords) <- "UCSC"
  gene_coords <- GenomeInfoDb::keepStandardChromosomes(
    gene_coords,
    pruning.mode = "coarse"
  )

  d1 <- Signac::CreateChromatinAssay(
    counts = dc,
    sep = c(":", "-"),
    verbose = TRUE,
    meta.data = d@meta.data
  )

  ## Use CORE, CNE, PBM, PBM_HLH, PHYLOFACTS, POLII, SPLICE, and PBM_HOMEO
  list_clt1 <- c(
    "CORE", "CNE", "PBM", "PBM_HLH",
    "PBM_HOMEO", "PHYLOFACTS", "POLII",
    "SPLICE"
  )
  list_clt <- setNames(
    lapply(
      list_clt1,
      function(x) {
        opt1 <- list(
          collection = x,
          tax_group = "vertebrates",
          all_versions = FALSE
        )
        return(opt1)
      }
    ),
    list_clt1
  )

  ## Query JASPAR and combine position frequency matrices
  list_pfm <- setNames(
    lapply(
      list_clt,
      function(y) {
        pfm1 <- TFBSTools::getMatrixSet(
          x = JASPAR2020, # nolint
          opts = y
        )
        return(pfm1)
      }
    ),
    list_clt1
  )

  list_pfm <- list_pfm[lengths(list_pfm) > 0]

  ## 4 JASPAR collections have PFMs
  list_pfm <- c(
    list_pfm[[1]],
    list_pfm[[2]],
    list_pfm[[3]],
    list_pfm[[4]]
  )

  ## Add motifs to chromatinassay
  d1 <- Signac::AddMotifs(
    object = d1,
    genome = BSgenome.Hsapiens.UCSC.hg38, # nolint
    pfm = list_pfm,
  )

  ## Save motif list
  tf_list <- Signac::GetMotifData(d1)
  tf_list <- data.frame(
    "ID" = colnames(tf_list),
    "name" = unlist(
      lapply(
        colnames(tf_list),
        function(x) {
          name(TFBSTools::getMatrixByID(JASPAR2020, ID = x))
        }
      )
    )
  )
  write.table(
    tf_list,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    file = "analysis/table.motif.names.jaspar2020.txt"
  )
  saveRDS(d1, "analysis/data.chromassay.wTFmotifs.rds")
  return(d1)

}