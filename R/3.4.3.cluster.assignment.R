#### Cell Type Prediction Score Distribution ####

## Return cell type predictions and assign highest ranked type to each cell

list.analysis[["cluster.proportions"]] <- dplyr::select(
  unique(
    dplyr::left_join(
      setNames(
        aggregate(
          list.analysis$pred.prop.summary[["Proportion"]],
          list(
            list.analysis$pred.prop.summary[["seurat.clusters"]]
          ),
          FUN = max
        ),
        c(
          "seurat.clusters",
          "Proportion"
        )
      ),
      list.analysis$pred.prop.summary[,c(
        "predicted.id",
        "Proportion"
      )],
      by = c(
        "Proportion"
      )
    )
  ),
  c(
    "seurat.clusters",
    "predicted.id",
    "Proportion"
  )
)


list.analysis[["cluster.assignments"]] <- dplyr::mutate(
  list.analysis$cluster.proportions,
  "predicted.id" = paste(
    seq(
      1:nrow(
        list.analysis$cluster.proportions
      )
    ),
    gsub(
      "FOXN4",
      "",
      gsub(
        "\\.",
        "",
        gsub(
          "^*..",
          "",
          list.analysis$cluster.proportions$predicted.id
        )
      )
    ),
    sep = "."
  ),
  "CellGroup" = as.factor(
    gsub(
      "FOXN4",
      "",
      gsub(
        "\\.",
        "",
        gsub(
          "^*..",
          "",
          list.analysis$cluster.proportions$predicted.id
        )
      )
    )
  )
)

### Change grouping columns to factors

list.analysis[["cluster.assignments"]][["predicted.id"]] <- factor(
  list.analysis[["cluster.assignments"]][["predicted.id"]],
  levels = c(
    gtools::mixedsort(
      list.analysis[["cluster.assignments"]][["predicted.id"]]
    )
  )
)

list.analysis[["cluster.assignments"]][["CellGroup"]] <- factor(
  list.analysis[["cluster.assignments"]][["CellGroup"]],
  levels = c(
    list.ct
  )
)


list.analysis[["cluster.assignments"]] <- dplyr::select(
  dplyr::left_join(
    list.analysis$`Predicted Clusters`@meta.data,
    list.analysis$cluster.assignments,
    by = "seurat.clusters"
  ),
  c("predicted.id.y",
    "CellGroup.y"
  )
)


## Add Cell Type and Cell Group columns to seurat object

list.analysis$`Predicted Clusters` <- AddMetaData(
  list.analysis$`Predicted Clusters`,
  list.analysis$cluster.assignments[["predicted.id.y"]],
  col.name = "CellType"
)

list.analysis$`Predicted Clusters` <- AddMetaData(
  list.analysis$`Predicted Clusters`,
  list.analysis$cluster.assignments[["CellGroup.y"]],
  col.name = "CellGroup"
)




