#### Load all libraries, plot themes, and import data files

#### Libraries ####

# Libraries

lib.list.gen <- c(
  "parallel","Seurat","ggplot2",
  "ggpubr","circlize","dplyr",
  "SingleCellExperiment","MAST","ComplexHeatmap","EnhancedVolcano","biomaRt","ggsci"
  )

invisible(
  lapply(
    lib.list.gen, 
    library, 
    character.only = TRUE
    )
  )


#### Color Schemes ####

col1 <- c(pal_npg("nrc")(10),
          pal_aaas("default")(10),
          pal_lancet("lanonc")(9),
          pal_jama("default")(6))

# For discrete scale plots

## Plots with more than 4 variables

col1a <- c(
  "#440154FF","#F0F921FF","#482173FF",
  "#FCD225FF","#433E85FF","#FDAD32FF",
  "#38598CFF","#F58C46FF","#2D708EFF",
  "#E76F5AFF","#25858EFF","#D5546EFF",
  "#1E9B8AFF","#C03A83FF","#2BB07FFF",
  "#A62098FF","#51C56AFF","#8707A6FF",
  "#85D54AFF","#6300A7FF","#C2DF23FF",
  "#3E049CFF","#FDE725FF","#0D0887FF"
  )

## Plots with 4 or less variables

col1b <- c(
  "#433E85FF","#25858EFF",
  "#51C56AFF","#FDE725FF"
  )


# For continuous scale plots

## Heatmaps

col2a <- viridis::viridis(12)

col3a <- c(
  "#395D9CFF","#91DBB4FF","#403872FF",
  "#240C4FFF","#ED6925FF","#DEF5E5FF",
  "#FAC62DFF","#AE305CFF","#3497A9FF",
  "#FCA50AFF","#000004FF","#2D1D38FF",
  "#DD513AFF","#F8850FFF","#45BDADFF",
  "#420A68FF","#5D126EFF","#FCFFA4FF",
  "#BBE7C8FF","#382A54FF","#40498EFF",
  "#3671A0FF","#38AAACFF","#3484A5FF",
  "#932667FF","#781C6DFF","#C73E4CFF",
  "#0B0405FF","#0C0826FF","#F2E661FF",
  "#60CEACFF","#1D111DFF"
  )

#### Plot themes ####

# Generic parameters shared among all plots

thm.gen <- theme(
  # Plot Title
  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 14
    ),
  
  # Panel
  panel.border = element_blank(),
  panel.background = element_blank(),
  panel.grid.major.y = element_line(
    colour = 'grey85'
    ),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  
  # Axes
  axis.ticks.y = element_blank(),
  axis.text.x = element_text(
    face = 'bold',
    size = 14,
    angle = 45,
    hjust = 1,
    vjust = 1
    ),
  axis.text.y = element_text(
    face = 'bold',
    size = 14
    ),
  axis.title.x = element_text(
    face = 'bold',
    size = 14
    ),
  axis.title.y = element_text(
    face='bold',
    size = 14
    ),
  
  # Strip
  strip.background = element_rect(
    fill = 'slategray2'
    ),
  strip.text = element_text(
    face = 'bold',
    size = 12
    ),
  # Margins
  plot.margin = unit(
    c(
      2,2,2,2
      ),
    "cm"
    )
  )


# legend parameters

## add legend title

thm.leg.title.y <- theme(
  legend.title = element_text(
    size = 16,
    face = 'bold'
    )
  )

thm.leg.title.n <- theme(
  legend.title = element_blank()
  )

thm.leg.main <- theme(
  legend.text = element_text(
    size = 14),
  legend.key.size = unit(
    0.4,
    'cm'
    ),
  legend.key = element_blank(),
  legend.justification = c(
    "right",
    "top"
    ),
  legend.box.just = "right"
  )


# Individual plot themes

## PCA, PLS-DA, UMAP

thm.mult <- thm.gen +
  thm.leg.title.y +
  thm.leg.main

## Boxplot, barplot, volcano plot, distribution plot

thm.univ <- thm.gen +
  thm.leg.title.n +
  thm.leg.main




