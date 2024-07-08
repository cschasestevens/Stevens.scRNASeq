#### Load all libraries, plot themes, and import data files

#### Libraries ####

# Libraries

lib.list.gen <- c(
  "parallel","Seurat","ggplot2",
  "ggpubr",
  "gg3D",
  "ggrepel",
  "shadowtext","RColorBrewer","circlize","ggsci","presto",
  "future"
  )

invisible(
  lapply(
    lib.list.gen, 
    library, 
    character.only = TRUE
    )
  )





#### Color Schemes ####

# For discrete scale plots

## Plots with more than 4 variables

col1 <- c(pal_npg("nrc")(10),
          pal_aaas("default")(10),
          pal_lancet("lanonc")(9),
          pal_jama("default")(6))

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
  sample(
    viridis::mako(16),
    16
    ),
  sample(
    viridis::inferno(16),
    16
    )
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




