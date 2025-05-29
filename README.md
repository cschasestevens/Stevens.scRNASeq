# Stevens.scRNASeq v3.02 (20250529)

Processing and Analysis of Single-cell RNA-Sequencing and Single-cell ATAC-Sequencing Datasets

## Description

Utilizes Seurat and Signac in tandem with various R packages to perform processing and analysis of single-cell RNA-Sequencing (scRNA-Seq) and single-cell ATAC-Sequencing (scATAC-Seq) datasets.    The methods included in this package provide a seamless workflow for commonly used Seurat and Signac functions, statistical methods, and visualization.    Most analyses can be run in parallel using intuitive functions to expedite time-consuming steps such as dataset integration and differential expression analysis.    The package is compatible with Windows, Linux, or WSL2. However, analyses conducted in Windows default to sequential processing due to inherent stability issues of parallel processing in Windows.    Version 3.00 adds functions for QC and analysis of Phenocycler segmentation results exported from QuPath.

* Important: Functions run in Linux or WSL2 require installation of additional dependencies and setup of a conda environment. Also, Phenocycler analysis is primarily conducted using QuPath and Python.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.4.3 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate
    * Imports: 
        * ggplot2,
        * viridis,
        * ggsci,
        * rtracklayer,
        * GenomeInfoDb,
        * Rsamtools,
        * SoupX,
        * scDblFinder,
        * SingleCellExperiment,
        * SummarizedExperiment,
        * Seurat,
        * SeuratObject,
        * BiocGenerics,
        * parallel,
        * decontX,
        * Signac,
        * BSgenome,
        * harmony,
        * future,
        * ggpubr,
        * plotly,
        * htmlwidgets,
        * ggrepel,
        * dplyr,
        * reshape2,
        * lazyeval,
        * magrittr,
        * circlize,
        * ComplexHeatmap,
        * grid,
        * gtools,
        * patchwork,
        * Azimuth,
        * BSgenome.Hsapiens.UCSC.hg38,
        * chromVAR,
        * JASPAR2020,
        * TFBSTools,
        * MAST,
        * EnhancedVolcano,
        * biomaRt,
        * topGO,
        * shadowtext,
        * stringr,
        * org.Hs.eg.db,
        * CellChat

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/Stevens.scRNASeq", 
  ref = "master", 
  build_vignettes = TRUE
)
```

## Help
* Browse vignettes by running the following:

```
browseVignettes("Stevens.scRNASeq")
```

* Access function documentation by running the following:

```
# Type function name after package name
?Stevens.scRNASeq::sc_top10_marker_heatmap
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
* 3.02
    * Added functions for file conversion from Seurat to .h5ad.
    * Added function for visualizing phenocycler cell maps.
    * Integrated phenocycler data format into existing functions for scRNA-Seq and multiome data analysis.
* 3.00
    * Added functions for QC and analysis of Phenocycler segmentation results.
* 2.4
    * Bug fixes for subsetting and reclustering multiome data.
    * Simplified parameters for UMAP and visualization.
    * Added function for Azimuth based cell type predictions.
    * Updated custom cell type prediction function.
    * Updated dot plot function to accept both gene expression and ATAC data.
* 2.3
    * Simplified plotting functions with default parameters.
    * Added batch plotting functions for analyses spanning multiple Seurat objects.
* 2.2
    * Integration of Cell Chat analysis and chord diagram plotting functions.
    * Revised scATAC-Seq coverage plot functions for greater plotting flexibility.
* 2.1
    * Added support for reclustering analysis of scRNA-Seq data
    * Added flexibility for heat map and UMAP functions
    * Addresses errors in package documentation
* 2.0
    * Added functions for processing 10X multiome data sets
    * Updated v1.2 functions for compatibility with Seurat v5.0 and higher
    * Provided additional flexibility for plot parameters
    * Added common multiome plot types, including coverage plots and gene track plots
* 1.2
    * Added functions for enrichment analysis (both Gene Ontology and custom gene sets)
* 1.0
    * Initial Release

## License

This project is licensed under the GNU General Public License Version 3 - see the LICENSE.md file for details

## Acknowledgments

* Seurat package: Hao, Y., Stuart, T., Kowalski, M.H. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024). https://doi.org/10.1038/s41587-023-01767-y
* ComplexHeatmap package: Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. <doi:10.1093/bioinformatics/btw313>
* Circlize package: Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
