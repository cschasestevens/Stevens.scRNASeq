# Stevens.scRNASeq v1.2

Streamlined Processing and Analysis of Single-cell RNA-Sequencing Datasets

## Description

Utilizes Seurat in tandem with various R packages to perform processing and analysis of single-cell RNA-Sequencing (scRNA-Seq) datasets. The methods included in this package provide a seamless workflow for commonly used Seurat functions, statistical methods, and visualization of scRNA-Seq data. Most analyses can be run in parallel using intuitive functions to expedite time-consuming steps such as dataset integration and differential expression analysis. The package is compatible with Windows, Linux, or WSL2. However, analyses conducted in Windows default to sequential processing due to inherent stability issues of parallel processing in Windows.

* Important: Seurat objects created in Seurat v5 may not integrate counts data properly. A forthcoming update will address this issue while retaining backwards compatibility with Seurat v4 objects.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.3.1 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr,
        * rmarkdown
    * Imports: 
        * ggplot2,
        * dplyr,
        * ggsci,
        * ggrepel,
        * gtools,
        * SoupX,
        * scDblFinder,
        * SingleCellExperiment,
        * SummarizedExperiment,
        * MAST,
        * viridis,
        * BiocGenerics,
        * parallel,
        * reshape2,
        * ggpubr,
        * Seurat,
        * SeuratObject,
        * future,
        * circlize,
        * ComplexHeatmap,
        * magrittr,
        * EnhancedVolcano,
        * lazyeval,
        * topGO,
        * org.Hs.eg.db,
        * biomaRt,
        * shadowtext

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/Stevens.scRNASeq", 
  ref = "master", 
  build_vignettes = T
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
?Stevens.scRNASeq::sc.predict.clusters
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
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
