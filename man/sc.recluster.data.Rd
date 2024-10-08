% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.2.recluster.R
\name{sc.recluster.data}
\alias{sc.recluster.data}
\title{Reclustering Analysis}
\usage{
sc.recluster.data(
  so,
  ct,
  g.list,
  md.list,
  h.w,
  h.h,
  fs.c,
  fs.r,
  parl,
  core.perc,
  slot1,
  slot2,
  run.qc
)
}
\arguments{
\item{so}{An object of class Seurat. Must contain the columns 'CellType' and 'CellGroup' in the metadata slot.}

\item{ct}{Character string vector containing the name(s) of up to 2 cell type groups present in the CellGroup column.
The provided Seurat object is then subsetted to include only clusters from the provided group names.}

\item{g.list}{Path to an existing text file containing genes for plotting expression.}

\item{md.list}{A vector of character strings indicating metadata columns for overlaying on a loadings plot.}

\item{h.w}{Numeric value for marker gene heatmap width (passed to ComplexHeatmap).}

\item{h.h}{Numeric value for marker gene heatmap height (passed to ComplexHeatmap).}

\item{fs.c}{Numeric value for marker gene heatmap column fontsize (passed to ComplexHeatmap).}

\item{fs.r}{Numeric value for marker gene heatmap row fontsize (passed to ComplexHeatmap).}

\item{parl}{Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.}

\item{core.perc}{Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.}

\item{slot1}{A character string indicating which assay should be used for dimension reduction.}

\item{slot2}{A character string indicating the reduction name for plotting gene expression UMAPs.}

\item{run.qc}{A logical indicating whether gene expression QC plots should be generated for the reclustered data.}
}
\value{
A reclustered Seurat Object with summary UMAP plots and a marker gene list.
}
\description{
Performs reclustering and basic visualization of a subsetted Seurat Object. Subsets are determined based on CellGroup.
}
\examples{

# p.umap <- sc.top10.marker.heatmap(d.annotated,"seurat.clusters",18,24,6,8)

}
