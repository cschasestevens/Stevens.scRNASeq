% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.2.umap.R
\name{sc.umap.standard}
\alias{sc.umap.standard}
\title{Standard UMAP Plot}
\usage{
sc.umap.standard(so, md.var)
}
\arguments{
\item{so}{An object of class Seurat.}

\item{md.var}{A character string indicating the clustering column for overlaying on a UMAP plot.}
}
\value{
A UMAP plot with points grouped by a specific metadata column.
}
\description{
Generates a single UMAP plot with a specific metadata overlay for a Seurat object.
}
\examples{

# p.umap <- sc.umap.standard(d.integrated,"col1")

}
