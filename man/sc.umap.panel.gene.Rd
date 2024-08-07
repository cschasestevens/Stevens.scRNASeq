% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.2.umap.R
\name{sc.umap.panel.gene}
\alias{sc.umap.panel.gene}
\title{Visualize Individual Gene Expression}
\usage{
sc.umap.panel.gene(so, md.var, g.name, col.scheme, col.names, leg.x, leg.y)
}
\arguments{
\item{so}{An object of class Seurat.}

\item{md.var}{A character string indicating the clustering column for overlaying on a UMAP plot.}

\item{g.name}{A character string indicating the gene name.}

\item{col.scheme}{The color scheme to be used for distinguishing between groups, provided as a vector.}

\item{col.names}{A vector of the same length as the provided color scheme for assigning colors to each group.}

\item{leg.x}{A numeric value indicating the placement of the figure legend on the x-axis.}

\item{leg.y}{A numeric value indicating the placement of the figure legend on the y-axis.}
}
\value{
A series of plots stored as a ggplot2 object for visualizing cluster gene expression.
}
\description{
Generates a series of plots to visualize individual gene expression per cluster.
}
\examples{

# p.umap <- sc.umap.panel(d.integrated,c("col1","col2","col3"))

}
