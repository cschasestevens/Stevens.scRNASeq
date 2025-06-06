% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.1.heatmaps.R
\name{sc_dotplot}
\alias{sc_dotplot}
\title{Dot Plot}
\usage{
sc_dotplot(
  p_type = "threshold",
  topn = 10,
  l_gene,
  filt_var = NULL,
  so,
  asy1 = "SCT",
  ct_filt = FALSE,
  ct = NULL,
  ct_var,
  split_var = FALSE,
  list_var = NULL,
  col1 = col_grad(),
  vline1 = FALSE,
  xint = NULL
)
}
\arguments{
\item{p_type}{Should a custom gene list or DGEA results object
be used for plotting?
Type "cstm" for custom gene lists and "threshold"
for using DGEA/DA results.}

\item{topn}{If p_type is "threshold," select the number of
top DEGs/TFs to plot.}

\item{l_gene}{A list of DGEA/DA results
or a vector of selected genes for plotting.}

\item{filt_var}{If plotting based on a threshold,
which variable should be used for filtering?}

\item{so}{An object of class Seurat.}

\item{asy1}{Assay to use for the selected Seurat object. Note
that the assay must match the names in the seleceted
list type (ex. DGEA/DA).}

\item{ct_filt}{Filter based on cell type?}

\item{ct}{A pattern provided as a character string for
matching a specific cell type or types if ct_filt is TRUE.}

\item{ct_var}{Character string indicating the name of the
cluster/cell type variable.}

\item{split_var}{Logical indicating whether the cluster variable
should be stratified by additional group variables.}

\item{list_var}{If split_var is TRUE, provide a vector of character strings
containing two group variables for stratifying plot points.}

\item{col1}{Gradient color scheme to use.}

\item{vline1}{Add a vertical line to plot?}

\item{xint}{If vline1 is TRUE, specify the x-intercept of the vertical line
to place on the dot plot.}
}
\value{
A dotplot and formatted data frame for the selected assay.
}
\description{
Generates a dot plot from a Seurat Object based on a DGEA/DA list based on
a threshold or custom gene list for a specific cell type and comparison.
}
\examples{

# dot1 <- sc_dotplot( # nolint
#   p_type = "cstm",
#   l_gene = c("IRF1", "IRF2", "STAT1"),
#   so = d,
#   asy1 = "SCT",
#   ct_var = "CellType",
#   split_var = TRUE,
#   list_var = "Group",
#   vline1 = FALSE
# )

}
