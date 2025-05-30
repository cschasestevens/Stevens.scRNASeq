% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.5.phenocycler.QC.R
\name{pc_qc_sum}
\alias{pc_qc_sum}
\title{Phenocycler Normalization Summary}
\usage{
pc_qc_sum(type1 = "norm", so, samp_name, mcc = 2)
}
\arguments{
\item{type1}{Data type; provide a character string indicating
if the data have been normalized.}

\item{so}{Input Seurat object.}

\item{samp_name}{Sample name.}

\item{mcc}{Number of cores to use.}
}
\value{
A summary of normalization performance for the
selected Seurat Object.
}
\description{
Helper function for generating a short summary and images for
a normalized Phenocycler dataset.
}
\examples{

# d <- pc_qc_sum(
#   type1 = "raw",
#   so = dser,
#   samp_name = "sample1"
# )

}
