% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.4.integrate.files.R
\name{sc_integrate_data}
\alias{sc_integrate_data}
\title{Batch Integration of scRNA-Seq and Multiome Datasets}
\usage{
sc_integrate_data(l_so, l_par, proc_mode, parl, core_perc)
}
\arguments{
\item{l_so}{List of processed Seurat objects to be integrated.}

\item{l_par}{List of processing parameters passed to function.
For multiome data integration, must have a gene annotation and
reference genome file.}

\item{proc_mode}{Processing mode to be used
(either "scRNA-Seq" or "multiome").}

\item{parl}{Logical indicating whether processing should be run
in parallel. Set to FALSE if running sequentially.}

\item{core_perc}{Proportion (as a numeric value) of available cores
to use if running in parallel. Set to 1 if running sequentially.}
}
\value{
A Seurat object containing integrated data for all samples
present in a scRNA-Seq or multiome experiment.
}
\description{
Integrates a list of samples processed as Seurat objects.
}
\examples{

# d_integrated <- sc_integrate_data(
#   l_so = list_data,
#   l_par = l_params,
#   proc_mode = "multiome",
#   parl = TRUE,
#   core_perc = 0.5
# )

}
