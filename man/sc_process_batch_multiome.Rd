% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1.3.process.files.R
\name{sc_process_batch_multiome}
\alias{sc_process_batch_multiome}
\title{Batch Processing of CellRanger Multiome Files}
\usage{
sc_process_batch_multiome(
  df_par,
  p_name,
  m_cell,
  m_feat,
  atac_high,
  rna_high,
  rna_low,
  perc_mt,
  nsm_high,
  tss_low,
  parl,
  core_perc
)
}
\arguments{
\item{df_par}{Data frame of parameters to use for data processing.}

\item{p_name}{Character string providing the project name.}

\item{m_cell}{Minimum cells per feature.}

\item{m_feat}{Minimum features per cell.}

\item{atac_high}{Upper limit for per cell ATAC counts.}

\item{rna_high}{Upper limit for per cell RNA counts.}

\item{rna_low}{Lower limit for per cell RNA counts.}

\item{perc_mt}{Upper limit for percentage of mitochondrial reads per cell.}

\item{nsm_high}{Upper limit for nucleosome signal. Used to filter
potential doublets, potential artifacts, etc.}

\item{tss_low}{Lower limit for TSS enrichment score. Used to filter
low quality cells.}

\item{parl}{Logical indicating whether processing should be run in
parallel (Linux and WSL2 only). Set to FALSE if running sequentially.}

\item{core_perc}{Percentage of available cores to use if running
in parallel (Linux and WSL2 only). Set to 1 if running sequentially.}
}
\value{
A processed list of 10X multiome sample files converted
into Seurat objects with a summary list of QC and processing details.
}
\description{
Processes a list of 10X multiome samples for downstream integration.
}
\examples{

# list_data <- sc_process_batch_multiome(
#  # parameter list
#  l_params[["param"]],
#  # Project name
#  "Hiro.CF.multiome",
#  # Minimum cells per feature
#  5,
#  # Minimum features per cell
#  200,
#  # ATAC counts upper limit
#  100000,
#  # RNA counts upper limit
#  25000,
#  # RNA counts lower limit
#  1000,
#  # Percentage mitochondrial reads upper limit
#  20,
#  # Nucleosome signal upper limit
#  2.5,
#  # TSS enrichment lower limit
#  2,
#  # Run in parallel?
#  TRUE,
#  # Core percentage
#  0.1
# )

}
