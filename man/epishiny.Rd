% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app.R
\name{epishiny}
\alias{epishiny}
\title{Main function to start EpiShiny}
\usage{
epishiny(...)
}
\description{
This application was built based on the python package
epiScanpy, the tool for single-cell epigenomic analysis (Danese et al. 2019).
The application also provide the basic analysis functions such as clustering
and plotting from the python package scanpy (Wolf, Angerer, and Theis 2018).
Therefore, this function supports the input for both scRNA-seq and scATAC-seq.
Note that The function requires no parameters and always returns NULL.
}
\examples{
\dontrun{
library(epiview)
epiview()
}

}
\references{
Danese, Anna, Maria L. Richter, David S. Fischer, Fabian J. Theis, and
Maria Colomé-Tatché. 2019. “EpiScanpy: Integrated Single-Cell Epigenomic
Analysis.” Preprint. Bioinformatics. https://doi.org/10.1101/648097.

Wolf, F. Alexander, Philipp Angerer, and Fabian J. Theis. 2018. “SCANPY:
Large-Scale Single-Cell Gene Expression Data Analysis.” Genome Biology 19
(1): 15. https://doi.org/10.1186/s13059-017-1382-0.
}
