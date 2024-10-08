% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{list_genotype_fields}
\alias{list_genotype_fields}
\title{List All Genotype Fields}
\usage{
list_genotype_fields(vcf)
}
\arguments{
\item{vcf}{path to a vcf file}
}
\value{
a vector of VCF fields
}
\description{
List All Genotype Fields
}
\examples{
vcf <- system.file("example.vcf", package = "Rtifact")
list_genotype_fields(vcf)
}
