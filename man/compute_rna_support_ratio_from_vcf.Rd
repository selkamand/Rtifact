% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expect_support.R
\name{compute_rna_support_ratio_from_vcf}
\alias{compute_rna_support_ratio_from_vcf}
\title{Compute proportion of DNA mutations with RNA support}
\usage{
compute_rna_support_ratio_from_vcf(
  vcf,
  dna_sample,
  rna_sample,
  min_alt_supporting_reads = 2,
  confidence = 0.95,
  field_AF = "AF",
  field_DP = "DP",
  field_AD = "AD",
  pass_only = TRUE,
  return_data = FALSE,
  verbose = TRUE,
  exclude_sex_chromosomes = TRUE
)
}
\arguments{
\item{vcf}{path to a multi sample VCF which describes mutations in both DNA and RNA samples.
These can be produced by running the \href{https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md#append-mode}{sage variant caller in append mode}.}

\item{dna_sample}{the name of the tumour DNA sample in the VCF}

\item{rna_sample}{the name of the tumour RNA sample in the VCF}

\item{min_alt_supporting_reads}{the minimum number of alt supporting reads
required to consider a DNA variant supported by RNA (numeric)}

\item{confidence}{For any particular variant, the theoretical percentage of
finding >=\code{min_alt_supporting_reads} above which we are comfortable
classifying the variant as 'expected' to have RNA support.}

\item{field_AF}{Genotype field that describes allele frequency uncorrected for tumour purity/ploidy. (string)}

\item{field_DP}{Genotype field that describes read depth. (string)}

\item{field_AD}{Genotype field that describes allelic depth. (string)}

\item{pass_only}{only compute rna support ratio using variants where FILTER=PASS in VCF (flag)}

\item{return_data}{return the data.frame supplied to \code{\link[=compute_rna_support_ratio]{compute_rna_support_ratio()}} instead of the computed result. Mainly used for debug (flag)}

\item{verbose}{verbose (flag)}

\item{exclude_sex_chromosomes}{Exclude variants from sex chromosomes. See \code{\link[=sex_chromosomes]{sex_chromosomes()}} for all values of 'CHROM' we would identify as sex chromosomes. (flag)}
}
\value{
A list of class \code{rna_support} containing the following components:
\describe{
\item{\strong{total_variants_evaluated}}{Integer. The total number of variants evaluated for RNA support.}
\item{\strong{expected_rna_support_variants}}{Integer. The number of variants expected to show RNA support based on their DNA characteristics (e.g., variant allele frequency, read depth).}
\item{\strong{proportion_with_observed_support}}{Numeric. The proportion of variants expected to show RNA support where RNA support was actually observed. Sometimes called 'Psupp'}
\item{\strong{proportion_unexpected_with_observed_support}}{Numeric. The proportion of variants not expected to show RNA support that nonetheless had RNA support observed.}
}
}
\description{
In what proportion of the DNA mutations which theoretically should have RNA support
do we observe that support.
Expectation of RNA support is based on DNA allele frequency and RNA depth.
Observed RNA support is based on presense of >= \code{min_alt_supporting_reads}.
}
\examples{
vcf <- system.file("example.vcf", package = "Rtifact")
compute_rna_support_ratio_from_vcf(vcf, dna_sample = "DNA", rna_sample = "RNA")
}
