
#' Probability a DNA mutation will have RNA support
#'
#' @param DNA_AF the unadjusted allele frequency of the variant (numeric)
#' @param RNA_DP the total number of RNA reads covering the Loci (numeric)
#' @param min_alt_supporting_reads the minimum number of alt supporting reads (numeric)
#'
#'
#' @return a vector of probabilities (proportions)
#' @export
#'
#' @examples
#' # Single Variant
#' probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 1, min_alt_supporting_reads = 1)
#'
#' # Multiple Variants
#' probability_of_rna_support(
#'   DNAvaf = c(0.5, 0.5, 0.1),
#'   RNAdepth = c(1, 2, 2),
#'   min_alt_supporting_reads = c(1, 1, 2)
#' )
#'
#' # Multiple variants using the same setting for 'min_alt_supporting_reads'
#' probability_of_rna_support(
#'   DNAvaf = c(0.5, 0.5),
#'   RNAdepth = c(1, 2),
#'   min_alt_supporting_reads = 1
#' )
probability_of_rna_support <- function(DNA_AF, RNA_DP, min_alt_supporting_reads = 2){
  if(any(min_alt_supporting_reads <= 0)) stop("min_alt_supporting_reads must be greater than zero")

  # Calculate the probability of observing at least min_alt_supporting_reads
  # given a particular RNAdepth and the DNA VAF
  pbinom(q=min_alt_supporting_reads-1, size=RNA_DP, prob = DNA_AF, lower.tail = FALSE)
}



artifact_annotate_dataframe <- function(
    data,
    col_DNA_AF = "DNA_AF",
    col_RNA_DP = "RNA_DP",
    min_alt_supporting_reads = 2,
    confidence = 0.95){

  data[["probability_of_rna_support"]] <- probability_of_rna_support(
    DNA_AF = data[[col_DNA_AF]],
    RNA_DP = data[[col_RNA_DP]],
    min_alt_supporting_reads = min_alt_supporting_reads
  )

  data[["expecting_rna_support"]] <- data[["probability_of_rna_support"]] > confidence

  return(data)
}


# returns Nsupported/Nexpected
#
#' Title
#'
#' @param data A dataframe representing all autosomal mutations present in 1 sample.
#' Each row should represent a different mutation. Must the following columns (order does NOT matter, just column names)
#' 1. \strong{col_DNA_AF} Raw allele frequency in DNA (not adjusted for purity/ploidy). AF in oncoanalyser.
#' 2. \strong{col_RNA_DP} Total depth of coverage at loci in RNA. DP in oncoanalyser.
#' 3. \strong{col_RNA_AD} Number of reads supporting ALT allele in RNA. AD in oncoanalyser.
#'
#'
#' @param col_DNA_AF name of DNA_AF column (string)
#' @param col_RNA_DP name of RNA_DP column (string)
#' @param col_RNA_AD name of RNA_AD column (string)
#' @param min_alt_supporting_reads the minimum number of alt supporting reads
#' required to consider a DNA variant supported by RNA (numeric)
#'
#' @param confidence For any particular variant, the theoretical percentage of
#' finding >`min_alt_supporting_reads` above which we are comfortable
#' classifying the variant as 'expected' to have RNA support.
#'
#' @return \strong{Psupp}: The proportion of variants in the sample where we
#' expected to see  RNA evidence for and actually saw that RNA evidence in practice.
#'
#' @export
#'
#' @examples
#'
#' Create a sample data frame representing mutations in a single sample
#' sample_data <- data.frame(
#'   DNA_AF = c(0.5, 0.3, 0.2, 0.1),  # DNA variant allele frequencies
#'   RNA_DP = c(100, 50, 20, 10),     # RNA read depths at each locus
#'   RNA_AD = c(48, 15, 4, 0)         # RNA reads supporting the alternate allele
#' )
#'
#' # Calculate Psupp for the sample
#' Psupp <- artifact(
#'   data = sample_data,
#'   col_DNA_AF = "DNA_AF",
#'   col_RNA_DP = "RNA_DP",
#'   col_RNA_AD = "RNA_AD",
#'   min_alt_supporting_reads = 2,
#'   confidence = 0.95
#' )
artifact <- function(
    data,
    col_DNA_AF = "DNA_AF",
    col_RNA_DP = "RNA_DP",
    col_RNA_AD = "RNA_AD",
    min_alt_supporting_reads = 2,
    confidence = 0.95
  ){

  prob_rna_support <- probability_of_rna_support(
    DNA_AF = data[[col_DNA_AF]],
    RNA_DP = data[[col_RNA_DP]],
    min_alt_supporting_reads = min_alt_supporting_reads
  )

  # Assertions
  cnames <- colnames(data)
  if(!col_DNA_AF %in% cnames) stop("Could not find column: [", col_DNA_AF, "] in `data`.")
  if(!col_RNA_DP %in% cnames) stop("Could not find column: [", col_RNA_DP, "] in `data`.")
  if(!col_RNA_AD %in% cnames) stop("Could not find column: [", col_RNA_AD, "] in `data`.")

  expecting_rna_support = prob_rna_support > confidence
  rna_supported = data[[col_RNA_AD]] >= min_alt_supporting_reads

  n_expecting_support = sum(expecting_rna_support)
  n_supported = sum(rna_supported)

  n_expectedly_supported = sum(expecting_rna_support & rna_supported)
  n_unxpectedly_supported = sum(!expecting_rna_support & rna_supported)

  Psupp <- n_expectedly_supported/n_expecting_support

  return(Psupp)
}
