
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
#' probability_of_rna_support(DNA_AF = 0.5, RNA_DP = 1, min_alt_supporting_reads = 1)
#'
#' # Multiple Variants
#' probability_of_rna_support(
#'   DNA_AF = c(0.5, 0.5, 0.1),
#'   RNA_DP = c(1, 2, 2),
#'   min_alt_supporting_reads = c(1, 1, 2)
#' )
#'
#' # Multiple variants using the same setting for 'min_alt_supporting_reads'
#' probability_of_rna_support(
#'   DNA_AF = c(0.5, 0.5),
#'   RNA_DP = c(1, 2),
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



#' Compute proportion of DNA mutations with RNA support
#'
#' In what proportion of the DNA mutations which theoretically should have RNA support
#' do we observe that support.
#' Expectation of RNA support is based on DNA allele frequency and RNA depth.
#' Observed RNA support is based on presense of >= `min_alt_supporting_reads`.
#'
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
#' finding >=`min_alt_supporting_reads` above which we are comfortable
#' classifying the variant as 'expected' to have RNA support.
#'
#' @return \strong{Psupp}: The proportion of variants in the sample where we
#' expected to see  RNA evidence for and observed that RNA evidence in practice.
#'
#' @export
#'
#' @examples
#'
#' # Create a sample data frame representing mutations in a single sample
#' sample_data <- data.frame(
#'   DNA_AF = c(0.5, 0.3, 0.2, 0.1),  # DNA variant allele frequencies
#'   RNA_DP = c(100, 50, 20, 10),     # RNA read depths at each locus
#'   RNA_AD = c(48, 15, 4, 0)         # RNA reads supporting the alternate allele
#' )
#'
#' # Calculate Psupp for the sample
#' Psupp <- compute_rna_support_ratio(
#'   data = sample_data,
#'   col_DNA_AF = "DNA_AF",
#'   col_RNA_DP = "RNA_DP",
#'   col_RNA_AD = "RNA_AD",
#'   min_alt_supporting_reads = 2,
#'   confidence = 0.95
#' )
compute_rna_support_ratio <- function(
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

#' Compute proportion of DNA mutations with RNA support
#'
#' @param vcf path to a multi sample VCF which describes mutations in both DNA and RNA samples.
#' These can be produced by running the \href{https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md#append-mode}{sage variant caller in append mode}.
#' @param dna_sample the name of the tumour DNA sample in the VCF
#' @param rna_sample the name of the tumour RNA sample in the VCF
#' @param field_AF Genotype field that describes allele frequency uncorrected for tumour purity/ploidy. (string)
#' @param field_DP Genotype field that describes read depth. (string)
#' @param field_AD Genotype field that describes allelic depth. (string)
#' @param pass_only only compute rna support ratio using variants where FILTER=PASS in VCF (flag)
#' @param return_data return the data.frame supplied to [compute_rna_support_ratio()] instead of the computed result. Mainly used for debug (flag)
#' @param verbose verbose (flag)
#' @inherit compute_rna_support_ratio_from_vcf description return
#' @inheritParams compute_rna_support_ratio
#' @export
#'
#' @examples
#' vcf <- system.file("example.vcf", package = "Rtifact")
#' compute_rna_support_ratio_from_vcf(vcf, dna_sample = "DNA", rna_sample = "RNA")
compute_rna_support_ratio_from_vcf <- function(
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
    verbose=TRUE#,
    #exclude_chromosomes = sex_chromosomes()
    ){

  # Assertions
  if(length(vcf) > 1) stop("To calculate expected RNA support from multiple VCFs, try `compute_rna_support_ratio_from_manifest()`")
  if(!file.exists(vcf)) stop("Could not find VCF file: [", vcf, "]")
  # Parse file
  vcf_obj <- VariantAnnotation::readVcf(file = vcf)

  #Genotype fields
  fields <- names(VariantAnnotation::geno(vcf_obj))

  if(!field_AF %in% fields) stop(
  "Could not find the allele-frequency field: [",field_AF ,"] in VCF header.
   To see a list of valid genotype fields, run `list_genotype_fields(",vcf,")`")

  if(!field_DP %in% fields) stop(
    "Could not find the depth field: [",field_DP ,"] in VCF header.
   To see a list of valid genotype fields, run `list_genotype_fields(",vcf,")`")

  if(!field_AD %in% fields) stop(
    "Could not find the allelic depth field: [",field_AD ,"] in VCF header.
   To see a list of valid genotype fields, run `list_genotype_fields(",vcf,")`")

  # Check dna_sample & rna_sample are in the VCF
  samples <- get_samples(vcf_obj)

  if(!dna_sample %in% samples) stop(
    "Failed to find DNA sample [", dna_sample,"] in VCF. \nValid samples include: [",
    paste0(samples, collapse = ", "),
    "]"
  )

  if(!rna_sample %in% samples) stop(
    "Failed to find RNA sample [", rna_sample,"] in VCF. \nValid samples include: [",
    paste0(samples, collapse = ", "),
    "]")

  # Filter for PASS variants
  if(pass_only){
    vcf_obj <- vcf_obj[is_pass(vcf_obj)]
  }

  # Exclude multiallelic loci
  vcf_obj <- vcf_obj[!is_multiallelic(vcf_obj)]
  n_multiallelics <- sum(is_multiallelic(vcf_obj))
  if(verbose & n_multiallelics > 0) message("Found: ", n_multiallelics, " multiallelic sites which will be ignored by the Rtifact algorithm")

  # Extract the genotype Information we need
  geno <- VariantAnnotation::geno(vcf_obj)

  dna_af <- vapply(geno[[field_AF]][,dna_sample], \(x) {x[1]}, FUN.VALUE = numeric(1))
  rna_depths <- vapply(geno[[field_DP]][,rna_sample], \(x) {x[1]}, FUN.VALUE = numeric(1))

  # Pulls out the 'alt' support (second value in AD). This works because there will be no multi-allelics - they have all been excluded
  rna_alt_support <- vapply(geno[[field_AD]][,rna_sample], \(x){ x[2] }, FUN.VALUE = numeric(1))

  # Checks
  if(any(is.na(rna_depths))) stop("Problem parsing the total depth  genotype field [", field_DP ,"]. Are you sure `field_DP` is set to an appropriate genotype field?")
  if(any(is.na(rna_alt_support))) stop("Problem parsing the RNA allelic depth genotype field [", field_AD ,"]. Are you sure `field_AD` is set to an appropriate genotype field?")
  if(any(is.na(dna_af))) stop("Problem parsing the DNA allele frequency field [", field_AF ,"]. Are you sure `field_AF` is set to an appropriate genotype field?")

  if(!all(is_wholenumber(rna_depths))) stop("Total depth genotype field must contain only whole numbers, however the field [",field_DP,"] breaks this rule. Are you sure `field_DP` is set to the right field?")
  if(!all(is_wholenumber(rna_alt_support))) stop("The RNA allelic depth genotype field must contain only whole-number read support for each allele, however the field [", field_AD ,"] breaks this rule. Are you sure `field_AD` is set to the right field?")


  # Turn into a dataframe compatible with 'artifact' function
  data = data.frame(
    DNA_AF = dna_af,
    RNA_DP = rna_depths,
    RNA_AD = rna_alt_support
  )

  if(return_data){
    return(data)
  }

  # Compute Degree of RNA support
  compute_rna_support_ratio(data, min_alt_supporting_reads = min_alt_supporting_reads, confidence = confidence)

}

#' Compute RNA support from multiple VCFs
#'
#' @param manifest path to tsv with the following columns
#' 1. vcf: path to a multi sample VCF which describes mutations in both DNA and RNA.
#' These are typically produced by calling the \href{https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md#append-mode}{sage variant caller in append mode}.
#' 2. dna_sample: the name of the tumour DNA sample in the VCF
#' 3. rna_sample: the name of the tumour RNA sample in the VCF
#' @param sep delimiter used in manifest file.
#' @return the manifest dataframe with an extra column 'rna_support_ratio' describing
#' \strong{Psupp}: The proportion of variants in the sample where we expected to
#' see  RNA evidence for and observed that RNA evidence in practice.
#' @export
#'
#' @examples
#' if(interactive()){
#'   manifest <- system.file("manifest.tsv", package = "Rtifact")
#'   write_simulated_vcfs("cohort")
#'   compute_rna_support_ratio_from_manifest(
#'     manifest,
#'     field_AF = "AF",
#'     field_DP = "DP",
#'     field_AD = "AD",
#'     min_alt_supporting_reads = 2,
#'     confidence = 0.95
#'   )
#' }
#'
#' @inherit compute_rna_support_ratio_from_vcf description
#' @inheritParams compute_rna_support_ratio_from_vcf
#' @inheritParams compute_rna_support_ratio
compute_rna_support_ratio_from_manifest <- function(
    manifest,
    min_alt_supporting_reads = 2,
    confidence = 0.95,
    field_AF = "AF",
    field_DP = "DP",
    field_AD = "AD",
    pass_only = TRUE,
    sep = "\t",
    verbose=TRUE
    ){
  if(!file.exists(manifest)) stop("Could not find manifest file: [", manifest,"]")

  df_manifest <- utils::read.delim(file = manifest, header = TRUE, sep = sep)

  cnames <- colnames(df_manifest)
  if(!"vcf" %in% cnames) stop("Manifest file must include the column: 'vcf'")
  if(!"dna_sample" %in% cnames) stop("Manifest file must include the column: 'dna_sample'")
  if(!"rna_sample" %in% cnames) stop("Manifest file must include the column: 'rna_sample'")
  if(!nrow(df_manifest) > 0) stop("Manifest is empty")

  df_manifest[["rna_support_ratio"]] <- vapply(
    X = seq_len(nrow(df_manifest)),
    FUN = \(i){
      sample <- df_manifest[i, , drop = FALSE]
      compute_rna_support_ratio_from_vcf(
        vcf = sample[["vcf"]],
        dna_sample =  sample[["dna_sample"]],
        rna_sample =  sample[["rna_sample"]],
        min_alt_supporting_reads = min_alt_supporting_reads,
        confidence = confidence,
        field_AF = field_AF,
        field_DP = field_DP,
        field_AD = field_AD,
        pass_only = TRUE,
        verbose = verbose
        )
      },
    FUN.VALUE = numeric(1)
  )

  return(df_manifest)

}


