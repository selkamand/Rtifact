

# Utilities from VCF filepath ---------------------------------------------


#' List All Genotype Fields
#'
#' @param vcf path to a vcf file
#'
#' @return a vector of VCF fields
#' @export
#'
#' @examples
#' vcf <- system.file("example.vcf", package = "Rtifact")
#' list_genotype_fields(vcf)
list_genotype_fields <- function(vcf){
  # Parse file
  vcf_obj <- VariantAnnotation::readVcf(file = vcf)

  #Genotype fields
  fields <- names(VariantAnnotation::geno(vcf_obj))

  fields_from_header <- VariantAnnotation::geno(
    VariantAnnotation::header(vcf_obj)
  )
  fields_with_descriptions <- paste0(rownames(fields_from_header), ": ", fields_from_header$Description)
  return(fields_with_descriptions)
}


# Utilies for VCF objects ---------------------------------------------------
get_samples <- function(vcf_obj){
  VariantAnnotation::samples(VariantAnnotation::header(vcf_obj))
}

is_multiallelic <- function(vcf_obj){
  lengths(VariantAnnotation::alt(vcf_obj)) > 1
}

is_pass <- function(vcf_obj){
  VariantAnnotation::filt(vcf_obj) == "PASS"
}


# Simulate ----------------------------------------------------------------

write_simulated_vcf <- function(outfile, nvariants=20){
  header <- '##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency per sample">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DNA	RNA'

  # Randomly generate the info we need
  n = nvariants
  pos = round(runif(n=n, min = 0, max = 1e9), digits = 0)
  chr = paste0("chr", sort(round(runif(n=n, min = 1, max = 22), digits = 0)))
  id="."
  ref = sample(size = n, c("A", "C", "T", "G"), replace = TRUE)
  alt = sample_alt(ref)
  qual = 20
  filter = sample(x=c("PASS", "FAIL"), size = n, replace = TRUE, prob = c(0.9, 0.1))
  info = "."
  format = "GT:AF:AD:DP"
  dna_gt = "0/1"
  rna_gt = "./."
  dna_af = round(runif(n=n, min = 0.1, max = 1), digits = 2)
  rna_af = round(runif(n=n, min = 0, max = 1), digits = 2)
  dna_dp = round(pmax(rnorm(n = n, mean = 20, sd = 6), 0), digits = 0)
  rna_dp = round(pmax(rnorm(n = n, mean = 5, sd = 6), 0), digits = 0)

  dna_ad_alt = round(dna_af * dna_dp, digits = 0)
  dna_ad_ref = dna_dp - dna_ad_alt
  dna_ad = paste(dna_ad_ref, dna_ad_alt, sep = ",")

  rna_ad_alt = round(rna_af * rna_dp, digits = 0)
  rna_ad_ref = rna_dp - rna_ad_alt
  rna_ad = paste(rna_ad_ref, rna_ad_alt, sep = ",")

  # Create Dataframe
  df <- data.frame(
    chrom = chr,
    pos = pos,
    id = id,
    ref = ref,
    alt = alt,
    qual = qual,
    filter = filter,
    info=info,
    format = format,
    dna_gt = paste(dna_gt, dna_af, dna_ad, dna_dp, sep = ":"),
    rna_gt = paste(rna_gt, rna_af, rna_ad, rna_dp, sep = ":")
    )

  # Write VCF
  file.create(outfile)
  writeLines(text = header, con = outfile)
  write.table(df, file = outfile, col.names = FALSE, quote = FALSE, append = TRUE, row.names = FALSE, sep = "\t")
}

sample_alt <- function(ref){
  vapply(
    X=ref,
    FUN = function(r){
        sample(size = 1, x = setdiff(c("A", "C", "T", "G"), r))
      },
    FUN.VALUE = character(1)
    )
}

#' Generates a bunch simulated VCFs
#'
#' Used to test running in cohort mode from manifest
#'
#' @param outfolder folder to write simulated VCFs to
#' @param prefix VCF name prefixes
#' @param nsamples how many VCFs to make
#' @param nvariants how many variants to write to each VCF
#' @param seed random seed
#' @return invisibly returns NULL. run for its side effects
#' @export
#'
#' @examples
#' if(interactive){
#' write_simulated_vcfs(outfolder="cohort")
#' }
write_simulated_vcfs <- function(outfolder = "cohort/", prefix = "simulated", nsamples = 20, nvariants = 20, seed = 111){

  if(dir.exists(outfolder)) stop("outfolder [", outfolder, "] already exists. Please remove and rerun `write_simulated_vcfs`")
  dir.create(outfolder)
  paths = paste0(outfolder, "/", prefix,'_', seq_len(nsamples), ".vcf")

  for (path in paths){
    with_temp_seed(seed = seed, {write_simulated_vcf(outfile = path, nvariants = nvariants)})
  }

  invisible(NULL)
}


with_temp_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- .Random.seed
  } else {
    old_seed <- NULL
  }

  set.seed(seed)
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  force(expr)
}

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
