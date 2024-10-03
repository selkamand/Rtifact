
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rtifact

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/Rtifact)](https://CRAN.R-project.org/package=Rtifact)

<!-- badges: end -->

Rtifact leverages matched whole-genome sequencing (WGS) and RNA
sequencing (RNA-seq) data to identify samples with an unexpectedly high
number of DNA mutations lacking RNA support, indicating potential
artificial mutations or DNA sample contamination.

## Installation

You can install the development version of Rtifact like so:

``` r
if (!require("remotes"))
    install.packages("remotes")

remotes::install_github("selkamand/Rtifact")
```

## Quick Start

### Run on a single VCF files

The starting point of most analysis will be somatic tumour VCFs where
DNA variants have also been genotyped in a matched RNA sample.

[Oncoanalyser](https://nf-co.re/oncoanalyser/) will automatically
produced this VCF so long as both DNA and RNA fastq files are supplied
(see `<tumor_dna_id>.purple.somatic.vcf.gz` file)

``` r
# Load library
library(Rtifact)

# Run from VCF
vcf <- system.file("example.vcf", package = "Rtifact")
compute_rna_support_ratio_from_vcf(
  vcf, 
  dna_sample = "DNA", # Name of tumour DNA sample 
  rna_sample = "RNA"  # Name of tumour RNA sample
)
#> RNA Support Evaluation Result:
#> ----------------------------
#> Total Variants:  3 
#> [Nexpected] Variants Expected to have RNA Support:  2 
#> [Psupp] Observed rate of support: 50.0% (1/2) 
#> ----------------------------
```

### Run on cohort of VCF files

``` r
# Run on a cohort using a manifest file
manifest <- system.file("manifest.tsv", package = "Rtifact")

# Lets create the VCFs we want to test on
write_simulated_vcfs("cohort")

compute_rna_support_ratio_from_manifest(
  manifest, 
  field_AF = "AF",
  field_DP = "DP",
  field_AD = "AD",
  min_alt_supporting_reads = 2,
  confidence = 0.95
)
#>                        vcf dna_sample rna_sample total_variants_evaluated
#> 1   cohort/simulated_1.vcf        DNA        RNA                       17
#> 2  cohort/simulated_10.vcf        DNA        RNA                       18
#> 3  cohort/simulated_11.vcf        DNA        RNA                       20
#> 4  cohort/simulated_12.vcf        DNA        RNA                       16
#> 5  cohort/simulated_13.vcf        DNA        RNA                       20
#> 6  cohort/simulated_14.vcf        DNA        RNA                       18
#> 7  cohort/simulated_15.vcf        DNA        RNA                       16
#> 8  cohort/simulated_16.vcf        DNA        RNA                       18
#> 9  cohort/simulated_17.vcf        DNA        RNA                       19
#> 10 cohort/simulated_18.vcf        DNA        RNA                       19
#> 11 cohort/simulated_19.vcf        DNA        RNA                       19
#> 12  cohort/simulated_2.vcf        DNA        RNA                       20
#> 13 cohort/simulated_20.vcf        DNA        RNA                       17
#> 14  cohort/simulated_3.vcf        DNA        RNA                       17
#> 15  cohort/simulated_4.vcf        DNA        RNA                       18
#> 16  cohort/simulated_5.vcf        DNA        RNA                       17
#> 17  cohort/simulated_6.vcf        DNA        RNA                       17
#> 18  cohort/simulated_7.vcf        DNA        RNA                       20
#> 19  cohort/simulated_8.vcf        DNA        RNA                       20
#> 20  cohort/simulated_9.vcf        DNA        RNA                       16
#>    expected_rna_support_variants proportion_unexpected_with_observed_support
#> 1                             11                                   0.3333333
#> 2                              7                                   0.4545455
#> 3                              2                                   0.2777778
#> 4                              6                                   0.4000000
#> 5                             10                                   0.7000000
#> 6                              8                                   0.0000000
#> 7                              8                                   0.2500000
#> 8                              3                                   0.5333333
#> 9                              5                                   0.5000000
#> 10                             5                                   0.3571429
#> 11                             7                                   0.3333333
#> 12                             8                                   0.2500000
#> 13                             5                                   0.4166667
#> 14                             5                                   0.2500000
#> 15                             8                                   0.3000000
#> 16                             6                                   0.3636364
#> 17                            10                                   0.4285714
#> 18                             8                                   0.4166667
#> 19                             6                                   0.5000000
#> 20                             5                                   0.0000000
#>    proportion_with_observed_support
#> 1                         0.7272727
#> 2                         0.7142857
#> 3                         1.0000000
#> 4                         0.8333333
#> 5                         0.9000000
#> 6                         0.7500000
#> 7                         0.7500000
#> 8                         0.6666667
#> 9                         1.0000000
#> 10                        1.0000000
#> 11                        0.8571429
#> 12                        0.8750000
#> 13                        0.8000000
#> 14                        1.0000000
#> 15                        0.8750000
#> 16                        1.0000000
#> 17                        0.9000000
#> 18                        0.8750000
#> 19                        0.6666667
#> 20                        1.0000000
```

### Run from dataframe

``` r
# Create a sample data frame representing mutations in a single sample
sample_data <- data.frame(
  DNA_AF = c(0.5, 0.3, 0.2, 0.1),  # DNA variant allele frequencies
  RNA_DP = c(100, 50, 20, 10),     # RNA read depths at each locus
  RNA_AD = c(48, 15, 4, 0)         # RNA reads supporting the alternate allele
)

# Calculate RNA support evaluation
compute_rna_support_ratio(
  data = sample_data,
  col_DNA_AF = "DNA_AF",
  col_RNA_DP = "RNA_DP",
  col_RNA_AD = "RNA_AD",
  min_alt_supporting_reads = 2,
  confidence = 0.95
)
#> RNA Support Evaluation Result:
#> ----------------------------
#> Total Variants:  4 
#> [Nexpected] Variants Expected to have RNA Support:  2 
#> [Psupp] Observed rate of support: 100.0% (2/2) 
#> ----------------------------
```

## The artifact algorithm.

<u>**For each sample**</u> (with a substantial number of mutations)

1.  **Identify Mutations**: Select all autosomal, exonic mutations

2.  **Calculate Expected RNA Support:** Determine which DNA mutations
    are expected to have RNA support (by default, ≥2 supporting reads)
    based on DNA variant allele frequency $\mathrm{DNA}_{VAF}$ and RNA
    depth $\mathrm{RNA}_{DP}$. Denote this count as $N_{expected}$. See
    section: \[**Classifying mutations based on whether RNA support is
    expected**\].

3.  **Assess Actual Support**: Count how many of the mutations that
    theoretically should have RNA support actually have support
    ($N_{supp}$)

4.  **Compute** $P_{supp}$: Calculate the proportion of expected
    variants we actually find support for.

$$
P_{supp} = \frac{N_{supp}}{N_{expected}}
$$

<u>**Across the cohort**</u>

Repeat analysis for each sample in a cohort and identify the
distribution of $P_{supp}$. This distribution is required to identify
whether each sample has a ‘typical’ or ‘atypical’ proportion of DNA
variants with RNA support

### Interpreting RNA support level

1.  Quantify outlieryness of $P{supp}$ for each sample
2.  Classify samples based on degree of RNA support (relative to the
    average). **RNA support level** should be **typical** or **low**.

**Typical RNA Support**:

- Indicates DNA mutations are likely genuine. At the very least, DNA
  mutations were not artificially induced during library prep.

**Low RNA Support** may suggest:

1.  Absence of tumor cells in the RNA sample.
2.  Contamination of genomic sample (but NOT the RNA sample) with
    non-patient DNA.
3.  Artificial DNA mutations from the library preparation

### Classifying mutations based on whether RNA support is expected

To determine if a DNA variant should ‘theoretically’ have RNA support,
we must consider:

1.  **DNA VAF:** The variant’s allele frequency in DNA (not corrected
    for purity/ploidy)
2.  **RNA Depth:** The read depth at the variant’s locus in RNA.
3.  **Mutation Type:** Certain mutations (e.g., protein-truncating
    variants) may not be expressed due to mechanisms like
    nonsense-mediated decay \[Not added yet\].

Assuming no allele-specific or cell-specific heterogeneity in
expression, a variant’s expected RNA support is calculated using the
binomial distribution:

**Probability Calculation:**

``` math
\displaylines{
\textbf{Probability Density Function}\\
P(X=r) = nCr \times p^r \times (1 − p)^{n-r} \\~\\
\underline{Where: } \\
\textbf{n}\text{ - Number of events (RNA depth)} \\
\textbf{r}\text{ - Number of required successes (Number of alt supporting RNA reads)} \\
\textbf{p}\text{ - Probability of one success (DNA}_{VAF}\text{)} \\
\textbf{nCr}\text{ – Number of combinations ("n choose r")} \\
\textbf{P(X=r)}\text{ – Probability of an exact number of successes happening. (i.e getting an exact number of Alt supporting reads in RNA)}
\\~\\
\textbf{Cumulative Distribution Function}\\
P(AD_{alt} \ge 2) = 1 - P(AD_{alt} = 0) - P(AD_{alt} = 1)
}
```

- $P(AD_{alt} \ge 2)$: Probability of observing at least 2 alternate
  allele reads.

- Calculated using the cumulative binomial probability function (e.g.,
  `pbinom()` in R).

**Classification Threshold:** Variants with $P(AD_{alt} \ge 2) > 95 \%$
are classified as **expecting RNA support**.

> **NOTE**
>
> The assumptions of no allele-specific or cell-specific heterogeneity
> in expression is obviously not safe. We do this nonetheless, since any
> inaccuracy of our theoretical RNA support calculation is not
> problematic so long as the residuals to our model are reasonably
> consistent across samples.

**Intuitive explanation**

Assuming no allele-specific expression and uniform expression across all
cells, a DNA variant with a <u>50% variant allele frequency (VAF)</u>
should be present in approximately 50% of the RNA reads at that locus.
Here’s how RNA depth influences the probability of detecting such a
variant in the RNA, assuming we need at least 2 read support to call it:

| RNA Depth | Likelihood of Observing RNA support ($\ge 2$ **variant**-supportin**g reads)** |
|----|----|
| 2 | 25% |
| 3 | 50% |
| 10 | 99% |
| High Depth (e.g., \>= 100 reads) | Approaching certainty |

By filtering variants based on these probabilities, we can identify
those very likely to have RNA support. If these high-probability
variants show low RNA support, it raises concerns about their
authenticity (or RNA sample purity)

### Note on Assumptions

We acknowledge that assuming no allele-specific or cell-specific
heterogeneity in expression is an oversimplification. However,
inaccuracies in our theoretical RNA support calculations are acceptable
as long as average deviation from expected degrees of RNA support is
reasonably consistent across samples. It is the distribution of
$P_{supp}$, not its absolute values which are used to identify outliers
indicative of potential artificial mutations.
