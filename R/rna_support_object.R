#' Create RNA Support Result
#'
#' This function creates an object of class `rna_support_result` which holds
#' the results of the RNA support evaluation, including the total number of
#' variants evaluated, expected RNA support variants, the proportion of variants
#' with observed support (expected & unexpected).
#'
#' @param total_variants Integer. The total number of variants evaluated.
#' @param expected_variants Integer. The number of variants expected to have RNA support.
#' @param observed_support_proportion Numeric. The proportion of expected variants with RNA support observed.
#' @param unexpected_support_proportion Numeric. The proportion of unexpected variants that still had RNA support.
#'
#' @return A list of class `rna_support` containing the following components:
#' \describe{
#'   \item{\strong{total_variants_evaluated}}{Integer. The total number of variants evaluated for RNA support.}
#'   \item{\strong{expected_rna_support_variants}}{Integer. The number of variants expected to show RNA support based on their DNA characteristics (e.g., variant allele frequency, read depth).}
#'   \item{\strong{proportion_with_observed_support}}{Numeric. The proportion of variants expected to show RNA support where RNA support was actually observed. Sometimes called 'Psupp'}
#'   \item{\strong{proportion_unexpected_with_observed_support}}{Numeric. The proportion of variants not expected to show RNA support that nonetheless had RNA support observed.}
#' }
#'
#' @export
#'
#' @examples
#' create_rna_support_result(
#'   total_variants = 100,
#'   expected_variants = 80,
#'   observed_support_proportion = 0.100,
#'   unexpected_support_proportion = 0.10
#' )
create_rna_support_result <- function(total_variants, expected_variants, observed_support_proportion, unexpected_support_proportion){

  result <- list(
    total_variants_evaluated = total_variants, # Total number of variants evaluated
    expected_rna_support_variants = expected_variants, # Number of variants expected to have RNA support
    proportion_unexpected_with_observed_support = unexpected_support_proportion, # Proportion of unexpected variants with observed RNA support
    proportion_with_observed_support = observed_support_proportion # Proportion of expected variants with observed RNA support
  )

  class(result) <- "rna_support"

  return(result)
}


#' Print method for rna_support object
#'
#' Custom print method to display the properties of an rna_support object.
#'
#' @param x An object of class `rna_support`.
#' @param ... Additional arguments passed to or from other methods (not used here).
#'
#' @keywords internal
#' @export
print.rna_support <- function(x, ...) {
  if (!inherits(x, "rna_support")) {
    stop("The object is not of class 'rna_support'")
  }

  cat("RNA Support Evaluation Result:\n")
  cat("----------------------------\n")
  cat("Total Variants: ", x$total_variants_evaluated, "\n")
  cat("[Nexpected] Variants Expected to have RNA Support: ", x$expected_rna_support_variants, "\n")
  cat(
    sprintf("[Psupp] Observed rate of support: %.1f%% (%d/%d)",
      x$proportion_with_observed_support * 100,
      round(x$proportion_with_observed_support * x$expected_rna_support_variants, digits = 0),
      x$expected_rna_support_variants
    ),
  "\n"
  )
  #cat("Observed Support (Unexpected): ", sprintf("%.1f%%", x$proportion_unexpected_with_observed_support * 100), "\n")
  cat("----------------------------\n")
}


#' Convert rna_support Object to Data Frame
#'
#' This function converts an object of class `rna_support` to a data frame.
#'
#' @param x An object of class `rna_support`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data frame with one row and columns corresponding to the components
#' of the `rna_support` object.
#'
#' @export
as.data.frame.rna_support <- function(x, ...) {
  if (!inherits(x, "rna_support")) {
    stop("The input object is not of class 'rna_support'.")
  }

  as.data.frame(x=unclass(x), stringsAsFactors = FALSE)
}
