test_that("multiplication works", {

  # Works as expected
  expect_equal(
    probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 1, min_alt_supporting_reads = 1), 0.5
  )

  expect_equal(
    probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 2, min_alt_supporting_reads = 1), 0.75
  )

  expect_equal(
    probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 2, min_alt_supporting_reads = 2), 0.25
  )
  expect_equal(
    probability_of_rna_support(DNAvaf = 0.1, RNAdepth = 2, min_alt_supporting_reads = 2), 0.01
  )

  # Works when vectors are supplied
  expect_equal(
    probability_of_rna_support(
      DNAvaf = c(0.5, 0.5, 0.1),
      RNAdepth = c(1, 2, 2),
      min_alt_supporting_reads = c(1, 1, 2)
    ),
    c(0.5, 0.75, 0.01)
    )

  # Works when a mix of vectors and integers supplied
  expect_equal(
    probability_of_rna_support(
      DNAvaf = c(0.5, 0.5),
      RNAdepth = c(1, 2),
      min_alt_supporting_reads = 1
    ),
    c(0.5, 0.75)
  )

  # Returns zero if min alt supporting reads is greater than depth
  expect_equal(
    probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 1, min_alt_supporting_reads = 2), 0
  )

  # Errors when min_alt_supporting_reads = 0
  expect_error(
    probability_of_rna_support(DNAvaf = 0.5, RNAdepth = 1, min_alt_supporting_reads = 0),
    "min_alt_supporting_reads must be greater than zero"
  )
})
