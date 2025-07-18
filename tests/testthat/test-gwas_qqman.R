library(testthat)
library(omixVizR)
library(ggplot2)

test_that("Sample data file exists", {
  sample_file <- system.file("extdata", "sample_gwas.assoc.linear", package = "omixVizR")
  expect_true(file.exists(sample_file), "Sample GWAS file should be present in inst/extdata.")
})

test_that("Function executes and returns a list of ggplot objects", {
  sample_file <- system.file("extdata", "sample_gwas.assoc.linear", package = "omixVizR")
  skip_if_not(file.exists(sample_file), "Sample file not found, skipping further tests.")
  
  plots <- plot_qqman(
    plink_assoc_file = sample_file,
    pheno_name = "SamplePheno",
    save_plot = FALSE
  )
  
  expect_type(plots, "list")
  expect_named(plots, c("manhattan_plot", "qq_plot"))
  expect_s3_class(plots$manhattan_plot, "ggplot")
  expect_s3_class(plots$qq_plot, "ggplot")
})
