library(testthat)
library(omixVizR)
library(ggplot2)

test_that("plot_gwas_power executes without errors and returns correct structure with binary traits", {
  # Run the function without saving the plot
  suppressMessages({power_results <- plot_gwas_power(
        trait_type = "bt",
        n_cases = 4324,
        n_controls = 93945,
        maf_levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
        or_range = seq(1.01, 2.00, 0.01),
        save_plot = FALSE
    )
  })
  
  # 1. Check the returned object structure
  expect_type(power_results, "list")
  expect_named(power_results, c("plot", "power_data"))
  
  # 2. Check the plot object
  expect_s3_class(power_results$plot, "ggplot")
  
  # 3. Check the data object
  expect_s3_class(power_results$power_data, "data.table")
})


test_that("plot_gwas_power executes without errors and returns correct structure with quantitative traits", {
  # Run the function without saving the plot
  suppressMessages({power_results_qt <- plot_gwas_power(
        trait_type = "qt",
        sd_trait = 0.09365788681305078,
        N = 10000,
        maf_levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
        effect_size = seq(0.01, 0.10, 0.001),
        save_plot = FALSE
    )
  })
  
  # 1. Check the returned object structure
  expect_type(power_results_qt, "list")
  expect_named(power_results_qt, c("plot", "power_data"))
  
  # 2. Check the plot object
  expect_s3_class(power_results_qt$plot, "ggplot")
  
  # 3. Check the data object
  expect_s3_class(power_results_qt$power_data, "data.table")
})
