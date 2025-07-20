#' @title Plot GWAS Statistical Power
#' @description
#' Generates a statistical power analysis plot for GWAS studies.
#' Supports binary (case-control) traits over a range of odds ratios and minor allele frequencies,
#' and quantitative traits over a range of effect sizes and minor allele frequencies.
#' This function uses the 'genpwr' package for calculations and creates a highly customized ggplot.
#'
#' @param trait_type Character string specifying trait type: "bt" for binary (case-control) or "qt" for quantitative traits. Default: `"bt"`.
#' @param n_cases Number of cases in the study (required if `trait_type = "bt"`).
#' @param n_controls Number of controls in the study (required if `trait_type = "bt"`).
#' @param sd_trait Numeric, standard deviation of the quantitative trait (required if `trait_type = "qt"`).
#' @param N Numeric, total sample size for quantitative traits (required if `trait_type = "qt"`).
#' @param maf_levels A numeric vector of Minor Allele Frequencies (MAFs) to test.
#'   Default: `c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50)`.
#' @param or_range A numeric vector specifying the sequence of Odds Ratios (ORs) to test.
#'   Default: `seq(1.01, 2.00, 0.001)`. Used when `trait_type = "bt"`.
#' @param effect_size A numeric vector specifying the sequence of effect sizes (beta) to test for quantitative traits.
#'   Default: `seq(0.01, 0.30, 0.001)`. Used when `trait_type = "qt"`.
#' @param alpha The significance level (alpha) for the power calculation.
#'   Default: `5e-8`.
#' @param plot_title A string for the plot title. Can include newlines (`\\n`).
#'   Default: A title generated from case/control numbers.
#' @param save_plot Logical, whether to save the plot to a file. If `FALSE`, the
#'   plot object is only returned. Default: `TRUE`.
#' @param output_graphics The file format for saving the plot. Currently supports
#'   "png" and "pdf". Default: "png".
#' @param width The width of the saved plot in inches. Default: 17.
#' @param height The height of the saved plot in inches. Default: 9.
#' @param dpi The resolution of the saved plot in dots per inch. Default: 600.
#'
#' @return A list containing two elements:
#'   \item{plot}{The ggplot object for the power plot.}
#'   \item{power_data}{A data.table containing the full results from the power analysis.}
#'
#' @author Zhen Lu <luzh29@mail2.sysu.edu.cn>
#'
#' @details
#' This function automates the process of calculating and visualizing GWAS power for both
#' binary (case-control) and quantitative traits. For binary traits, it analyzes power across
#' odds ratios, while for quantitative traits, it analyzes power across effect sizes.
#' It highlights the minimum OR/effect size required to achieve 80% power for the lowest and
#' third-lowest MAF levels, adding dashed lines and color-coded labels for clarity.
#'
#' @examples
#' \donttest{
#'   # Binary trait example (case-control)
#'   power_results_bt <- plot_gwas_power(
#'     trait_type = "bt",
#'     n_cases = 4324,
#'     n_controls = 93945,
#'     save_plot = FALSE
#'   )
#'
#'   # Quantitative trait example
#'   power_results_qt <- plot_gwas_power(
#'     trait_type = "qt",
#'     sd_trait = 0.09365788681305078,
#'     N = 10000,
#'     maf_levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
#'     effect_size = seq(0.01, 0.10, 0.001),
#'     save_plot = FALSE
#'   )
#'
#'   # Access the ggplot object and data
#'   # print(power_results_bt$plot)
#'   # print(power_results_bt$power_data)
#' }
#'
#' @rdname plot_gwas_power
#' @export
#' @importFrom data.table as.data.table :=
#' @importFrom genpwr genpwr.calc
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual expand_limits
#'   scale_x_continuous scale_y_continuous labs theme_classic geom_hline
#'   geom_segment theme ggsave margin expansion element_text
#' @importFrom ggtext element_markdown
#' @importFrom ggsci pal_npg
#' @importFrom showtext showtext_auto showtext_opts
#' @importFrom sysfonts font_add
#' @importFrom scales label_number
#' @importFrom dplyr filter arrange case_when if_else
#' @importFrom magrittr %>% %<>%
#' @importFrom utils head
#' @importFrom stats setNames
#' @section Font Information:
#' The MetroSans font included in this package is sourced from
#' \url{https://fontshub.pro/font/metro-sans-download#google_vignette}.
#' It is intended for academic research and non-commercial use only. For commercial use, please contact the font copyright holder.
#'
#' The font files are included in the package's inst/extdata directory and are automatically loaded for plotting.
#'
plot_gwas_power <- function(
  trait_type = "bt",
  n_cases= NULL,
  n_controls= NULL,
  sd_trait= NULL,
  N= NULL,
  maf_levels = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.50),
  or_range = seq(1.01, 2.00, 0.001),
  effect_size= seq(0.01, 0.30, 0.001),
  alpha = 5e-8,
  plot_title = NULL,
  save_plot = TRUE,
  output_graphics = "png",
  width = 17,
  height = 9,
  dpi = 600
) {

  if (trait_type == "bt" & (is.null(n_cases) | is.null(n_controls))) {
    stop("For binary traits, please provide n_cases and n_controls.")
  } else if (trait_type == "qt" & (is.null(sd_trait) | is.null(N))) {
    stop("For quantitative traits, please provide sd_trait and N.")
  }

  # --- 1. Setup: Fonts and Colors ---
  font_dir <- system.file("extdata", package = "omixVizR")
  if (!"MetroSans" %in% sysfonts::font_families()) {
    sysfonts::font_add(
      family = "MetroSans",
      regular = file.path(font_dir, "MetroSans-Regular.ttf"),
      bold = file.path(font_dir, "MetroSans-Bold.ttf"),
      bolditalic = file.path(font_dir, "MetroSans-BoldItalic.ttf")
    )
  }
  showtext::showtext_auto(enable = TRUE)
  showtext::showtext_opts(dpi = 600)

  output_dir = tempdir()
  output_file <- file.path(output_dir, paste0("gwas_power_plot.", output_graphics))

  plot_colors <- ggsci::pal_npg("nrc")(length(maf_levels))

  # --- 2. Power Calculation ---
  if (trait_type == "bt") {
    total_n <- n_cases + n_controls
    case_rate <- n_cases / total_n

    power_data <- genpwr::genpwr.calc(
      calc = "power", model = "logistic", ge.interaction = NULL,
      N = total_n, Case.Rate = case_rate, k=NULL,
      MAF = maf_levels, OR = or_range,
      Alpha = alpha,
      True.Model = "Additive", Test.Model = "Additive"
    )
  } else if (trait_type == "qt") {
    power_data <- genpwr::genpwr.calc(
      calc = "power", model = "linear", ge.interaction = NULL,
      N = N, sd_y = sd_trait, k=NULL,
      MAF = maf_levels, ES= effect_size,
      Alpha = alpha,
      True.Model = "Additive", Test.Model = "Additive"
    )
  }

  # --- 3. Data Preparation for Plotting ---
  maf_labels <- paste("MAF =", format(maf_levels))
  power_data$MAF %<>% factor(levels = maf_levels, labels = maf_labels)

  # --- 4. Dynamic Calculation for Annotations ---
  # Find minimum effect size for 80% power at specific MAFs for annotations
  data.table::setDT(power_data)
  if (trait_type == "bt") {
    power_data[, effect_size_value := OR]
  } else if (trait_type == "qt") {
    power_data[, effect_size_value := ES]
  }
  min_or_maf1 <- power_data %>%
      dplyr::filter(`Power_at_Alpha_5e-08` >= 0.80) %>%
      dplyr::arrange(MAF, effect_size_value) %>%
      utils::head(1)
  min_or_maf3 <- power_data %>%
    dplyr::filter(`Power_at_Alpha_5e-08` >= 0.80, MAF == maf_labels[3]) %>%
    dplyr::arrange(effect_size_value) %>%
    utils::head(1)

  # --- 5. Create Plot ---
  # Set default title if not provided
  if (is.null(plot_title)) {
    if (trait_type == "bt") {
      plot_title <- sprintf(
        "Cases / Controls = %s / %s\nGWAS Significance Level = %s",
        format(n_cases, big.mark = ","),
        format(n_controls, big.mark = ","),
        format(alpha, scientific = TRUE)
      )
    } else if (trait_type == "qt") {
      plot_title <- sprintf(
        "Sample Size = %s\nGWAS Significance Level = %s",
        format(N, big.mark = ","),
        format(alpha, scientific = TRUE)
      )
    }
  }

  if (trait_type == "qt") {
    x_axis_label <- expression("Effect Size (" * beta * ")")
  } else {
    x_axis_label <- "Odds Ratio (OR)"
  }

  get_axis_step <- function(value) {
    thresholds <- c(0.5, 1, 2, 3)
    steps <- c(0.01, 0.1, 0.2, 0.5, 1.0)
    
    idx <- findInterval(value, thresholds)+1
    return(steps[idx])
  }
  detect_decimal_places <- function(x) {
    x_char <- format(x, scientific = FALSE)
    decimal_places <- sapply(x_char, function(num) {
      if (grepl("\\.", num)) {
        decimal_part <- sub(".*\\.", "", num)
        decimal_part <- sub("0+$", "", decimal_part)
        nchar(decimal_part)
      } else {
        0
      }
    })
    return(max(decimal_places))
  }

  # Define dynamic breaks and labels for the x-axis
  if (trait_type == "bt") {
    min_eff <- floor(min(as.numeric(or_range), na.rm = TRUE))
    max_eff <- ceiling(max(as.numeric(or_range), na.rm = TRUE))
    eff_range <- max(or_range) - min(or_range)
    by_eff= get_axis_step(eff_range)
    x_breaks <- sort(unique(c(seq(min_eff, max_eff, ifelse(by_eff==0.01, 0.1, by_eff)), min_or_maf3$effect_size_value, min_or_maf1$effect_size_value)))
    decimal_places= detect_decimal_places(max(or_range))
    multiplier <- 10^decimal_places
    max_limit=  ceiling(max(or_range) * multiplier) / multiplier
    xaxis_limit= c(min_eff, max_limit)
  } else if (trait_type == "qt") {
    min_eff <- floor(min(as.numeric(effect_size), na.rm = TRUE))
    max_eff <- ceiling(max(as.numeric(effect_size), na.rm = TRUE))
    eff_range <- max(effect_size) - min(effect_size)
    by_eff= get_axis_step(eff_range)
    x_breaks <- sort(unique(c(seq(min_eff, max_eff, ifelse(by_eff==0.01, 0.1, by_eff)), min_or_maf3$effect_size_value, min_or_maf1$effect_size_value)))
    decimal_places= detect_decimal_places(max(effect_size))
    multiplier <- 10^decimal_places
    max_limit=  ceiling(max(effect_size) * multiplier) / multiplier
    xaxis_limit= c(min_eff, max_limit)
  }
  
  x_labels <- function(x) {
    dplyr::case_when(
      x == min_or_maf3$effect_size_value ~ paste0("<span style='color:", plot_colors[3], "'>", sprintf("%.2f", x), "</span>"),
      x == min_or_maf1$effect_size_value ~ paste0("<span style='color:", plot_colors[1], "'>", sprintf("%.2f", x), "</span>"),
      TRUE ~ paste0("<span style='color:black'>", scales::label_number(accuracy = 0.01)(x), "</span>")
    )
  }

  p <- ggplot2::ggplot(power_data, ggplot2::aes(x = effect_size_value, y = `Power_at_Alpha_5e-08`)) +
    ggplot2::geom_line(ggplot2::aes(color = MAF), linewidth = 1.28) +
    ggplot2::scale_color_manual(
      values = stats::setNames(plot_colors, levels(power_data$MAF)),
      name = "Minor Allele Frequency (MAF)",
      guide = ggplot2::guide_legend(override.aes = list(linewidth = 2.8))
    ) +
    ggplot2::expand_limits(x = xaxis_limit, y = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1.0, 0.2),
      expand = ggplot2::expansion(mult = c(0, 0.08), add = c(0, 0)),
      labels = function(x) dplyr::if_else(x == 0, "0", sprintf("%.2f", x))
    ) +
    ggplot2::labs(
      x = x_axis_label,
      y = "GWAS Statistical Power",
      title = plot_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_hline(yintercept = 0.80, linewidth = 0.88, linetype = "dashed") +
    # Dashed lines for specific ORs
    ggplot2::geom_segment(
      data = min_or_maf3,
      ggplot2::aes(x = effect_size_value, xend = effect_size_value, y = 0, yend = 0.80),
      linewidth = 0.88, linetype = "dashed", color = plot_colors[3]
    ) +
    ggplot2::geom_segment(
      data = min_or_maf1,
      ggplot2::aes(x = effect_size_value, xend = effect_size_value, y = 0, yend = 0.80),
      linewidth = 0.88, linetype = "dashed", color = plot_colors[1]
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "MetroSans", size = 21),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 18, colour = "black"),
      axis.text.x = ggtext::element_markdown(size = 18),
      legend.position = "inside",
      legend.position.inside = c(0.88, 0.40),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 24, vjust = -0.8),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(t = 15, r = 10, b = 10, l = 10, unit = "pt")
    )
  if (trait_type == "qt"){
    p= p + ggplot2::theme(axis.title.x = ggplot2::element_text(family = "sans", face = "bold"))
  }

  # --- 6. Save Plot (if requested) ---
  if (save_plot) {
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width, height = height, dpi = dpi
    )
    message("GWAS power plot saved to: ", output_file)
  }

  # --- 7. Return Results ---
  return(invisible(list(plot = p, power_data = data.table::as.data.table(power_data))))
}
