#' @title Plot GWAS QQ and Manhattan Plots
#' @description Create GWAS QQ & Manhattan Plots.
#' @param plink_assoc_file Path to the PLINK association file.
#' @param pheno_name Phenotype name.
#' @param maf_filter Minor allele frequency filter, Default: NULL
#' @param output_graphics Output graphics format, Default: 'png'
#' @param save_plot Logical, whether to save plots to files. If FALSE, plots are only displayed. Default: TRUE
#' @param lambda1_qq_pos A numeric vector of length 2 specifying the `c(hjust, vjust)` for the lambda text in the QQ plot. Default: `c(2.1, -5.5)`.
#' @param lambda2_qq_pos A numeric vector of length 2 specifying the `c(hjust, vjust)` for the SNP count (N) text in the QQ plot. Default: `c(1.565, -4.0)`.
#' @return A list containing the ggplot objects for the Manhattan and QQ plots.
#' @author Zhen Lu <luzh29@mail2.sysu.edu.cn>
#' @author Yanhong Liu <liuyh275@mail2.sysu.edu.cn>
#' @author Siyang Liu <liusy99@mail.sysu.edu.cn>
#' @details
#' This function reads a PLINK association file and generates Manhattan and QQ plots for the GWAS results.
#' @examples
#' \donttest{
#'   sample_file <- system.file("extdata", "sample_gwas.assoc.linear", package = "omixVizR")
#'   
#'   # Check if the file exists before running the example
#'   if (file.exists(sample_file)) {
#'     # Run the function with the sample data
#'     plots <- plot_qqman(
#'       plink_assoc_file = sample_file,
#'       pheno_name = "SamplePheno",
#'       save_plot = FALSE
#'     )
#'     # You can then access the plots like this:
#'     # print(plots$manhattan_plot)
#'     # print(plots$qq_plot)
#'   } else {
#'     message("Sample file not found, skipping example.")
#'   }
#' }
#' @seealso
#'  \href{https://CRAN.R-project.org/package=lulab.utils}{lulab.utils}
#' @rdname plot_qqman
#' @export
#' @importFrom data.table fread .N := shift fifelse
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom magrittr %>% %<>%
#' @importFrom showtext showtext_auto
#' @importFrom ggtext element_markdown
#' @importFrom ggrepel geom_text_repel
#' @importFrom purrr map
#' @importFrom ggbreak scale_x_break
#' @importFrom grid unit
#' @importFrom scales pretty_breaks
#' @importFrom sysfonts font_families font_add
#' @importFrom graphics hist
#' @importFrom stats median qbeta qnorm
#' @section Font Information:
#' The MetroSans font included in this package is sourced from
#' \url{https://fontshub.pro/font/metro-sans-download#google_vignette}.
#' It is intended for academic research and non-commercial use only. For commercial use, please contact the font copyright holder.
#'
#' The font files are included in the package's inst/extdata directory and are automatically loaded for plotting.
plot_qqman = function(
  plink_assoc_file,
  pheno_name,
  maf_filter = NULL,
  output_graphics = "png",
  save_plot = TRUE,
  lambda1_qq_pos = c(2.1, -5.5),
  lambda2_qq_pos = c(1.565, -4.0)
) {
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
  gwasresults = data.table::fread(plink_assoc_file, header = TRUE)

  required_cols <- c("CHR", "BP", "P")
  if (!all(required_cols %in% colnames(gwasresults))) {
    stop("The input file must contain columns named 'CHR', 'BP', and 'P'.")
  }
  gwasresults[, BP := as.numeric(BP)]
  gwasresults[, P := as.numeric(P)]

  dat1 <- gwasresults[CHR != "NA" & !is.na(CHR), ] %>%
    .[, CHR := as.character(CHR)] %>%
    .[CHR == "X", CHR := "23"] %>%
    .[CHR == "Y", CHR := "24"] %>%
    .[, CHR := as.numeric(CHR)] %>%
    data.table::setorder(CHR, BP) %>%
    .[P != 0, ]

  fig_ylim <- ifelse(max(-log10(dat1$P), na.rm = TRUE) < 8, 12, max(-log10(dat1$P), na.rm = TRUE) + 4)

  if (data.table::uniqueN(dat1$CHR) == 1) {
    plotData <- dat1[, `:=`(
      maxBP = max(BP),
      halfBP = max(BP) / 2,
      xlabBP = max(BP) / 2,
      xaxis = BP
    )]
    dat3 <- plotData
  } else {
    dat2 <- dat1[, .(min = min(BP) - 1e7, max = max(BP) + 1e7), by = CHR] %>%
      data.table::melt(id.vars = "CHR", value.name = "BP") %>%
      .[, .(CHR, BP)] %>%
      .[, P := 10] %>%
      {data.table::rbindlist(list(dat1, .), fill = TRUE)} %>%
      data.table::setorder(CHR, BP)
    dat3 <- dat2[, .(maxBP = max(BP)), by = CHR] %>%
      .[, `:=`(
        accum = cumsum(maxBP),
        CHR = CHR + 1,
        halfBP = maxBP / 2
      )] %>%
      .[, .SD, .SDcols = c("CHR", "accum", "halfBP")] %>%
      {data.table::rbindlist(list(data.table::data.table(CHR = 1, accum = 0), .), fill = TRUE)} %>%
      data.table::setorder(CHR) %>%
      .[, xlabBP := data.table::shift(halfBP, 1, type = "lead") + accum] %>%
      .[CHR != max(dat2$CHR) + 1, ]
    plotData <- dat3[dat2, on = "CHR"] %>%
      .[, `:=`(
        xaxis = BP + accum,
        accum = NULL
      )]
  }

  # Manhattan plot
  plotData <- plotData[, highlight := data.table::fifelse(P <= 5e-8, "yes", "no", "no")][
    !is.na(SNP) & !is.na(P) & is.finite(-log10(P)),
  ]

  find_empty_y <- function(y, binwidth = 0.25, min_empty_bins = 3) {
    rng <- range(y, na.rm = TRUE)
    brks <- seq(floor(rng[1]), ceiling(rng[2]), by = binwidth)
    if (length(brks) < 3) return(NULL)
    h <- hist(y, breaks = brks, plot = FALSE)
    is_zero <- h$counts == 0
    r <- rle(is_zero)
    gap_id <- which(r$values & r$lengths >= min_empty_bins)
    if (!length(gap_id)) return(NULL)
    gaps <- purrr::map(gap_id, function(i) {
      start_bin <- sum(r$lengths[seq_len(i - 1)]) + 1
      end_bin <- start_bin + r$lengths[i] - 1
      c(brks[start_bin], brks[end_bin + 1])
    })
    unlist(gaps, use.names = FALSE)
  }
  gap_vec <- find_empty_y(-log10(plotData$P), binwidth = 30, min_empty_bins = 3)

  plot_manhattan <- ggplot2::ggplot(plotData) +
    ggplot2::geom_point(ggplot2::aes(x = xaxis, y = -log10(P), color = as.factor(CHR)), size = 0.2) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), color = "#DC0000FF", linetype = 2) +
    ggplot2::scale_x_continuous(breaks = dat3$xlabBP, labels = dat3$CHR, expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, fig_ylim),
      breaks = function(limits) {
        pretty_breaks <- scales::pretty_breaks()(limits)
        pretty_breaks[pretty_breaks < fig_ylim]
      }
    ) +
    ggplot2::geom_point(data = plotData[highlight == "yes", ], ggplot2::aes(x = xaxis, y = -log10(P)), color = "#D55E00", size = 0.5) +
    ggrepel::geom_text_repel(
      data = plotData[highlight == "yes", ],
      ggplot2::aes(x = xaxis, y = -log10(P), label = SNP), color = "#D55E00", size = 6, nudge_y = 1.0, fontface = "italic",
      max.overlaps = 20,
      segment.color = "#D55E00" # #009E73
    ) +
    ggplot2::scale_color_manual(values = rep(c("#3B5488", "#53BBD5"), 12)) +
    ggplot2::theme_classic(base_size = 25) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "none",
      text = ggplot2::element_text(family = "MetroSans"),
      strip.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 13, face = "bold"),
      axis.title = ggtext::element_markdown(size = 16.5, face = "bold"),
      axis.ticks.length = grid::unit(0.12, "cm"),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      axis.ticks = ggplot2::element_line(linewidth = 0.8),
      axis.text.y.right = ggplot2::element_blank(),
      axis.ticks.y.right = ggplot2::element_blank(),
      axis.line.y.right = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Chromosome",
      y = "-log<sub>10</sub>(<i>P</i>)"
    )

  if (!is.null(gap_vec)) {
    pb <- ggplot2::ggplot_build(plot_manhattan)
    y_breaks <- pb$layout$panel_params[[1]]$y$get_breaks()
    gap_mat <- matrix(gap_vec, ncol = 2, byrow = TRUE)
    inside_gap <- function(y) any(y > gap_mat[, 1] & y < gap_mat[, 2])
    y_breaks <- y_breaks[!purrr::map_lgl(y_breaks, inside_gap)]
    plot_manhattan <- plot_manhattan +
      ggbreak::scale_y_break(
        breaks = gap_vec, scales = 0.7,
        expand = c(0, 0), ticklabels = y_breaks
      )
  }

  if (save_plot) {
    output_file <- file.path(output_dir, paste0(pheno_name, "_Manhattan_plot.", output_graphics))
    ggplot2::ggsave(
      filename = output_file,
      plot = plot_manhattan,
      width = 16, height = 4, dpi = 600
    )
    message("Manhattan plot saved to: ", output_file)
  }
  # QQ plot
  qqPlotData <- plotData[, .(P)]
  qqPlotData[, P_v2 := P]
  qqPlotData[, y := -log10(P_v2)]
  data.table::setorder(qqPlotData, -y)
  qqPlotData[, a := 1:.N]
  qqPlotData[, x := -log10(a / .N)]
  qqPlotData[, upper := -log10(qbeta(0.025, 1:.N, (.N):1))]
  qqPlotData[, lower := -log10(qbeta(0.975, 1:.N, (.N):1))]
  lambda <- round(median(qnorm(qqPlotData$P / 2)^2, na.rm = TRUE) / 0.454, 3)
  N <- nrow(qqPlotData)
  lambdaData <- data.table::data.table(
    label1 = paste0("lambda[GC] == ", lambda),
    label2 = paste0("N[italic(SNP)] == ", N)
  )

  plot_qq <- ggplot2::ggplot(qqPlotData) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y), size = 0.2) +
    ggplot2::geom_text(
      data = lambdaData,
      ggplot2::aes(x = Inf, y = -Inf, label = label1),
      hjust = lambda1_qq_pos[1], vjust = lambda1_qq_pos[2], size = 8,
      fontface = "bold",
      parse = TRUE
    ) +
    ggplot2::geom_text(
      data = lambdaData,
      ggplot2::aes(x = Inf, y = -Inf, label = label2),
      hjust = lambda2_qq_pos[1], vjust = lambda2_qq_pos[2], size = 8,
      fontface = "bold",
      parse = TRUE
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      fill = "gray20", alpha = 0.3
    ) +
    ggplot2::theme_bw(base_size = 25) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      text = ggplot2::element_text(family = "MetroSans"),
      strip.text = ggplot2::element_blank(),
      aspect.ratio = 1,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 13, face = "bold"),
      axis.title = ggtext::element_markdown(size = 16.5, face = "bold"),
      axis.ticks.length = grid::unit(0.09, "cm"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black")
    ) +
    ggplot2::labs(
      x = "Expected -log<sub>10</sub>(<i>P</i>)",
      y = "Observed -log<sub>10</sub>(<i>P</i>)"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, fig_ylim)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, fig_ylim)) +
    ggplot2::coord_fixed()

  if (save_plot) {
    output_file <- file.path(output_dir, paste0(pheno_name, "_QQ_plot.", output_graphics))
    ggplot2::ggsave(
      filename = output_file,
      plot = plot_qq,
      width = 8, height = 8, dpi = 600
    )
    message("QQ plot saved to: ", output_file)
  }

  return(invisible(list(manhattan_plot = plot_manhattan, qq_plot = plot_qq)))
}