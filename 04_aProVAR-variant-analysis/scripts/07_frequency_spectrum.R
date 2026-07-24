##### 07_frequency_spectrum.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Summarize the aPRoVAR allele-frequency spectrum using explicit singleton and
# doubleton classes plus numeric frequency bins, overall and by variant origin.

project_dir <- Sys.getenv("APROVAR_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  script_path <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
  )
  project_dir <- dirname(dirname(normalizePath(script_path, mustWork = FALSE)))
  Sys.setenv(APROVAR_PROJECT_DIR = project_dir)
}

source(file.path(project_dir, "R", "00_config.R"))
source(file.path(project_dir, "R", "01_functions.R"))

require_packages(c("dplyr", "ggplot2", "scales", "tidyr", "writexl"))
create_output_directories()
assert_file_exists(ANALYSIS_READY_RDS, "Analysis-ready checkpoint")

analysis_ready <- read_rds_parallel(ANALYSIS_READY_RDS)
variant_table <- analysis_ready$analysis_table

frequency_breaks <- c(
  0,
  0.01,
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30,
  0.35,
  0.40,
  0.45,
  0.50
)

frequency_labels <- c(
  "[0, 0.01)",
  "[0.01, 0.05)",
  "[0.05, 0.10)",
  "[0.10, 0.15)",
  "[0.15, 0.20)",
  "[0.20, 0.25)",
  "[0.25, 0.30)",
  "[0.30, 0.35)",
  "[0.35, 0.40)",
  "[0.40, 0.45)",
  "[0.45, 0.50)"
)

spectrum_data <- variant_table |>
  dplyr::mutate(
    numeric_bin = cut(
      af_local,
      breaks = frequency_breaks,
      labels = frequency_labels,
      include.lowest = TRUE,
      right = FALSE
    ),
    frequency_bin = dplyr::case_when(
      variant_freq_class == "singleton" ~ "Singleton",
      variant_freq_class %in% c("doubleton", "doubleton_hom") ~ "Doubleton",
      TRUE ~ as.character(numeric_bin)
    )
  ) |>
  dplyr::filter(!is.na(frequency_bin))

frequency_bin_levels <- c("Singleton", "Doubleton", frequency_labels)

overall_spectrum <- spectrum_data |>
  dplyr::count(frequency_bin, name = "n_variants") |>
  tidyr::complete(
    frequency_bin = frequency_bin_levels,
    fill = list(n_variants = 0L)
  ) |>
  dplyr::mutate(
    frequency_bin = factor(frequency_bin, levels = frequency_bin_levels)
  ) |>
  dplyr::arrange(frequency_bin)

origin_spectrum <- spectrum_data |>
  dplyr::count(
    variant_origin,
    frequency_bin,
    name = "n_variants"
  ) |>
  tidyr::complete(
    variant_origin = c("known", "absent", "novel"),
    frequency_bin = frequency_bin_levels,
    fill = list(n_variants = 0L)
  ) |>
  dplyr::mutate(
    frequency_bin = factor(frequency_bin, levels = frequency_bin_levels),
    Origin = format_origin(variant_origin)
  ) |>
  dplyr::arrange(Origin, frequency_bin)

writexl::write_xlsx(
  list(
    overall = overall_spectrum,
    by_variant_origin = origin_spectrum |>
      dplyr::select(-Origin)
  ),
  file.path(TABLES_DIR, "07_allele_frequency_spectrum.xlsx")
)

overall_plot <- ggplot2::ggplot(
  overall_spectrum,
  ggplot2::aes(x = frequency_bin, y = n_variants)
) +
  ggplot2::geom_col(fill = "#0072B2", width = 0.78) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::labs(
    x = "aPRoVAR alternative allele frequency",
    y = "Number of variants (pseudo-log scale)",
    title = "Allele-frequency spectrum in aPRoVAR"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(color = "black"),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  overall_plot,
  "Allele_frequency_spectrum_overall",
  width = 10,
  height = 5
)

origin_plot <- ggplot2::ggplot(
  origin_spectrum,
  ggplot2::aes(
    x = frequency_bin,
    y = n_variants,
    fill = Origin
  )
) +
  ggplot2::geom_col(position = "dodge", width = 0.78) +
  ggplot2::scale_fill_manual(values = ORIGIN_COLORS) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::labs(
    x = "aPRoVAR alternative allele frequency",
    y = "Number of variants (pseudo-log scale)",
    fill = NULL,
    title = "Allele-frequency spectrum by population-database status"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(color = "black"),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  origin_plot,
  "Allele_frequency_spectrum_by_variant_origin",
  width = 11,
  height = 5.5
)

message("Allele-frequency spectrum analysis complete.")

