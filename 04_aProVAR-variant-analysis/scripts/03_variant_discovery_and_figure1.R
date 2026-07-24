##### 03_variant_discovery_and_figure1.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Quantify cumulative variant discovery across randomly ordered individuals and
# generate the frequency-class summaries used in aPRoVAR Figure 1.

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

require_packages(c("dplyr", "ggplot2", "purrr", "scales", "tidyr", "writexl"))
create_output_directories()
assert_file_exists(ANALYSIS_READY_RDS, "Analysis-ready checkpoint")

message("Reading analysis-ready variants and the presence matrix...")
analysis_ready <- read_rds_parallel(ANALYSIS_READY_RDS)
variant_table <- analysis_ready$analysis_table
presence_matrix <- analysis_ready$presence_matrix

if (nrow(variant_table) != nrow(presence_matrix)) {
  stop(
    "The analysis table and presence matrix do not contain the same number of sites.",
    call. = FALSE
  )
}

set.seed(RANDOM_SEED)
sample_order <- sample(colnames(presence_matrix))
number_of_samples <- length(sample_order)

# Determine the first randomly ordered individual in whom each site appears.
# Processing row blocks avoids creating a second full 500k-by-1,010 matrix.
message("Calculating the first observation of each variant site...")
block_size <- 25000L
first_seen <- rep(NA_integer_, nrow(presence_matrix))

for (block_start in seq.int(1L, nrow(presence_matrix), by = block_size)) {
  block_end <- min(block_start + block_size - 1L, nrow(presence_matrix))
  rows <- block_start:block_end
  ordered_block <- presence_matrix[rows, sample_order, drop = FALSE]
  has_variant <- rowSums(ordered_block) > 0L
  first_in_block <- max.col(ordered_block, ties.method = "first")
  first_in_block[!has_variant] <- NA_integer_
  first_seen[rows] <- first_in_block
}

variant_table <- variant_table |>
  dplyr::mutate(
    first_seen = first_seen,
    broad_frequency_class = dplyr::case_when(
      variant_freq_class %in% c(
        "singleton",
        "doubleton",
        "doubleton_hom",
        "rare"
      ) ~ "rare",
      variant_freq_class == "polymorphic" ~ "polymorphic",
      TRUE ~ NA_character_
    ),
    cumulative_frequency_class = dplyr::case_when(
      variant_freq_class %in% c(
        "singleton",
        "doubleton",
        "doubleton_hom"
      ) ~ "singletons",
      variant_freq_class == "rare" ~ "rare",
      variant_freq_class == "polymorphic" ~ "polymorphic",
      TRUE ~ NA_character_
    )
  )

cumulative_count <- function(mask) {
  discovery_positions <- variant_table$first_seen[mask]
  discovery_positions <- discovery_positions[!is.na(discovery_positions)]
  c(0L, cumsum(tabulate(discovery_positions, nbins = number_of_samples)))
}

cumulative_results <- data.frame(
  n_individuals = 0:number_of_samples,
  total = cumulative_count(rep(TRUE, nrow(variant_table))),
  known = cumulative_count(variant_table$variant_origin == "known"),
  absent = cumulative_count(variant_table$variant_origin == "absent"),
  novel = cumulative_count(variant_table$variant_origin == "novel"),
  known_singletons = cumulative_count(
    variant_table$variant_origin == "known" &
      variant_table$cumulative_frequency_class == "singletons"
  ),
  known_rare = cumulative_count(
    variant_table$variant_origin == "known" &
      variant_table$cumulative_frequency_class == "rare"
  ),
  known_polymorphic = cumulative_count(
    variant_table$variant_origin == "known" &
      variant_table$cumulative_frequency_class == "polymorphic"
  ),
  absent_singletons = cumulative_count(
    variant_table$variant_origin == "absent" &
      variant_table$cumulative_frequency_class == "singletons"
  ),
  absent_rare = cumulative_count(
    variant_table$variant_origin == "absent" &
      variant_table$cumulative_frequency_class == "rare"
  ),
  absent_polymorphic = cumulative_count(
    variant_table$variant_origin == "absent" &
      variant_table$cumulative_frequency_class == "polymorphic"
  ),
  novel_singletons = cumulative_count(
    variant_table$variant_origin == "novel" &
      variant_table$cumulative_frequency_class == "singletons"
  ),
  novel_rare = cumulative_count(
    variant_table$variant_origin == "novel" &
      variant_table$cumulative_frequency_class == "rare"
  ),
  novel_polymorphic = cumulative_count(
    variant_table$variant_origin == "novel" &
      variant_table$cumulative_frequency_class == "polymorphic"
  ),
  check.names = FALSE
)

frequency_counts <- variant_table |>
  dplyr::filter(!is.na(variant_freq_class)) |>
  dplyr::count(
    variant_origin,
    variant_freq_class,
    name = "n_variants"
  ) |>
  tidyr::complete(
    variant_origin = c("known", "absent", "novel"),
    variant_freq_class = c(
      "singleton",
      "doubleton",
      "doubleton_hom",
      "rare",
      "polymorphic"
    ),
    fill = list(n_variants = 0L)
  )

broad_frequency_counts <- variant_table |>
  dplyr::filter(!is.na(broad_frequency_class)) |>
  dplyr::count(
    variant_origin,
    broad_frequency_class,
    name = "n_variants"
  ) |>
  tidyr::complete(
    variant_origin = c("known", "absent", "novel"),
    broad_frequency_class = c("rare", "polymorphic"),
    fill = list(n_variants = 0L)
  )

rare_composition <- variant_table |>
  dplyr::filter(
    variant_freq_class %in% c(
      "singleton",
      "doubleton",
      "doubleton_hom",
      "rare"
    )
  ) |>
  dplyr::mutate(
    rare_component = dplyr::case_when(
      variant_freq_class == "singleton" ~ "Singleton",
      variant_freq_class %in% c("doubleton", "doubleton_hom") ~ "Doubleton",
      variant_freq_class == "rare" ~ "AC > 2 and AF < 1%"
    )
  ) |>
  dplyr::count(variant_origin, rare_component, name = "n_variants") |>
  dplyr::group_by(variant_origin) |>
  dplyr::mutate(
    proportion = n_variants / sum(n_variants),
    label = paste0(
      scales::comma(n_variants),
      "\n(",
      scales::percent(proportion, accuracy = 0.1),
      ")"
    )
  ) |>
  dplyr::ungroup()

writexl::write_xlsx(
  list(
    cumulative_discovery = cumulative_results,
    detailed_frequency_counts = frequency_counts,
    rare_vs_polymorphic = broad_frequency_counts,
    rare_composition = rare_composition
  ),
  file.path(TABLES_DIR, "03_variant_discovery_and_frequency_counts.xlsx")
)

# Cumulative variant discovery ------------------------------------------------

cumulative_plot_data <- cumulative_results |>
  tidyr::pivot_longer(
    cols = c(total, known, absent, novel),
    names_to = "variant_origin",
    values_to = "n_variants"
  ) |>
  dplyr::mutate(
    variant_origin = factor(
      variant_origin,
      levels = c("total", "known", "absent", "novel"),
      labels = c("Total", "Known", "Absent", "Novel")
    )
  )

cumulative_colors <- c("Total" = "#D55E00", ORIGIN_COLORS)

cumulative_plot <- ggplot2::ggplot(
  cumulative_plot_data,
  ggplot2::aes(
    x = n_individuals,
    y = n_variants,
    color = variant_origin
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_color_manual(values = cumulative_colors) +
  ggplot2::scale_x_continuous(
    breaks = unique(c(scales::breaks_pretty()(0:number_of_samples), number_of_samples)),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  ) +
  ggplot2::labs(
    x = "Cumulative number of individuals",
    y = "Cumulative number of variants",
    color = NULL,
    title = "Cumulative discovery of genetic variants across aPRoVAR"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  cumulative_plot,
  "Figure1_cumulative_variant_discovery",
  width = 9,
  height = 4.5
)

# Rare versus polymorphic sites -----------------------------------------------

broad_plot_data <- broad_frequency_counts |>
  dplyr::mutate(
    Origin = format_origin(variant_origin),
    `Frequency class` = factor(
      broad_frequency_class,
      levels = c("rare", "polymorphic"),
      labels = c("Rare (<1%)", "Polymorphic (≥1%)")
    )
  )

broad_plot <- ggplot2::ggplot(
  broad_plot_data,
  ggplot2::aes(
    x = `Frequency class`,
    y = n_variants,
    fill = Origin
  )
) +
  ggplot2::geom_col(position = "dodge", width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(n_variants)),
    position = ggplot2::position_dodge(width = 0.72),
    vjust = -0.25,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = ORIGIN_COLORS) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of variants",
    fill = NULL,
    title = "Rare and polymorphic variants by population-database status"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  broad_plot,
  "Figure1_rare_vs_polymorphic_variants",
  width = 9,
  height = 4.5
)

# Composition of rare sites ---------------------------------------------------

rare_plot_data <- rare_composition |>
  dplyr::mutate(
    Origin = format_origin(variant_origin),
    rare_component = factor(
      rare_component,
      levels = c("Singleton", "Doubleton", "AC > 2 and AF < 1%")
    )
  )

rare_plot <- ggplot2::ggplot(
  rare_plot_data,
  ggplot2::aes(x = 2, y = proportion, fill = rare_component)
) +
  ggplot2::geom_col(color = "white", width = 1) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3
  ) +
  ggplot2::coord_polar(theta = "y") +
  ggplot2::xlim(0.5, 2.5) +
  ggplot2::facet_wrap(~ Origin, nrow = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      "Singleton" = "#56B4E9",
      "Doubleton" = "#CC79A7",
      "AC > 2 and AF < 1%" = "#0072B2"
    )
  ) +
  ggplot2::labs(
    fill = NULL,
    title = "Composition of rare variants"
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  rare_plot,
  "Figure1_rare_variant_composition",
  width = 10,
  height = 4.5
)

discovery_results <- list(
  sample_order = sample_order,
  cumulative_results = cumulative_results,
  frequency_counts = frequency_counts,
  broad_frequency_counts = broad_frequency_counts,
  rare_composition = rare_composition
)

write_rds_parallel(discovery_results, DISCOVERY_RESULTS_RDS)
message("Variant-discovery analysis complete.")

