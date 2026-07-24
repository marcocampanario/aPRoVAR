##### 08_article_figures.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Produce all figures from aPRoVAR 1.0 (Campanário and Janke et al, 2026)

# After running all the pipeline, and having all objects loaded.

#### Figure 1a ####
# Allele frequency spectrum of APROVAR-1010-WES

figure1a_data <- readxl::read_excel(file.path(TABLES_DIR, "07_allele_frequency_spectrum.xlsx"), sheet = 2)

old_levels <- if (is.factor(origin_spectrum$frequency_bin)) {
  levels(origin_spectrum$frequency_bin)
} else {
  unique(as.character(origin_spectrum$frequency_bin))
}

collapsed_levels <- unique(
  ifelse(
    tolower(old_levels) %in%
      c("singleton", "singletons", "doubleton", "doubletons"),
    "Singletons + Doubletons",
    old_levels
  )
)

origin_spectrum_collapsed <- origin_spectrum |>
  dplyr::mutate(
    frequency_bin = dplyr::if_else(
      tolower(as.character(frequency_bin)) %in%
        c("singleton", "singletons", "doubleton", "doubletons"),
      "Singletons + Doubletons",
      as.character(frequency_bin)
    )
  ) |>
  dplyr::group_by(frequency_bin, Origin) |>
  dplyr::summarise(
    n_variants = sum(n_variants),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    frequency_bin = factor(
      frequency_bin,
      levels = collapsed_levels
    )
  )

figure1a_plot <- ggplot2::ggplot(
  origin_spectrum_collapsed,
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
    labels = scales::label_number(big.mark = ","),
    breaks = c(0, 10, 1000, 10000, 100000, 300000)
  ) +
  ggplot2::labs(
    x = "aPRoVAR alternative allele frequency",
    y = "Number of variants (log scale)",
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
figure1a_plot

save_plot_formats(
  figure1a_plot,
  "Figure1a",
  width = 11,
  height = 5.5
)

#### Figure 1b ####
# Cumulative curve of variants per individual of APROVAR-1010-WES

figure1b_data <- readxl::read_excel(file.path(TABLES_DIR, "03_variant_discovery_and_frequency_counts.xlsx"), sheet = 1)

figure1b_plot_data <- figure1b_data |>
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

figure1b_plot <- ggplot2::ggplot(
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
    breaks = unique(c(scales::breaks_pretty()(0:number_of_samples), 1000)),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    breaks = c(0, 100000, 200000, 300000, 400000, 500000, 550000),
    expand = ggplot2::expansion(mult = c(0, 0.05))
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
    axis.text.x = ggplot2::element_text(color = "black", angle = 45),
    axis.text.y = ggplot2::element_text(color = "black"),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

figure1b_plot

save_plot_formats(
  figure1b_plot,
  "Figure1b",
  width = 9,
  height = 4
)

#### Figure 1c ####
# Barplot of rare and polymorphic genetic variants of aPRoVAR 1.0

figure1c_data <- readxl::read_excel(
  file.path(
    TABLES_DIR,
    "03_variant_discovery_and_frequency_counts.xlsx"
  ),
  sheet = 3
) |>
  dplyr::transmute(
    status = factor(
      dplyr::recode(
        variant_origin,
        known = "Known",
        absent = "Absent",
        novel = "Novel"
      ),
      levels = c("Known", "Absent", "Novel")
    ),
    classe = factor(
      dplyr::recode(
        broad_frequency_class,
        rare = "Rare",
        polymorphic = "Polymorphic"
      ),
      levels = c("Rare", "Polymorphic")
    ),
    n = n_variants,
    fill_key = paste(status, classe, sep = "_")
  ) |>
  dplyr::arrange(status, classe)

figure1c_pallette <- c("Known_Polymorphic"     = "#0072B2",
  "Known_Rare"     = "#0072B2",
  "Absent_Polymorphic"    = "goldenrod",
  "Absent_Rare"    = "goldenrod",
  "Novel_Polymorphic"     = "#009E73",
  "Novel_Rare"     = "#009E73"
)

figure1c_plot_scale_real <- ggplot2::ggplot(
  figure1c_data,
  ggplot2::aes(
    x = classe,
    y = n,
    fill = fill_key
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::label_number(big.mark = ",")(n)
    ),
    vjust = -0.25,
    size = 3.2
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(status),
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  ggplot2::scale_fill_manual(
    values = figure1c_pallette,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    breaks = c(0, 100000, 200000, 300000, 400000),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of variants (log scale)",
    fill = NULL,
    title = "Distribution of rare and polymorphic variants",
    subtitle = "Stratified by Known, Absent, and Novel categories"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11)
  )

figure1c_plot_scale_log <- ggplot2::ggplot(
  figure1c_data,
  ggplot2::aes(
    x = classe,
    y = n,
    fill = fill_key
  )
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::label_number(big.mark = ",")(n)
    ),
    vjust = -0.25,
    size = 3.2
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(status),
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  ggplot2::scale_fill_manual(
    values = figure1c_pallette,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number(big.mark = ","),
    breaks = c(0, 10, 100, 10000, 100000, 400000),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of variants (log scale)",
    fill = NULL,
    title = "Distribution of rare and polymorphic variants",
    subtitle = "Stratified by Known, Absent, and Novel categories"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11)
  )

save_plot_formats(
  figure1c_plot_scale_real,
  "Figure1c_scale_real",
  width = 9,
  height = 6
)

save_plot_formats(
  figure1c_plot_scale_log,
  "Figure1c_scale_log",
  width = 9,
  height = 6
)

#### Figure 1d ####
# Donut plots of rare variant composition of aPRoVAR 1.0

figure1d_data <- readxl::read_excel(file.path(TABLES_DIR, "03_variant_discovery_and_frequency_counts.xlsx"), sheet = 4)

figure1d_plot_data <- figure1d_data |>
  dplyr::mutate(
    Origin = format_origin(variant_origin),
    rare_component = factor(
      rare_component,
      levels = c("Singleton", "Doubleton", "AC > 2 and AF < 1%")
    )
  )

rare_components <- c(
  "Singleton",
  "Doubleton",
  "AC > 2 and AF < 1%"
)

figure1d_known_plot <- make_rare_donut(
  figure1d_plot_data,
  "Known"
)

figure1d_absent_plot <- make_rare_donut(
  figure1d_plot_data,
  "Absent"
)

figure1d_novel_plot <- make_rare_donut(
  figure1d_plot_data,
  "Novel"
)

save_plot_formats(
  figure1d_known_plot,
  "Figure1d_known",
  width = 5,
  height = 5
)

save_plot_formats(
  figure1d_absent_plot,
  "Figure1d_absent",
  width = 5,
  height = 5
)

save_plot_formats(
  figure1d_novel_plot,
  "Figure1d_novel",
  width = 5,
  height = 5
)

#### Figure 3 ####
# Word clouds of genes harboring P/LP variants in aPRoVAR 1.0

figure1d_data <- readxl::read_excel(file.path(TABLES_DIR, "06_wordcloud_data.xlsx"), sheet = 1)

wordcloud_plots <- list(
  Known = make_gene_wordcloud(wordcloud_data, "Known"),
  Absent = make_gene_wordcloud(wordcloud_data, "Absent"),
  Novel = make_gene_wordcloud(wordcloud_data, "Novel")
)

for (origin_name in names(wordcloud_plots)) {
  current_plot <- wordcloud_plots[[origin_name]]
  
  if (!is.null(current_plot)) {
    save_plot_formats(
      current_plot,
      paste0(
        "Figure3_",
        tolower(origin_name)
      ),
      width = 10,
      height = 10
    )
  }
}
