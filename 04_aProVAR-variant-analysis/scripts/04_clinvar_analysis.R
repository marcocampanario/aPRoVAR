##### 04_clinvar_analysis.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Extract and consolidate ClinVar germline assertions, summarize their
# distribution across Known, Absent, and Novel variants, and generate gene-level
# tables and word clouds for likely pathogenic/pathogenic variants.

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

require_packages(c(
  "dplyr",
  "ggplot2",
  "purrr",
  "readr",
  "scales",
  "stringr",
  "tidyr"
))
create_output_directories()
assert_file_exists(ANALYSIS_READY_RDS, "Analysis-ready checkpoint")

message("Reading analysis-ready variants...")
analysis_ready <- read_rds_parallel(ANALYSIS_READY_RDS)
variant_table <- analysis_ready$analysis_table

message("Extracting ClinVar assertions, variant types, genes, and phenotypes...")
annotated_variants <- variant_table |>
  dplyr::mutate(
    clinvar_classification = purrr::map_chr(
      variants,
      extract_clinvar_classification
    ),
    has_clinvar = !is.na(clinvar_classification),
    clinvar_category = purrr::map_chr(
      clinvar_classification,
      collapse_clinvar_category
    ),
    variant_type_raw = purrr::map_chr(variants, extract_variant_type),
    variant_type = purrr::map_chr(
      variant_type_raw,
      standardize_variant_type
    ),
    genes = purrr::map_chr(variants, extract_genes),
    clinvar_phenotypes = purrr::map_chr(
      variants,
      extract_clinvar_phenotypes
    ),
    is_plp = has_clinvar &
      clinvar_category == "Likely pathogenic/Pathogenic"
  )

clinvar_variants <- annotated_variants |>
  dplyr::filter(has_clinvar)

clinvar_summary <- clinvar_variants |>
  dplyr::count(
    variant_origin,
    clinvar_category,
    name = "n_variants"
  ) |>
  dplyr::arrange(
    factor(variant_origin, c("known", "absent", "novel")),
    dplyr::desc(n_variants)
  )

plp_summary <- annotated_variants |>
  dplyr::filter(is_plp) |>
  dplyr::count(
    variant_origin,
    variant_type,
    variant_freq_class,
    name = "n_variants"
  ) |>
  dplyr::arrange(
    factor(variant_origin, c("known", "absent", "novel")),
    variant_type,
    variant_freq_class
  )

other_clinvar_terms <- clinvar_variants |>
  dplyr::filter(clinvar_category == "Other") |>
  dplyr::count(clinvar_classification, sort = TRUE, name = "n_variants")

# Create one gene row per unique P/LP variant-gene pair before counting. This
# prevents duplicate transcript annotations from inflating the word-cloud size.
plp_gene_counts <- annotated_variants |>
  dplyr::filter(is_plp, !is.na(genes), nzchar(genes)) |>
  dplyr::select(site_index, variant_origin, genes) |>
  tidyr::separate_rows(genes, sep = ";") |>
  dplyr::mutate(genes = trimws(genes)) |>
  dplyr::filter(nzchar(genes)) |>
  dplyr::distinct(site_index, variant_origin, genes) |>
  dplyr::count(variant_origin, genes, name = "n_variants", sort = TRUE)

# Apply a local, versionable Ensembl-to-HGNC map when supplied. Live biomaRt
# queries are intentionally excluded because they change across database builds.
if (file.exists(GENE_SYMBOL_MAP_TSV)) {
  gene_map <- readr::read_tsv(
    GENE_SYMBOL_MAP_TSV,
    show_col_types = FALSE,
    progress = FALSE
  )

  assert_columns(
    gene_map,
    c("ensembl_gene_id", "hgnc_symbol"),
    "Gene-symbol map"
  )

  plp_gene_counts <- plp_gene_counts |>
    dplyr::left_join(
      gene_map |>
        dplyr::distinct(ensembl_gene_id, .keep_all = TRUE),
      by = c("genes" = "ensembl_gene_id")
    ) |>
    dplyr::mutate(
      genes = dplyr::if_else(
        !is.na(hgnc_symbol) & nzchar(hgnc_symbol),
        hgnc_symbol,
        genes
      )
    ) |>
    dplyr::select(-hgnc_symbol) |>
    dplyr::group_by(variant_origin, genes) |>
    dplyr::summarise(n_variants = sum(n_variants), .groups = "drop")
}

unresolved_ensembl_genes <- plp_gene_counts |>
  dplyr::filter(stringr::str_starts(genes, "ENSG"))

plp_gene_counts_with_symbols <- plp_gene_counts |>
  dplyr::filter(!stringr::str_starts(genes, "ENSG"))

write_tsv(
  clinvar_summary,
  file.path(TABLES_DIR, "04_clinvar_category_counts.tsv")
)

write_tsv(
  plp_summary,
  file.path(TABLES_DIR, "04_plp_counts_by_origin_type_and_frequency.tsv")
)

write_tsv(
  plp_gene_counts,
  file.path(TABLES_DIR, "04_plp_gene_counts.tsv")
)

write_tsv(
  unresolved_ensembl_genes,
  file.path(TABLES_DIR, "04_plp_unresolved_ensembl_genes.tsv")
)

write_tsv(
  other_clinvar_terms,
  file.path(REPORTS_DIR, "04_unconsolidated_clinvar_terms.tsv")
)

write_workbook(
  list(
    ClinVar_summary = clinvar_summary,
    PLP_summary = plp_summary,
    PLP_gene_counts = plp_gene_counts,
    Known_PLP = drop_list_columns(
      annotated_variants |>
        dplyr::filter(variant_origin == "known", is_plp)
    ),
    Absent_PLP = drop_list_columns(
      annotated_variants |>
        dplyr::filter(variant_origin == "absent", is_plp)
    ),
    Novel_PLP = drop_list_columns(
      annotated_variants |>
        dplyr::filter(variant_origin == "novel", is_plp)
    )
  ),
  file.path(TABLES_DIR, "04_clinvar_analysis.xlsx")
)

# ClinVar assertion distribution ---------------------------------------------

clinvar_plot_data <- clinvar_summary |>
  dplyr::mutate(
    Origin = format_origin(variant_origin),
    clinvar_category = factor(
      clinvar_category,
      levels = names(CLINVAR_COLORS)
    )
  )

clinvar_plot <- ggplot2::ggplot(
  clinvar_plot_data,
  ggplot2::aes(
    x = stats::reorder(clinvar_category, n_variants),
    y = n_variants,
    fill = clinvar_category
  )
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(n_variants)),
    hjust = -0.15,
    size = 2.8
  ) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ Origin, scales = "free_y") +
  ggplot2::scale_fill_manual(values = CLINVAR_COLORS, drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    expand = ggplot2::expansion(mult = c(0, 0.14))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of variants",
    fill = NULL,
    title = "ClinVar germline assertions in aPRoVAR"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(color = "black"),
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

save_plot_formats(
  clinvar_plot,
  "ClinVar_assertions_by_variant_origin",
  width = 12,
  height = 7
)

# P/LP gene word clouds --------------------------------------------------------

phenotype_adjustment_file <- file.path(
  RAW_DATA_DIR,
  "genes_for_AF_phenotype_adjustment.txt"
)

assert_file_exists(
  phenotype_adjustment_file,
  "Phenotype-adjustment gene list"
)

# The file is expected to contain one gene symbol per line.
phenotype_adjustment_genes <- readr::read_lines(
  phenotype_adjustment_file,
  progress = FALSE
) |>
  trimws()

phenotype_adjustment_genes <- unique(
  toupper(
    phenotype_adjustment_genes[
      nzchar(phenotype_adjustment_genes) &
        !startsWith(phenotype_adjustment_genes, "#")
    ]
  )
)

# Apply the phenotype-adjustment exclusion only to the word-cloud data.
set.seed(RANDOM_SEED)

wordcloud_data <- plp_gene_counts_with_symbols |>
  dplyr::filter(
    !toupper(genes) %in% phenotype_adjustment_genes
  ) |>
  dplyr::mutate(
    Origin = format_origin(variant_origin),
    cor_categoria = dplyr::case_when(
      n_variants == 3L ~ "red",
      n_variants == 2L ~ "orange",
      n_variants == 1L ~ "darkgrey",
      TRUE ~ "black"
    )
  ) |>
  dplyr::group_by(Origin) |>
  dplyr::arrange(
    dplyr::desc(n_variants),
    genes,
    .by_group = TRUE
  ) |>
  dplyr::mutate(
    angle = 90 * sample(
      c(0, 1),
      dplyr::n(),
      replace = TRUE,
      prob = c(0.70, 0.30)
    )
  ) |>
  dplyr::ungroup()

writexl::write_xlsx(
  list(
    wordcloud_data = wordcloud_data,
    plp_gene_counts_all = plp_gene_counts,
    plp_gene_counts_with_hugo = plp_gene_counts_with_symbols
  ),
  file.path(TABLES_DIR, "06_wordcloud_data.xlsx")
)

if (
  nrow(wordcloud_data) > 0L &&
  requireNamespace("ggwordcloud", quietly = TRUE)
) {
  set.seed(RANDOM_SEED)
  
  wordcloud_data <- wordcloud_data |>
    dplyr::group_by(Origin) |>
    dplyr::mutate(
      angle = 90 * sample(
        c(0, 1),
        dplyr::n(),
        replace = TRUE,
        prob = c(0.75, 0.25)
      )
    ) |>
    dplyr::ungroup()
  
  wordcloud_colors <- c(
    "1 variant" = "darkgrey",
    "2 variants" = "orange",
    "3 variants" = "red",
    ">3 variants" = "black"
  )
  
  # Use common size limits so that word sizes remain comparable across plots.
  word_size_limits <- range(
    wordcloud_data$n_variants,
    na.rm = TRUE
  )
  
  make_gene_wordcloud <- function(data, origin_name) {
    plot_data <- data |>
      dplyr::filter(Origin == origin_name)
    
    if (nrow(plot_data) == 0L) {
      return(NULL)
    }
    
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        label = genes,
        size = n_variants,
        color = count_category,
        angle = angle
      )
    ) +
      ggwordcloud::geom_text_wordcloud_area(
        grid_size = 2,
        eccentricity = 1,
        rm_outside = FALSE
      ) +
      ggplot2::scale_size_area(
        max_size = 20,
        limits = word_size_limits
      ) +
      ggplot2::scale_color_manual(
        values = wordcloud_colors,
        breaks = names(wordcloud_colors),
        drop = FALSE
      ) +
      ggplot2::guides(
        size = "none"
      ) +
      ggplot2::labs(
        color = "P/LP variants per gene",
        title = paste0(
          origin_name,
          " genes harboring P/LP variants"
        )
      ) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      )
  }
  
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
          "ClinVar_PLP_gene_wordcloud_",
          tolower(origin_name)
        ),
        width = 6,
        height = 6
      )
    }
  }
} else {
  warning(
    "Gene word clouds were not generated because ggwordcloud is unavailable, ",
    "no eligible HGNC-symbol rows remained, or every gene was excluded by ",
    "the phenotype-adjustment list."
  )
}


clinvar_results <- list(
  annotated_variants = annotated_variants,
  clinvar_summary = clinvar_summary,
  plp_summary = plp_summary,
  plp_gene_counts = plp_gene_counts,
  other_clinvar_terms = other_clinvar_terms
)

write_rds_parallel(clinvar_results, CLINVAR_RESULTS_RDS)
message("ClinVar analysis complete.")
