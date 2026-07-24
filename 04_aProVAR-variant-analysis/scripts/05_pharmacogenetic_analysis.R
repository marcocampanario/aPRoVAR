##### 05_pharmacogenetic_analysis.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Match aPRoVAR variants to ClinPGx variant-gene relationships, summarize the
# nine manuscript-prioritized pharmacogenes, and export clinically relevant
# ClinVar assertions for pharmacogenetic interpretation.

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
  "readr",
  "scales",
  "tidyr"
))
create_output_directories()
assert_file_exists(CLINVAR_RESULTS_RDS, "ClinVar-results checkpoint")
assert_file_exists(CLINPGX_RELATIONSHIPS_TSV, "ClinPGx relationships table")

clinvar_results <- read_rds_parallel(CLINVAR_RESULTS_RDS)
annotated_variants <- clinvar_results$annotated_variants

relationships <- readr::read_tsv(
  CLINPGX_RELATIONSHIPS_TSV,
  show_col_types = FALSE,
  progress = FALSE
)

assert_columns(
  relationships,
  c("Entity1_type", "Entity1_name", "Entity2_type", "Entity2_name"),
  "ClinPGx relationships table"
)

# Retain only true ClinPGx associations and extract every dbSNP rsID reported
# in Entity1_name. The relationship type and Entity2 are intentionally not
# restricted because relevant rsID associations are not limited to
# Variant-Gene records.
associated_rsids <- relationships |>
  dplyr::filter(Association == "associated") |>
  dplyr::transmute(
    dbsnp_id = stringr::str_extract_all(
      as.character(Entity1_name),
      stringr::regex("\\brs[0-9]+\\b", ignore_case = TRUE)
    )
  ) |>
  tidyr::unnest_longer(dbsnp_id) |>
  dplyr::mutate(dbsnp_id = tolower(trimws(dbsnp_id))) |>
  dplyr::filter(!is.na(dbsnp_id), nzchar(dbsnp_id)) |>
  dplyr::distinct()

# Extract and normalize the rsIDs observed in the aPRoVAR dataset.
variant_rsids <- annotated_variants |>
  dplyr::select(site_index, dbsnp_id) |>
  dplyr::transmute(
    site_index,
    dbsnp_id = stringr::str_extract_all(
      as.character(dbsnp_id),
      stringr::regex("\\brs[0-9]+\\b", ignore_case = TRUE)
    )
  ) |>
  tidyr::unnest_longer(dbsnp_id) |>
  dplyr::mutate(dbsnp_id = tolower(trimws(dbsnp_id))) |>
  dplyr::filter(!is.na(dbsnp_id), nzchar(dbsnp_id)) |>
  dplyr::distinct()

matched_associated_rsids <- variant_rsids |>
  dplyr::inner_join(associated_rsids, by = "dbsnp_id") |>
  dplyr::distinct(site_index, dbsnp_id)

all_pharmacogenetic_variants <- annotated_variants |>
  dplyr::semi_join(matched_associated_rsids, by = "site_index")

all_pharmacogenetic_variants_clinvar_variants <- all_pharmacogenetic_variants |>
  dplyr::filter(has_clinvar)

priority_matches <- matched_relationships |>
  dplyr::filter(pharmacogene %in% PRIORITY_PHARMACOGENES)

priority_pharmacogenetic_variants <- annotated_variants |>
  dplyr::semi_join(priority_matches, by = "site_index") |>
  dplyr::left_join(
    priority_matches |>
      dplyr::group_by(site_index) |>
      dplyr::summarise(
        priority_pharmacogene = paste(
          sort(unique(pharmacogene)),
          collapse = ";"
        ),
        .groups = "drop"
      ),
    by = "site_index"
  )

priority_clinvar_variants <- priority_pharmacogenetic_variants |>
  dplyr::filter(has_clinvar)

clinically_relevant_priority_variants <- priority_clinvar_variants |>
  dplyr::filter(
    clinvar_category %in% c(
      "Drug response",
      "Likely pathogenic/Pathogenic"
    )
  )

pharmacogene_counts <- priority_matches |>
  dplyr::distinct(site_index, pharmacogene) |>
  dplyr::count(pharmacogene, name = "n_variants", sort = TRUE) |>
  tidyr::complete(
    pharmacogene = PRIORITY_PHARMACOGENES,
    fill = list(n_variants = 0L)
  ) |>
  dplyr::arrange(factor(pharmacogene, PRIORITY_PHARMACOGENES))

all_pharmacogenetic_variants_clinvar_summary <- all_pharmacogenetic_variants_clinvar_variants |>
  dplyr::count(clinvar_category, name = "n_variants", sort = TRUE)

priority_clinvar_summary <- priority_clinvar_variants |>
  dplyr::count(clinvar_category, name = "n_variants", sort = TRUE)

overall_summary <- data.frame(
  metric = c(
    "Unique aPRoVAR variants matched to any ClinPGx pharmacogene",
    "Unique aPRoVAR variants in the nine priority pharmacogenes",
    "Priority variants with a ClinVar assertion",
    "Priority variants with a drug-response or P/LP assertion"
  ),
  value = c(
    dplyr::n_distinct(all_pharmacogenetic_variants$site_index),
    dplyr::n_distinct(priority_pharmacogenetic_variants$site_index),
    dplyr::n_distinct(priority_clinvar_variants$site_index),
    dplyr::n_distinct(clinically_relevant_priority_variants$site_index)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  overall_summary,
  file.path(TABLES_DIR, "05_pharmacogenetic_summary.tsv")
)

write_workbook(
  list(
    Summary = overall_summary,
    Priority_gene_counts = pharmacogene_counts,
    Priority_ClinVar_categories = priority_clinvar_summary,
    All_ClinVar_categories = all_pharmacogenetic_variants_clinvar_summary,
    Priority_variants = drop_list_columns(priority_pharmacogenetic_variants),
    Drug_response_or_PLP = drop_list_columns(
      clinically_relevant_priority_variants
    )
  ),
  file.path(TABLES_DIR, "05_pharmacogenetic_variants.xlsx")
)

# All pharmacogenes
if (nrow(all_pharmacogenetic_variants_clinvar_summary) > 0L) {
  pharmacogenetic_plot <- ggplot2::ggplot(
    all_pharmacogenetic_variants_clinvar_summary,
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
      size = 3.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = CLINVAR_COLORS) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of variants",
      fill = NULL,
      title = "ClinVar assertions in all pharmacogenes"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.y = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  save_plot_formats(
    pharmacogenetic_plot,
    "All_Pharmacogenetic_ClinVar_assertions",
    width = 9,
    height = 6
  )
}

# Prioritized pharmacogenes
if (nrow(priority_clinvar_summary) > 0L) {
  pharmacogenetic_plot <- ggplot2::ggplot(
    priority_clinvar_summary,
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
      size = 3.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = CLINVAR_COLORS) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of variants",
      fill = NULL,
      title = "ClinVar assertions in the nine prioritized pharmacogenes"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.y = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  save_plot_formats(
    pharmacogenetic_plot,
    "Prioritized_Pharmacogenetic_ClinVar_assertions",
    width = 9,
    height = 6
  )
}


pharmacogenetic_results <- list(
  associated_rsids = associated_rsids,
  all_pharmacogenetic_variants = all_pharmacogenetic_variants,
  priority_pharmacogenetic_variants = priority_pharmacogenetic_variants,
  clinically_relevant_priority_variants = clinically_relevant_priority_variants,
  matched_associated_rsids = matched_associated_rsids,
  priority_matches = priority_matches,
  pharmacogene_counts = pharmacogene_counts,
  pharmacogenes_clinvar_summary = all_pharmacogenetic_variants_clinvar_summary,
  priority_clinvar_summary = priority_clinvar_summary,
  overall_summary = overall_summary
)

write_rds_parallel(
  pharmacogenetic_results,
  PHARMACOGENETIC_RESULTS_RDS
)
message("Pharmacogenetic analysis complete.")

