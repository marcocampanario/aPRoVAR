##### 06_plof_analysis.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Identify predicted loss-of-function variants from high-impact MANE Select
# RefSeq consequences, retain REVEL as a separate missense-prioritization field,
# and build the ClinVar/pLOF summary tables reported for aPRoVAR.

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

require_packages(c("dplyr", "purrr", "tidyr"))
create_output_directories()
assert_file_exists(CLINVAR_RESULTS_RDS, "ClinVar-results checkpoint")

clinvar_results <- read_rds_parallel(CLINVAR_RESULTS_RDS)
variant_table <- clinvar_results$annotated_variants

message("Extracting selected transcript impacts and REVEL scores...")
plof_annotated_variants <- variant_table |>
  dplyr::mutate(
    canonical_ensembl_impact = purrr::map_chr(
      variants,
      extract_selected_transcript_impact,
      source = "Ensembl",
      require_canonical = TRUE,
      require_mane = FALSE
    ),
    mane_ensembl_impact = purrr::map_chr(
      variants,
      extract_selected_transcript_impact,
      source = "Ensembl",
      require_canonical = TRUE,
      require_mane = TRUE
    ),
    mane_refseq_impact = purrr::map_chr(
      variants,
      extract_selected_transcript_impact,
      source = "RefSeq",
      require_canonical = TRUE,
      require_mane = TRUE
    ),
    revel_max = purrr::map_dbl(variants, extract_revel_max),
    is_high_revel = !is.na(revel_max) & revel_max >= REVEL_HIGH_THRESHOLD,
    is_plof = !is.na(mane_refseq_impact) &
      grepl("(^|;)high(;|$)", mane_refseq_impact),
    legacy_high_impact_or_revel = is_plof | is_high_revel,
    is_plof_plp = is_plof & is_plp
  )

plof_annotated_variants <- plof_annotated_variants |>
  dplyr::mutate(
    is_plof_plp_pd = legacy_high_impact_or_revel & is_plp
  )

# The original exploratory script used high MANE impact OR high REVEL under the
# pLOF label. REVEL predicts deleterious missense effects and is not itself a
# loss-of-function definition; both flags are retained here, but pLOF is based
# only on high-impact MANE Select RefSeq consequences.
# For aPRoVAR 1.0, we used the tag pLoF/pD for legacy_high_impact_or_revel AND is_plp
definition_audit <- data.frame(
  definition = c(
    "Strict pLOF: high-impact MANE Select RefSeq consequence",
    paste0("High REVEL score (≥ ", REVEL_HIGH_THRESHOLD, ")"),
    "Legacy composite: strict pLOF OR high REVEL",
    "GVs strictly pLOF or with REVEL ≥ 0.773, and P/LP"
  ),
  n_variants = c(
    sum(plof_annotated_variants$is_plof, na.rm = TRUE),
    sum(plof_annotated_variants$is_high_revel, na.rm = TRUE),
    sum(
      plof_annotated_variants$legacy_high_impact_or_revel,
      na.rm = TRUE
    ),
    sum(
      plof_annotated_variants$is_plof_plp_pd,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

count_variant_set <- function(data, set_name) {
  data |>
    dplyr::count(
      variant_origin,
      variant_type,
      variant_freq_class,
      name = "n_variants"
    ) |>
    dplyr::mutate(variant_set = set_name, .before = 1L)
}

variant_sets <- list(
  `Any ClinVar assertion` = plof_annotated_variants |>
    dplyr::filter(has_clinvar),
  `Likely pathogenic/Pathogenic` = plof_annotated_variants |>
    dplyr::filter(is_plp),
  `pLOF/pD` = plof_annotated_variants |>
    dplyr::filter(legacy_high_impact_or_revel),
  `pLOF/pD and Likely pathogenic/Pathogenic` = plof_annotated_variants |>
    dplyr::filter(is_plof_plp_pd)
)

summary_counts <- dplyr::bind_rows(
  purrr::imap(variant_sets, count_variant_set)
) |>
  dplyr::arrange(
    factor(
      variant_origin,
      c("known", "absent", "novel")
    ),
    variant_set,
    variant_type,
    variant_freq_class
  )

summary_totals <- summary_counts |>
  dplyr::group_by(variant_set, variant_origin, variant_type) |>
  dplyr::summarise(
    n_variants = sum(n_variants),
    .groups = "drop"
  )

write_tsv(
  definition_audit,
  file.path(REPORTS_DIR, "06_plof_definition_audit.tsv")
)

write_tsv(
  summary_counts,
  file.path(TABLES_DIR, "06_clinvar_plof_counts_by_frequency.tsv")
)

write_tsv(
  summary_totals,
  file.path(TABLES_DIR, "06_clinvar_plof_totals.tsv")
)

variant_export_sheets <- c(
  lapply(variant_sets, drop_list_columns),
  list(
    Summary_counts = summary_counts,
    Summary_totals = summary_totals,
    Definition_audit = definition_audit
  )
)

names(variant_export_sheets)[seq_along(variant_sets)] <- c(
  "ClinVar",
  "PLP",
  "pLOF_pD",
  "pLOF_pD_PLP"
)

write_workbook(
  variant_export_sheets,
  file.path(TABLES_DIR, "06_clinvar_plof_variant_sets.xlsx")
)

# Detailed annotation exports for potentially novel pLOF/pD P/LP variants --------

novel_plof_plp <- plof_annotated_variants |>
  dplyr::filter(variant_origin == "novel", is_plof_plp_pd)

novel_plof_plp_transcripts <- purrr::map2_dfr(
  novel_plof_plp$site_id,
  novel_plof_plp$variants,
  function(site_id, variant) {
    transcripts <- extract_transcript_table(variant)

    if (is.null(transcripts) || nrow(transcripts) == 0L) {
      return(data.frame())
    }

    dplyr::mutate(transcripts, site_id = site_id, .before = 1L)
  }
)

novel_plof_plp_consequences <- if (
  nrow(novel_plof_plp_transcripts) > 0L &&
    "consequence" %in% names(novel_plof_plp_transcripts)
) {
  novel_plof_plp_transcripts |>
    dplyr::select(site_id, consequence) |>
    tidyr::unnest_longer(consequence) |>
    dplyr::distinct()
} else {
  data.frame(site_id = character(), consequence = character())
}

write_workbook(
  list(
    Novel_pLOF_PLP_sites = drop_list_columns(novel_plof_plp),
    Transcripts = drop_list_columns(novel_plof_plp_transcripts),
    Consequences = novel_plof_plp_consequences
  ),
  file.path(TABLES_DIR, "06_novel_plof_plp_annotations.xlsx")
)

plof_results <- list(
  annotated_variants = plof_annotated_variants,
  summary_counts = summary_counts,
  summary_totals = summary_totals,
  definition_audit = definition_audit,
  novel_plof_pd_plp = novel_plof_plp,
  novel_plof_pd_plp_transcripts = novel_plof_plp_transcripts,
  novel_plof_pd_plp_consequences = novel_plof_plp_consequences
)

write_rds_parallel(plof_results, PLOF_RESULTS_RDS)
message("pLOF analysis complete.")
