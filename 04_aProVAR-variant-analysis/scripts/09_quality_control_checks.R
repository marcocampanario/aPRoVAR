##### 09_quality_control_checks.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Run final structural and internal-consistency checks on the modular aPRoVAR
# outputs and record the R session information used for reproducibility.

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

require_packages(c("dplyr", "purrr"))
create_output_directories()

required_checkpoints <- c(
  PREPARED_DATA_RDS,
  ANALYSIS_READY_RDS,
  DISCOVERY_RESULTS_RDS,
  CLINVAR_RESULTS_RDS,
  PHARMACOGENETIC_RESULTS_RDS,
  PLOF_RESULTS_RDS
)

missing_checkpoints <- required_checkpoints[!file.exists(required_checkpoints)]

if (length(missing_checkpoints) > 0L) {
  stop(
    "The following checkpoint(s) are missing: ",
    paste(missing_checkpoints, collapse = ", "),
    call. = FALSE
  )
}

analysis_ready <- read_rds_parallel(ANALYSIS_READY_RDS)
variant_table <- analysis_ready$analysis_table
presence_matrix <- analysis_ready$presence_matrix
annovar_alignment <- analysis_ready$annovar_alignment

analysis_site_count <- nrow(variant_table)
analysis_sample_count <- ncol(presence_matrix)
all_fmissing_zero <- all(variant_table$F_MISSING == 0)
matrix_rows_match <- analysis_site_count == nrow(presence_matrix)
valid_origins <- all(
  variant_table$variant_origin %in% c("known", "absent", "novel")
)
valid_frequency_classes <- all(
  is.na(variant_table$variant_freq_class) |
    variant_table$variant_freq_class %in% c(
      "singleton",
      "doubleton",
      "doubleton_hom",
      "rare",
      "polymorphic"
    )
)
unique_site_index <- !anyDuplicated(variant_table$site_index)
valid_annovar_alignment <-
  !is.null(annovar_alignment) &&
  isTRUE(annovar_alignment$same_chromosome_rate == 1) &&
  isTRUE(
    annovar_alignment$exact_key_match_rate >=
      annovar_alignment$minimum_exact_match
  ) &&
  isTRUE(
    annovar_alignment$standard_key_match_rate >=
      annovar_alignment$minimum_standard_match
  )
known_count <- sum(variant_table$variant_origin == "known")
absent_count <- sum(variant_table$variant_origin == "absent")
novel_count <- sum(variant_table$variant_origin == "novel")

rm(analysis_ready, variant_table, presence_matrix, annovar_alignment)
invisible(gc())

discovery_results <- read_rds_parallel(DISCOVERY_RESULTS_RDS)
cumulative_results <- discovery_results$cumulative_results
final_cumulative_count <- tail(cumulative_results$total, 1L)
rm(discovery_results, cumulative_results)
invisible(gc())

clinvar_results <- read_rds_parallel(CLINVAR_RESULTS_RDS)
annotated_variants <- clinvar_results$annotated_variants
clinvar_row_count <- nrow(annotated_variants)
clinvar_count <- sum(annotated_variants$has_clinvar)
plp_count <- sum(annotated_variants$is_plp)
rm(clinvar_results, annotated_variants)
invisible(gc())

plof_results <- read_rds_parallel(PLOF_RESULTS_RDS)
plof_variants <- plof_results$annotated_variants
plof_row_count <- nrow(plof_variants)
plof_count <- sum(plof_variants$is_plof)
plof_plp_count <- sum(plof_variants$is_plof_plp)
rm(plof_results, plof_variants)
invisible(gc())

qc_checks <- data.frame(
  check = c(
    "All analysis sites have F_MISSING = 0",
    "Analysis table and presence matrix have equal row counts",
    "Presence matrix has the configured number of samples",
    "Variant-origin classes are limited to known/absent/novel",
    "Frequency classes use only the documented levels",
    "Final cumulative count equals the number of analysis sites",
    "ClinVar table preserves the number of analysis sites",
    "pLOF table preserves the number of analysis sites",
    "Site index is unique in the analysis table",
    "ANNOVAR row alignment passed every configured safeguard"
  ),
  passed = c(
    all_fmissing_zero,
    matrix_rows_match,
    analysis_sample_count == EXPECTED_SAMPLE_COUNT,
    valid_origins,
    valid_frequency_classes,
    final_cumulative_count == analysis_site_count,
    clinvar_row_count == analysis_site_count,
    plof_row_count == analysis_site_count,
    unique_site_index,
    valid_annovar_alignment
  ),
  stringsAsFactors = FALSE
)

qc_metrics <- data.frame(
  metric = c(
    "Analysis-ready sites",
    "Samples",
    "Known sites",
    "Absent sites",
    "Novel sites",
    "Sites with any ClinVar assertion",
    "Likely pathogenic/pathogenic sites",
    "Strict pLOF sites",
    "Strict pLOF and P/LP sites"
  ),
  value = c(
    analysis_site_count,
    analysis_sample_count,
    known_count,
    absent_count,
    novel_count,
    clinvar_count,
    plp_count,
    plof_count,
    plof_plp_count
  ),
  stringsAsFactors = FALSE
)

write_tsv(qc_checks, file.path(REPORTS_DIR, "09_quality_control_checks.tsv"))
write_tsv(qc_metrics, file.path(REPORTS_DIR, "09_quality_control_metrics.tsv"))
writeLines(capture.output(sessionInfo()), file.path(REPORTS_DIR, "sessionInfo.txt"))

if (!all(qc_checks$passed)) {
  failed_checks <- qc_checks$check[!qc_checks$passed]
  stop(
    "Quality control failed: ",
    paste(failed_checks, collapse = "; "),
    call. = FALSE
  )
}

message("All final quality-control checks passed.")
