##### 02_classify_variants.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Classify aPRoVAR sites by population-database status and allele-frequency
# class, retain complete-genotype sites, and create the analysis-ready checkpoint.

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

require_packages(c("dplyr", "purrr", "tibble"))
create_output_directories()
assert_file_exists(PREPARED_DATA_RDS, "Prepared-data checkpoint")

message("Reading the prepared-data checkpoint...")
prepared_data <- read_rds_parallel(PREPARED_DATA_RDS)
site_table <- prepared_data$site_table
allele_table <- prepared_data$allele_table

assert_columns(
  allele_table,
  c(
    "af_gnomad_genome",
    "af_gnomad_exome",
    "af_gnomad",
    "af_1000g",
    "af_abraom",
    "dbsnp_id"
  ),
  "Prepared ALT-allele table"
)

message("Classifying ALT alleles by population-database status...")
allele_table <- allele_table |>
  dplyr::mutate(
    is_known = is_known_from_population_af(
      gnomad_genome_af = af_gnomad_genome,
      gnomad_exome_af = af_gnomad_exome,
      one_kg_af = af_1000g,
      abraom_af = af_abraom
    ),
    variant_origin = purrr::map2_chr(
      is_known,
      dbsnp_id,
      classify_origin
    )
  )

message("Aggregating ALT-allele annotations to the original site level...")
site_origin <- allele_table |>
  dplyr::group_by(site_index) |>
  dplyr::summarise(
    is_known = any(is_known, na.rm = TRUE),
    dbsnp_id = collapse_unique(dbsnp_id, separator = ","),
    af_gnomad_genome = list(af_gnomad_genome),
    af_gnomad_exome = list(af_gnomad_exome),
    af_gnomad = list(af_gnomad),
    af_1000g = list(af_1000g),
    af_abraom = list(af_abraom),
    af_local_annovar = list(af_local_annovar),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    variant_origin = purrr::map2_chr(
      is_known,
      dbsnp_id,
      classify_origin
    ),
    is_novel = variant_origin == "novel"
  )

message("Counting called alleles and assigning frequency classes...")
site_table <- site_table |>
  dplyr::left_join(site_origin, by = "site_index") |>
  dplyr::mutate(
    AF_used = purrr::map2(AF, AF_EXCL, choose_analysis_af),
    uses_adjusted_af = !purrr::map_lgl(AF_EXCL, is_empty_value),
    af_local = purrr::map_dbl(AF_used, function(value) {
      value <- suppressWarnings(as.numeric(value))
      value <- value[is.finite(value)]
      if (length(value) == 0L) NA_real_ else max(value)
    }),
    allele_summary = purrr::map(samples, count_site_alleles),
    allele_number = purrr::map_int(
      allele_summary,
      ~ if (is.null(.x$allele_number)) NA_integer_ else .x$allele_number
    ),
    alt_counts = purrr::map(allele_summary, "alt_counts"),
    total_alt_count = purrr::map_int(
      allele_summary,
      ~ if (is.null(.x$total_alt_count)) NA_integer_ else .x$total_alt_count
    ),
    is_doubleton_hom = purrr::map_lgl(
      allele_summary,
      ~ isTRUE(.x$is_doubleton_hom)
    ),
    variant_freq_class = purrr::pmap_chr(
      list(alt_counts, allele_number, is_doubleton_hom),
      classify_frequency_from_counts
    )
  ) |>
  dplyr::select(-allele_summary)

# The manuscript's descriptive analyses use only sites with complete genotype
# calls across all 1,010 individuals.
analysis_table <- site_table |>
  dplyr::filter(F_MISSING == 0)

origin_summary <- analysis_table |>
  dplyr::count(variant_origin, name = "n_variants") |>
  dplyr::arrange(factor(variant_origin, c("known", "absent", "novel")))

origin_frequency_summary <- analysis_table |>
  dplyr::count(
    variant_origin,
    variant_freq_class,
    name = "n_variants"
  ) |>
  dplyr::arrange(
    factor(variant_origin, c("known", "absent", "novel")),
    factor(
      variant_freq_class,
      c("singleton", "doubleton", "doubleton_hom", "rare", "polymorphic")
    )
  )

write_tsv(
  origin_summary,
  file.path(TABLES_DIR, "02_variant_origin_counts.tsv")
)

write_tsv(
  origin_frequency_summary,
  file.path(TABLES_DIR, "02_variant_origin_by_frequency_class.tsv")
)

classification_report <- c(
  "aPRoVAR VARIANT CLASSIFICATION REPORT",
  "=====================================",
  "",
  paste("Total sites before the complete-genotype filter:", nrow(site_table)),
  paste("Sites with F_MISSING = 0:", nrow(analysis_table)),
  "",
  "Population-database status among F_MISSING = 0 sites:",
  paste0(
    "- ",
    origin_summary$variant_origin,
    ": ",
    origin_summary$n_variants
  )
)

write_report(
  classification_report,
  file.path(REPORTS_DIR, "02_variant_classification_report.txt")
)

analysis_ready <- list(
  header = prepared_data$header,
  sample_names = prepared_data$sample_names,
  analysis_table = analysis_table,
  allele_table = allele_table,
  presence_matrix = prepared_data$presence_matrix[
    as.character(analysis_table$site_index),
    ,
    drop = FALSE
  ],
  origin_summary = origin_summary,
  origin_frequency_summary = origin_frequency_summary,
  annovar_alignment = prepared_data$annovar_alignment,
  annovar_alignment_metrics = prepared_data$annovar_alignment_metrics
)

message("Saving the analysis-ready checkpoint...")
write_rds_parallel(analysis_ready, ANALYSIS_READY_RDS)
message("Variant classification complete: ", ANALYSIS_READY_RDS)
