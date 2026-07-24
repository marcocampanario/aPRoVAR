##### 01_prepare_variant_data.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Read the parsed aPRoVAR multisample variant object, extract VCF-level metrics,
# split multiallelic sites, validate and align ABRaOM/ANNOVAR annotations, and
# construct the sample-by-site presence matrix used downstream.

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

require_packages(c("data.table", "dplyr", "purrr", "tibble"))
create_output_directories()

assert_file_exists(HEADER_RDS, "Header RDS file")
assert_file_exists(VARIANT_RDS_GZ, "Parsed variant RDS file")
assert_file_exists(ABRAOM_ANNOVAR_TSV, "ABRaOM ANNOVAR table")

message("Reading the aPRoVAR header and parsed variant data...")
header <- readRDS(HEADER_RDS)
site_table <- read_rds_parallel(VARIANT_RDS_GZ)
site_table <- tibble::as_tibble(site_table)

assert_columns(
  site_table,
  c(
    "chromosome",
    "position",
    "refAllele",
    "altAlleles",
    "samples",
    "variants",
    "vcfInfo"
  ),
  "Parsed variant table"
)

sample_names <- extract_sample_names(header)

if (length(sample_names) != EXPECTED_SAMPLE_COUNT) {
  warning(
    "The header contains ",
    length(sample_names),
    " samples; the configured expectation is ",
    EXPECTED_SAMPLE_COUNT,
    "."
  )
}

message("Extracting F_MISSING, AF, and AF_EXCL from vcfInfo...")
site_table <- site_table |>
  dplyr::mutate(
    site_index = dplyr::row_number(),
    F_MISSING = purrr::map_dbl(
      vcfInfo,
      ~ suppressWarnings(
        as.numeric(extract_vcf_info(.x, "F_MISSING"))
      )
    ),
    AF = purrr::map(
      vcfInfo,
      ~ parse_af_vector(extract_vcf_info(.x, "AF"))
    ),
    AF_EXCL = purrr::map(
      vcfInfo,
      ~ parse_af_vector(extract_vcf_info(.x, "AF_EXCL"))
    ),
    site_id = make_variant_id(
      chromosome,
      position,
      refAllele,
      altAlleles
    ),
    is_biallelic = lengths(altAlleles) == 1L
  ) |>
  dplyr::select(-vcfInfo)

if (anyDuplicated(site_table$site_id)) {
  warning(
    "Duplicate site IDs were detected. site_index will remain the primary key."
  )
}

message("Splitting multiallelic sites into one row per ALT allele...")
allele_table <- purrr::map_dfr(seq_len(nrow(site_table)), function(row_index) {
  alt_alleles <- as.character(unlist(site_table$altAlleles[[row_index]]))

  if (length(alt_alleles) == 0L) {
    return(tibble::tibble())
  }

  matched_annotations <- match_variant_annotations_to_alts(
    site_table$variants[[row_index]],
    length(alt_alleles)
  )

  tibble::tibble(
    site_index = site_table$site_index[[row_index]],
    allele_index = seq_along(alt_alleles),
    chromosome = as.character(site_table$chromosome[[row_index]]),
    position = as.integer(site_table$position[[row_index]]),
    refAllele = as.character(site_table$refAllele[[row_index]]),
    altAllele = alt_alleles,
    variant = matched_annotations
  )
}) |>
  dplyr::mutate(
    annovar_row_index = dplyr::row_number(),
    variant_id = paste(
      chromosome,
      position,
      refAllele,
      altAllele,
      sep = "-"
    ),
    af_gnomad_genome = purrr::map_dbl(
      variant,
      extract_population_all_af,
      annotation_name = "gnomad"
    ),
    af_gnomad_exome = purrr::map_dbl(
      variant,
      extract_population_all_af,
      annotation_name = "gnomad-exome"
    ),
    af_gnomad = purrr::map2_dbl(
      af_gnomad_genome,
      af_gnomad_exome,
      select_preferred_gnomad_af
    ),
    af_1000g = purrr::map_dbl(
      variant,
      extract_population_all_af,
      annotation_name = "oneKg"
    ),
    dbsnp_id = purrr::map_chr(variant, extract_dbsnp_id)
  ) |>
  dplyr::select(-variant)

message("Reading ABRaOM and local ANNOVAR frequency annotations...")
annovar_table <- data.table::fread(
  ABRAOM_ANNOVAR_TSV,
  sep = "\t",
  quote = "",
  data.table = FALSE,
  showProgress = FALSE
)

assert_columns(
  annovar_table,
  c("Chr", "Start", "Ref", "Alt", "abraom_freq", "Otherinfo1"),
  "ABRaOM ANNOVAR table"
)

message("Validating same-row correspondence with the ANNOVAR table...")
annovar_alignment <- validate_annovar_row_alignment(
  allele_table = allele_table,
  annovar_table = annovar_table,
  minimum_exact_match = ANNOVAR_MIN_EXACT_ROW_MATCH,
  minimum_standard_match = ANNOVAR_MIN_STANDARD_ROW_MATCH
)

message(
  sprintf(
    paste0(
      "ANNOVAR alignment validated: %.3f%% exact and %.3f%% standard-",
      "normalized same-row keys."
    ),
    100 * annovar_alignment$exact_key_match_rate,
    100 * annovar_alignment$standard_key_match_rate
  )
)

annovar_frequency <- annovar_table |>
  dplyr::transmute(
    annovar_row_index = dplyr::row_number(),
    af_abraom = parse_annovar_numeric(abraom_freq),
    af_local_annovar = parse_annovar_numeric(Otherinfo1)
  )

# ANNOVAR rewrites many INDEL and complex-allele representations, so direct
# CHR-POS-REF-ALT matching is incomplete. The row index is used only after the
# strict correspondence checks above establish that both tables derive from
# the same split VCF in the same order.
allele_table <- allele_table |>
  dplyr::left_join(annovar_frequency, by = "annovar_row_index")

annovar_alignment_metrics <- data.frame(
  metric = c(
    "ALT-allele records",
    "ANNOVAR records",
    "Same-row chromosome matches",
    "Same-row chromosome match percent",
    "Same-row exact variant-key matches",
    "Same-row exact variant-key match percent",
    "Same-row standard-normalized key matches",
    "Same-row standard-normalized key match percent"
  ),
  value = c(
    annovar_alignment$allele_records,
    annovar_alignment$annovar_records,
    annovar_alignment$same_chromosome_count,
    100 * annovar_alignment$same_chromosome_rate,
    annovar_alignment$exact_key_match_count,
    100 * annovar_alignment$exact_key_match_rate,
    annovar_alignment$standard_key_match_count,
    100 * annovar_alignment$standard_key_match_rate
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  annovar_alignment_metrics,
  file.path(REPORTS_DIR, "01_annovar_alignment_metrics.tsv")
)

message("Constructing the variant-presence matrix...")
presence_matrix <- make_presence_matrix(
  sample_list = site_table$samples,
  sample_names = sample_names,
  row_ids = as.character(site_table$site_index)
)

preparation_metrics <- data.frame(
  metric = c(
    "Samples in header",
    "Total unsplit variant sites",
    "Split ALT-allele records",
    "Sites with F_MISSING equal to zero",
    "Sites with F_MISSING greater than 1%",
    "Sites with F_MISSING greater than 10%",
    "Sites with phenotype-adjusted AF_EXCL",
    "ALT alleles aligned to ANNOVAR",
    "ALT alleles with non-missing ABRaOM AF"
  ),
  value = c(
    length(sample_names),
    nrow(site_table),
    nrow(allele_table),
    sum(site_table$F_MISSING == 0, na.rm = TRUE),
    sum(site_table$F_MISSING > 0.01, na.rm = TRUE),
    sum(site_table$F_MISSING > 0.10, na.rm = TRUE),
    sum(!purrr::map_lgl(site_table$AF_EXCL, is_empty_value)),
    annovar_alignment$annovar_records,
    sum(!is.na(allele_table$af_abraom))
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  preparation_metrics,
  file.path(REPORTS_DIR, "01_preparation_metrics.tsv")
)

prepared_data <- list(
  header = header,
  sample_names = sample_names,
  site_table = site_table,
  allele_table = allele_table,
  presence_matrix = presence_matrix,
  preparation_metrics = preparation_metrics,
  annovar_alignment = annovar_alignment,
  annovar_alignment_metrics = annovar_alignment_metrics
)

message("Saving the prepared-data checkpoint...")
write_rds_parallel(prepared_data, PREPARED_DATA_RDS)
message("Preparation complete: ", PREPARED_DATA_RDS)
