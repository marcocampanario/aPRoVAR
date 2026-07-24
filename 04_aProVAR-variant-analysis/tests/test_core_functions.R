##### test_core_functions.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Run lightweight unit checks for core genotype, frequency, ClinVar, and
# statistical helper functions without requiring the full aPRoVAR dataset.

project_dir <- Sys.getenv("APROVAR_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) != 1L) {
    stop(
      "Set APROVAR_PROJECT_DIR or execute this test file with Rscript.",
      call. = FALSE
    )
  }

  project_dir <- dirname(dirname(
    normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  ))
}

Sys.setenv(APROVAR_PROJECT_DIR = project_dir)
source(file.path(project_dir, "R", "00_config.R"))
source(file.path(project_dir, "R", "01_functions.R"))

# Genotype parsing and allele counting ----------------------------------------

sample_data <- data.frame(
  genotype = c("0/0", "0/1", "1/1", "./."),
  stringsAsFactors = FALSE
)

allele_summary <- count_site_alleles(sample_data)

stopifnot(
  identical(allele_summary$allele_number, 6L),
  identical(allele_summary$total_alt_count, 3L),
  identical(allele_summary$alt_counts, 3L),
  !allele_summary$is_doubleton_hom,
  !genotype_has_alt("0/0"),
  genotype_has_alt("0/1"),
  genotype_has_alt("2|0"),
  !genotype_has_alt("./.")
)

# Allele-frequency classes ----------------------------------------------------

stopifnot(
  classify_frequency_from_counts(1L, 2020L) == "singleton",
  classify_frequency_from_counts(2L, 2020L) == "doubleton",
  classify_frequency_from_counts(2L, 2020L, TRUE) == "doubleton_hom",
  classify_frequency_from_counts(3L, 2020L) == "rare",
  classify_frequency_from_counts(21L, 2020L) == "polymorphic"
)

# Population-database origin --------------------------------------------------

stopifnot(
  classify_origin(TRUE, NA_character_) == "known",
  classify_origin(FALSE, "rs123") == "absent",
  classify_origin(FALSE, NA_character_) == "novel"
)

# ANNOVAR row-alignment keys ---------------------------------------------------

stopifnot(
  make_standard_annovar_key_from_vcf("chr1", 358078L, "G", "GACGC") ==
    "1|358078|-|ACGC",
  make_standard_annovar_key_from_vcf("chr1", 451369L, "TGAA", "T") ==
    "1|451370|GAA|-",
  make_standard_annovar_key_from_vcf("chr2", 100L, "A", "G") ==
    "2|100|A|G"
)

synthetic_alleles <- data.frame(
  chromosome = c("chr1", "chr1", "chr2"),
  position = c(100L, 200L, 300L),
  refAllele = c("A", "TGAA", "C"),
  altAllele = c("G", "T", "CT"),
  stringsAsFactors = FALSE
)

synthetic_annovar <- data.frame(
  Chr = c("chr1", "chr1", "chr2"),
  Start = c(100L, 201L, 300L),
  Ref = c("A", "GAA", "-"),
  Alt = c("G", "-", "T"),
  stringsAsFactors = FALSE
)

synthetic_alignment <- validate_annovar_row_alignment(
  synthetic_alleles,
  synthetic_annovar,
  minimum_exact_match = 1 / 3,
  minimum_standard_match = 1
)

stopifnot(
  synthetic_alignment$same_chromosome_rate == 1,
  synthetic_alignment$exact_key_match_rate == 1 / 3,
  synthetic_alignment$standard_key_match_rate == 1,
  inherits(
    try(
      validate_annovar_row_alignment(
        synthetic_alleles[-1L, ],
        synthetic_annovar
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      validate_annovar_row_alignment(
        synthetic_alleles,
        transform(synthetic_annovar, Chr = c("chr9", "chr1", "chr2")),
        minimum_exact_match = 0,
        minimum_standard_match = 0
      ),
      silent = TRUE
    ),
    "try-error"
  )
)

# Nested annotation extraction ------------------------------------------------

synthetic_variant <- list(
  gnomad = data.frame(allAf = c(0.01, 0.03)),
  `gnomad-exome` = data.frame(allAf = 0.05),
  dbsnp = c("rs123", "rs123"),
  revel = list(score = c(0.4, 0.9)),
  transcripts = data.frame(
    source = c("RefSeq", "Ensembl"),
    isCanonical = c(TRUE, TRUE),
    isManeSelect = c(TRUE, TRUE),
    impact = c("high", "moderate"),
    hgnc = c("GENE1", "GENE1"),
    stringsAsFactors = FALSE
  )
)

stopifnot(
  extract_population_all_af(synthetic_variant, "gnomad") == 0.02,
  extract_population_all_af(synthetic_variant, "gnomad-exome") == 0.05,
  select_preferred_gnomad_af(0.02, 0.05) == 0.02,
  select_preferred_gnomad_af(NA_real_, 0.05) == 0.05,
  select_preferred_gnomad_af(0, 0.05) == 0,
  is_known_from_population_af(0, 0.05, NA_real_, NA_real_),
  !is_known_from_population_af(0, 0, NA_real_, NA_real_),
  extract_dbsnp_id(synthetic_variant) == "rs123",
  extract_revel_max(synthetic_variant) == 0.9,
  extract_selected_transcript_impact(synthetic_variant) == "high",
  extract_genes(synthetic_variant) == "GENE1"
)

# ANNOVAR numeric parsing -----------------------------------------------------

stopifnot(
  identical(
    parse_annovar_numeric(c(".", "0.1", "", NA_character_)),
    c(NA_real_, 0.1, NA_real_, NA_real_)
  ),
  identical(
    parse_annovar_numeric(c(NA_real_, 0.2)),
    c(NA_real_, 0.2)
  )
)

# ClinVar consolidation -------------------------------------------------------

stopifnot(
  collapse_clinvar_category("Pathogenic;Likely pathogenic") ==
    "Likely pathogenic/Pathogenic",
  collapse_clinvar_category("Benign") == "Likely benign/Benign",
  collapse_clinvar_category("drug response") == "Drug response",
  collapse_clinvar_category("risk factor") == "Risk factor",
  collapse_clinvar_category("Unmapped assertion") == "Other"
)

# Allele-count tests -----------------------------------------------------------

fisher_result <- safe_fisher_test(10, 200, 5, 200)
lrt_result <- safe_lrt(10, 200, 5, 200)

stopifnot(
  is.finite(fisher_result[["p_value"]]),
  is.finite(fisher_result[["odds_ratio"]]),
  is.finite(lrt_result[["statistic"]]),
  is.finite(lrt_result[["p_value"]]),
  abs(safe_lrt(10, 200, 10, 200)[["statistic"]]) < 1e-12,
  is.na(safe_fisher_test(201, 200, 5, 200)[["p_value"]])
)

message("All core-function tests passed.")
