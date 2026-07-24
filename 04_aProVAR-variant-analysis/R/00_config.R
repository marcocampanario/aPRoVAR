##### 00_config.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Define input files, output directories, analysis constants, and figure settings
# used throughout the aPRoVAR variant annotation workflow.

# Resolve the project directory from an explicit environment variable whenever
# possible. run_all.R sets this variable before sourcing the analysis modules.
PROJECT_DIR <- Sys.getenv("APROVAR_PROJECT_DIR", unset = "")

if (!nzchar(PROJECT_DIR)) {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) == 1L) {
    current_file <- normalizePath(
      sub("^--file=", "", file_arg),
      mustWork = FALSE
    )
    PROJECT_DIR <- dirname(dirname(current_file))
  } else {
    PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)
  }
}

PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = FALSE)

# Directory layout ------------------------------------------------------------

DATA_DIR <- Sys.getenv(
  "APROVAR_DATA_DIR",
  unset = file.path(PROJECT_DIR, "data")
)

RAW_DATA_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DATA_DIR <- file.path(DATA_DIR, "processed")

RESULTS_DIR <- Sys.getenv(
  "APROVAR_RESULTS_DIR",
  unset = file.path(PROJECT_DIR, "results")
)

TABLES_DIR <- file.path(RESULTS_DIR, "tables")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
REPORTS_DIR <- file.path(RESULTS_DIR, "reports")

# Primary inputs --------------------------------------------------------------

# The defaults reproduce the final May 2026 checkpoint used in the original
# monolithic script. Override any path with the corresponding environment
# variable instead of editing analysis code.
HEADER_RDS <- Sys.getenv(
  "APROVAR_HEADER_RDS",
  unset = file.path(
    RAW_DATA_DIR,
    "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_header.rds"
  )
)

VARIANT_RDS_GZ <- Sys.getenv(
  "APROVAR_VARIANT_RDS_GZ",
  unset = file.path(
    RAW_DATA_DIR,
    "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_inicial.rds.gz"
  )
)

ABRAOM_ANNOVAR_TSV <- Sys.getenv(
  "APROVAR_ABRAOM_ANNOVAR_TSV",
  unset = file.path(
    RAW_DATA_DIR,
    "multi_vcf_XY_rerun_annovar.hg38_multianno_colunas-ok.txt"
  )
)

CLINPGX_RELATIONSHIPS_TSV <- Sys.getenv(
  "APROVAR_CLINPGX_RELATIONSHIPS_TSV",
  unset = file.path(RAW_DATA_DIR, "relationships.tsv")
)

# Optional two-column table with Ensembl gene IDs and HGNC symbols. It avoids a
# live biomaRt query and therefore keeps the workflow reproducible offline.
GENE_SYMBOL_MAP_TSV <- Sys.getenv(
  "APROVAR_GENE_SYMBOL_MAP_TSV",
  unset = file.path(RAW_DATA_DIR, "ensembl_to_hgnc.tsv")
)

# Intermediate checkpoints ----------------------------------------------------

PREPARED_DATA_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "01_prepared_variant_data.rds.gz"
)

ANALYSIS_READY_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "02_analysis_ready_variants.rds.gz"
)

DISCOVERY_RESULTS_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "03_variant_discovery_results.rds.gz"
)

CLINVAR_RESULTS_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "04_clinvar_results.rds.gz"
)

PHARMACOGENETIC_RESULTS_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "05_pharmacogenetic_results.rds.gz"
)

PLOF_RESULTS_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "06_plof_results.rds.gz"
)

POPULATION_COMPARISON_RDS <- file.path(
  PROCESSED_DATA_DIR,
  "07_population_frequency_comparison.rds.gz"
)

# Cohort and classification constants ----------------------------------------

EXPECTED_SAMPLE_COUNT <- as.integer(
  Sys.getenv("APROVAR_EXPECTED_SAMPLE_COUNT", unset = "1010")
)

POLYMORPHIC_AF_THRESHOLD <- as.numeric(
  Sys.getenv("APROVAR_POLYMORPHIC_AF_THRESHOLD", unset = "0.01")
)

REVEL_HIGH_THRESHOLD <- as.numeric(
  Sys.getenv("APROVAR_REVEL_HIGH_THRESHOLD", unset = "0.773")
)

RANDOM_SEED <- as.integer(
  Sys.getenv("APROVAR_RANDOM_SEED", unset = "42")
)

# The existing ANNOVAR table was generated from the same split VCF and retains
# its row order, but ANNOVAR rewrites many INDEL and complex-allele keys. Row-
# index alignment is allowed only after these minimum same-row match rates pass.
ANNOVAR_MIN_EXACT_ROW_MATCH <- as.numeric(
  Sys.getenv("APROVAR_ANNOVAR_MIN_EXACT_ROW_MATCH", unset = "0.90")
)

ANNOVAR_MIN_STANDARD_ROW_MATCH <- as.numeric(
  Sys.getenv("APROVAR_ANNOVAR_MIN_STANDARD_ROW_MATCH", unset = "0.98")
)

PIGZ_THREADS <- as.integer(
  Sys.getenv("APROVAR_PIGZ_THREADS", unset = "16")
)

# ABRaOM allele counts are not present in the ANNOVAR frequency-only input.
# Set APROVAR_ABRAOM_AN only when a defensible allele number is available for
# the analyzed locus. Fisher tests are left missing when it is not supplied.
ABRAOM_ALLELE_NUMBER <- suppressWarnings(
  as.numeric(Sys.getenv("APROVAR_ABRAOM_AN", unset = NA_character_))
)

# When AF_EXCL is used but its matching allele number is unavailable, the
# statistical comparison script can infer integer counts from this denominator.
# Such rows are explicitly flagged as inferred in the output.
LOCAL_ADJUSTED_ALLELE_NUMBER <- suppressWarnings(
  as.numeric(
    Sys.getenv(
      "APROVAR_LOCAL_ADJUSTED_AN",
      unset = as.character(2L * EXPECTED_SAMPLE_COUNT)
    )
  )
)

# Pharmacogenes prioritized in the aPRoVAR manuscript.
PRIORITY_PHARMACOGENES <- c(
  "CYP2C9",
  "CYP2D6",
  "CYP3A5",
  "VKORC1",
  "SLCO1B1",
  "TPMT",
  "NUDT15",
  "DPYD",
  "UGT1A1"
)

# External population fields used by Illumina Nirvana/gnomAD annotations.
EXTERNAL_POPULATION_FIELDS <- data.frame(
  population = c(
    "gnomAD NFE",
    "gnomAD AFR",
    "gnomAD EAS",
    "gnomAD AMR"
  ),
  af_field = c("nfeAf", "afrAf", "easAf", "amrAf"),
  ac_field = c("nfeAc", "afrAc", "easAc", "amrAc"),
  an_field = c("nfeAn", "afrAn", "easAn", "amrAn"),
  stringsAsFactors = FALSE
)

# Figure settings -------------------------------------------------------------

FIGURE_WIDTH <- 9
FIGURE_HEIGHT <- 5
FIGURE_DPI <- 600

ORIGIN_COLORS <- c(
  "Known" = "#0072B2",
  "Absent" = "#E69F00",
  "Novel" = "#009E73"
)

CLINVAR_COLORS <- c(
  "Likely pathogenic/Pathogenic" = "#D55E00",
  "Likely benign/Benign" = "#0072B2",
  "Uncertain significance" = "#56B4E9",
  "Conflicting classifications of pathogenicity" = "#E69F00",
  "Drug response" = "#009E73",
  "Affects a non-disease phenotype" = "#CC79A7",
  "Protective" = "#F0E442",
  "Low penetrance for Mendelian diseases" = "#999999",
  "Not provided" = "#444444",
  "GWAS hits" = "#E69F00",
  "Risk factor" = "#56B4E9",
  "Other" = "grey50"
)
