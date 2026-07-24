library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(biomaRt)

# Load the project functions
source(file.path(project_dir, "R", "00_config.R"))
source(file.path(project_dir, "R", "01_functions.R"))

# Load the analysis table if it is not already available
if (!exists("analysis_table")) {
  analysis_ready <- read_rds_parallel(ANALYSIS_READY_RDS)
  analysis_table <- analysis_ready$analysis_table
}

# Extract all Ensembl gene IDs present in the dataset
ensembl_ids <- tibble(
  gene = purrr::map_chr(
    analysis_table$variants,
    extract_genes
  )
) |>
  tidyr::separate_rows(gene, sep = ";") |>
  dplyr::mutate(
    gene = trimws(gene),
    # Remove an optional Ensembl version suffix
    ensembl_gene_id = stringr::str_remove(
      gene,
      "\\.[0-9]+$"
    )
  ) |>
  dplyr::filter(
    stringr::str_detect(
      ensembl_gene_id,
      "^ENSG[0-9]+$"
    )
  ) |>
  dplyr::distinct(ensembl_gene_id)

nrow(ensembl_ids)

ensembl_mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Query in batches to avoid overly large requests
id_batches <- split(
  ensembl_ids$ensembl_gene_id,
  ceiling(seq_along(ensembl_ids$ensembl_gene_id) / 500)
)

gene_map_query <- purrr::map_dfr(
  id_batches,
  function(ids) {
    biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol"
      ),
      filters = "ensembl_gene_id",
      values = ids,
      mart = ensembl_mart
    )
  }
) |>
  dplyr::filter(
    !is.na(hgnc_symbol),
    nzchar(hgnc_symbol)
  ) |>
  dplyr::distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

gene_map <- ensembl_ids |>
  dplyr::left_join(
    gene_map_query,
    by = "ensembl_gene_id"
  ) |>
  dplyr::mutate(
    hgnc_symbol = dplyr::coalesce(
      hgnc_symbol,
      ""
    )
  ) |>
  dplyr::arrange(ensembl_gene_id)

c(
  ensembl_ids = nrow(gene_map),
  mapped = sum(nzchar(gene_map$hgnc_symbol)),
  unresolved = sum(!nzchar(gene_map$hgnc_symbol))
)

readr::write_tsv(
  gene_map,
  GENE_SYMBOL_MAP_TSV,
  na = ""
)
