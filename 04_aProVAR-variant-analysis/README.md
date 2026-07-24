# aPRoVAR variant annotation and descriptive analysis

This directory contains the modular, manuscript-oriented version of the original Campanário & Janke *et al.* (2026) analyses. The workflow uses one script per analytical block, and produces explicit checkpoints so that downstream modules can be rerun without repeating the full JSON/Nirvana parsing stage.

## Workflow

| Order | Script | Main purpose |
|---:|---|---|
| 1 | `scripts/01_prepare_variant_data.R` | Read the parsed multisample object, extract VCF metrics, split multiallelic sites, join ABRaOM annotations, and build the presence matrix. |
| 2 | `scripts/02_classify_variants.R` | Assign Known, Absent, or Novel status; calculate exact genotype-based allele counts; classify singleton, doubleton, rare, and polymorphic sites; retain `F_MISSING = 0`. |
| 3 | `scripts/03_variant_discovery.R` | Calculate cumulative variant discovery and aPRoVAR variant composition. |
| 4 | `scripts/04_clinvar_analysis.R` | Consolidate ClinVar germline assertions, identify P/LP variants, count genes, and generate word clouds. |
| 5 | `scripts/05_pharmacogenetic_analysis.R` | Match ClinPGx relationships and summarize the nine priority pharmacogenes. |
| 6 | `scripts/06_plof_analysis.R` | Identify strict MANE Select RefSeq high-impact pLOF variants, pD variants with REVEL score ≥ 0.773, and build ClinVar/pLOF/pD tables. |
| 7 | `scripts/07_frequency_spectrum.R` | Generate the allele-frequency spectrum overall and by variant origin. |
| 8 | `scripts/08_article_figures.R` | Reproduce all Campanário & Janke *et al.* (2026) figures. |
| 9 | `scripts/09_quality_control_checks.R` | Validate cross-module row counts, classifications, checkpoints, and output provenance. |

Shared paths and constants are defined in `R/00_config.R`; reusable functions are defined in `R/01_functions.R`.

## Required inputs

Place the following files in `data/raw/`, or override their paths with environment variables:

1. `multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_header.rds`
2. `multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_inicial.rds.gz`
3. `multi_vcf_XY_rerun_annovar.hg38_multianno_colunas-ok.txt`
4. `relationships.tsv` from the ClinPGx relationship export
5. `ensembl_to_hgnc.tsv` from Ensembl v116 BioMaRt query (Optional*)

`*` An optional `ensembl_to_hgnc.tsv` file may be supplied with columns `ensembl_gene_id` and `hgnc_symbol`. Our workflow deliberately avoids live biomaRt requests so that gene-symbol resolution remains versionable and reproducible.

The ANNOVAR table must contain `Chr`, `Start`, `Ref`, `Alt`, `abraom_freq`, and `AF`. The parsed variant table must contain `chromosome`, `position`, `refAllele`, `altAlleles`, `samples`, `variants`, and `vcfInfo`.

## Installation

Install the required R packages once:

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "openxlsx", "purrr", "readr",
  "scales", "stringr", "tibble", "tidyr", "writexl"
))

# Optional figure packages
install.packages(c("ggVennDiagram", "ggwordcloud"))
```

`pigz` is optional. When available, the workflow uses it for parallel RDS decompression and compression; otherwise it falls back to base R.

## Running the workflow

From this project directory:

```bash
Rscript run_all.R
```

Each module can also be executed independently after its prerequisite checkpoint exists:

```bash
Rscript scripts/04_clinvar_analysis.R
```

Outputs are written to:

- `data/processed/`: reusable RDS checkpoints;
- `results/tables/`: TSV and Excel source tables;
- `results/figures/`: publication-resolution PNG/TIFF figures;
- `results/reports/`: metrics, diagnostic terms, provenance, and `sessionInfo()`.

## Path overrides

The workflow does not use `setwd()`. Configure paths from the shell when files are stored elsewhere:

```bash
export APROVAR_DATA_DIR=/path/to/aprovar/data
export APROVAR_RESULTS_DIR=/path/to/aprovar/results
export APROVAR_ABRAOM_ANNOVAR_TSV=/path/to/multianno.txt
export APROVAR_CLINPGX_RELATIONSHIPS_TSV=/path/to/relationships.tsv
Rscript run_all.R
```

All supported variables are documented in `R/00_config.R`.

## Analytical definitions

- **Known:** observed with AF > 0 in gnomAD Genomes, gnomAD Exomes, 1000 Genomes, or ABraOM. Both gnomAD resources contribute to population-presence classification; when one gnomAD frequency is required for reporting, Genomes is preferred and Exomes is used as the fallback.
- **Absent:** absent from those population-frequency resources but assigned a dbSNP identifier.
- **Novel:** absent from the population-frequency resources and without a dbSNP identifier.
- **Singleton:** alternate allele count (AC) = 1.
- **Doubleton:** AC = 2; homozygous doubletons are retained as a traceable subclass.
- **Rare:** AC > 2 and AF < 1%.
- **Polymorphic:** AF ≥ 1%.
- **pLOF:** predicted loss-of-function, based on high-impact consequence on a canonical MANE Select RefSeq transcript.
- **pD:** predicted damaging, based on REVEL score ≥ 0.773.


Frequency classes are derived from parsed genotypes. This handles haploid calls and chromosome-specific allele numbers more safely.

## Important count provenance

`AF_EXCL` is preferred over `AF` when present because it represents the phenotype-adjusted aPRoVAR frequency. If the corresponding adjusted allele number is unavailable, the population-comparison module uses `APROVAR_LOCAL_ADJUSTED_AN` and explicitly flags `local_counts_inferred = TRUE`.

ANNOVAR removes VCF anchor alleles and may simplify or reposition INDEL and complex-allele representations. Consequently, direct `CHR-POS-REF-ALT` matching is incomplete for this input. The workflow aligns the ANNOVAR table by source row index only after confirming equal row counts and complete chromosome agreement. These parameters are specific for our data set, based on prior comparison of records between aPRoVAR and ABraOM. The measured alignment metrics are written to `results/reports/01_annovar_alignment_metrics.tsv`; the workflow stops if any safeguard fails.

## Validation

Run the lightweight helper-function tests with:

```bash
Rscript tests/test_core_functions.R
```

The full workflow ends with structural quality-control checks. Because the raw aPRoVAR inputs are not included here, syntax and unit checks do not replace a complete run on the source dataset.
