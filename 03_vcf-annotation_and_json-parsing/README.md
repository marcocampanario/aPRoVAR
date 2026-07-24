# VCF annotation and Nirvana JSON parsing

Scripts by Tiago Minuzzi Freire da Fontoura Gomes ![blob-pkmn-ghastly.png](blob-pkmn-ghastly.png)

This directory contains the workflow used to parse the large JSON output produced by Illumina Nirvana and convert its position-, variant-, sample-, and gene-level annotations into reusable R objects for downstream aPRoVAR analyses.

The input to this block is the phenotype-adjusted multisample VCF generated in `../02_allele-frequency-phenotype-adjustment/` after annotation with Nirvana against GRCh38 resources.

## Workflow

| Order | Script or step | Main purpose |
|---:|---|---|
| 1 | `01_run_Nirvana.sh` | Annotate AF-adjusted multisample VCF and produce a compressed JSON file containing header, position, variant, sample, transcript, population, and gene annotations. |
| 2 | `02_read_json.R` | Split the Nirvana JSON into header, position, and gene sections; parse position and gene records in parallel; construct R tibbles; and save reusable compressed RDS checkpoints. |

## Required input

The R script expects the following file in its working directory:

```text
multi_gvcf_FINAL_FMISSING_AF_contextual.json.gz
```

This must be a line-oriented Nirvana JSON output with:

- the metadata header on the first line;
- a `positions` array;
- an optional `genes` array separated by the literal line `],"genes":[`;
- the closing JSON line at the end of the file.

## Software requirements

- R;
- `pigz`;
- R packages:
  - `jsonlite`;
  - `tidyverse`;
  - `data.table`;
  - `purrr`;
  - `furrr`;
  - `future`;
  - `progressr`;
  - `dplyr`;
  - `tidyr`.

Install the required R packages once:

```r
install.packages(c(
  "jsonlite", "tidyverse", "data.table", "purrr",
  "furrr", "future", "progressr", "dplyr", "tidyr"
))
```

`pigz` must be available on the shell `PATH`, because the script uses pipe connections for parallel compression and decompression.

## Resource configuration

The deposited script uses:

```r
plan(multisession, workers = 35)
```

and:

```text
pigz -p 16
```

These values reflect the original computational environment. Reduce the number of R workers and pigz threads if the available machine has fewer CPU cores or insufficient memory. Nirvana JSON files may be very large, and the workflow temporarily stores raw JSON strings and parsed R objects simultaneously.

## Running the parser

Place the Nirvana JSON file in the working directory and run:

```bash
Rscript 02_read_json.R
```

The script executes the following stages:

1. reads all lines of the JSON file;
2. separates the header, `positions`, and optional `genes` sections;
3. removes trailing commas from individual JSON records;
4. saves the unparsed character vectors as an intermediate compressed RDS;
5. parses position records in parallel with `furrr::future_map()` and `jsonlite::fromJSON()`;
6. parses gene records in parallel;
7. extracts the Nirvana header;
8. constructs `tbl_inicial`, containing site-, sample-, and variant-level fields;
9. constructs `tbl_genes`, containing gene-level identifiers and database annotations;
10. saves the final R objects using parallel pigz compression.

## Output objects

### Intermediate files

```text
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes.rds.gz
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes_positions_list.rds.gz
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_genes_list.rds.gz
```

These files contain, respectively:

- the extracted header and unparsed position/gene JSON strings;
- the parsed position-record list;
- the parsed gene-record list.

### Final files

```text
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_header.rds
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_inicial.rds.gz
multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_genes.rds.gz
```

`header.rds` contains the parsed Nirvana metadata header.

`tbl_inicial.rds.gz` contains one record per top-level Nirvana position with:

| Field | Content |
|---|---|
| `chromosome` | Chromosome name |
| `position` | Genomic position |
| `refAllele` | Reference allele |
| `altAlleles` | Alternate-allele list |
| `vcfInfo` | VCF INFO annotations, including `AF`, `AF_EXCL`, and `F_MISSING` when available |
| `quality` | Site quality |
| `filters` | VCF filter status |
| `cytogeneticBand` | Cytogenetic-band annotation |
| `samples` | Nested sample and genotype information |
| `variants` | Nested allele- and transcript-level annotations |

`tbl_genes.rds.gz` contains:

| Field | Content |
|---|---|
| `name` | Gene symbol |
| `ensemblGeneId` | Ensembl gene identifier |
| `hgncId` | HGNC identifier |
| `ncbiGeneId` | NCBI Gene identifier |
| `gnomAD` | Nested gnomAD gene annotation |
| `omim` | Nested OMIM annotation |
| `clingenGeneValidity` | ClinGen gene-validity annotation |
| `clingenDosageSensitivityMap` | ClinGen dosage-sensitivity annotation |
| `cosmic` | COSMIC gene annotation |

The `header.rds` and `tbl_inicial.rds.gz` files are required inputs for the modular analytical R package in `../04_aProVAR-variant-analysis/`.
