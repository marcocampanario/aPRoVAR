# aPRoVAR: Arquivo Paranaense Online de Variantes Genéticas

This repository contains the bioinformatics workflows used to generate, process, annotate, and analyze the APROVAR-1010-WES dataset described by Campanário & Janke *et al.* (2026).

aPRoVAR is a regional genomic resource comprising germline whole-exome sequencing data from 1,010 individuals from Paraná, Southern Brazil. The repository is organized into four consecutive analytical blocks, from raw sequencing reads to the manuscript-oriented descriptive analyses.

## Analytical workflow

| Block | Directory | Main purpose |
|---:|---|---|
| 01 | [`01_bam-processing_and_variant-call/`](01_bam-processing_and_variant-call/) | Process FASTQ files with DRAGEN, generate indexed BAMs, call small variants in VCF and GVCF files, and generate multisample VCF. |
| 02 | [`02_allele-frequency-phenotype-adjustment/`](02_allele-frequency-phenotype-adjustment/) | Calculate full-cohort allele frequencies and recalculate frequencies in phenotype-associated genes after excluding the matching recruitment group. |
| 03 | [`03_vcf-annotation_and_json-parsing/`](03_vcf-annotation_and_json-parsing/) | Parse Nirvana JSON annotations and convert position-, sample-, variant-, transcript-, and gene-level data into reusable R objects. |
| 04 | [`04_aProVAR-variant-analysis/`](04_aProVAR-variant-analysis/) | Classify and summarize variants, evaluate ClinVar and pharmacogenetic annotations, identify pLOF/pD variants, and reproduce article figures and tables. |

```text
FASTQ
  |
  v
01 — DRAGEN mapping and small-variant calling
  |
  v
Multisample VCF
  |
  v
02 — Phenotype-adjusted allele-frequency estimation
  |
  v
VCF containing AF, F_MISSING, and AF_EXCL
  |
  v
03 — Nirvana annotation and JSON parsing
  |
  v
Parsed RDS objects with Nirvana JSON annotations
  |
  v
04 — Variant classification and descriptive analyses
  |
  v
Manuscript tables, figures, checkpoints, and QC reports
```

## Repository structure

```text
aPRoVAR/
├── 01_bam-processing_and_variant-call/
│   ├── 01_DRAGEN_bam_processing_and_variant_call.sh
│   └── README.md
├── 02_allele-frequency-phenotype-adjustment/
│   ├── 01-anotar_informacao_missingness_campoINFOandAF.sh
│   ├── 02-escrever-BED-a-partir-de-lista-de-genes.R
│   ├── 03-recalcular_AF_vies_fenotipo.sh
│   ├── 04-anotar_vcf_original_com_AFs_recalculadas.sh
│   └── README.md
├── 03_vcf-annotation_and_json-parsing/
│   ├── 01_run_Nirvana.R
│   ├── 02_read_json.R
│   └── README.md
├── 04_aProVAR-variant-analysis/
│   ├── R/
│   ├── scripts/
│   ├── tests/
│   ├── results/
│   └── README.md
├── LICENSE
└── README.md
```

## Block 01: BAM processing and variant calling

The first block processes each sample with DRAGEN Bio-IT Platform v4.3.6. Reads are aligned to GRCh38, coordinate-sorted, duplicate-marked, and written as indexed BAM files. Small variants are called within the Twist Illumina Exome 2.0 Plus target regions, producing per-sample VCFs and GVCFs, and a multisample VCF after joint genotyping.

CNV and SV calling were disabled. Sample sex is supplied through the FASTQ sample sheet.

See the [Block 01 README](01_bam-processing_and_variant-call/README.md) for the required sample-sheet format, DRAGEN resources, configuration variables, and outputs.

## Block 02: Phenotype-adjusted allele frequencies

The APROVAR-1010-WES cohort contains individuals recruited through studies of COVID-19, sepsis, and breast cancer. To minimize phenotype-driven frequency bias, variants in genes associated with each condition are recalculated after excluding individuals recruited because of the matching phenotype.

The workflow:

1. adds full-cohort `AF` and `F_MISSING`;
2. converts phenotype-specific gene lists into BED intervals;
3. recalculates `AC`, `AN`, and `AF` after phenotype-specific sample exclusions;
4. annotates the recalculated value into the original VCF as `AF_EXCL`.

See the [Block 02 README](02_allele-frequency-phenotype-adjustment/README.md) for the `config.txt` format, commands, phenotype-gene evidence, references, and interpretation of `AF_EXCL`.

## Block 03: VCF annotation and JSON parsing

The phenotype-adjusted multisample VCF is annotated with Illumina Nirvana. The included R parser separates the large Nirvana JSON into header, position, and gene sections; parses position and gene records in parallel; and constructs compressed RDS objects for downstream analyses.

The current repository snapshot includes the JSON parser but not the command or script used to run Nirvana. A compatible Nirvana JSON file is therefore required to reproduce this block.

See the [Block 03 README](03_vcf-annotation_and_json-parsing/README.md) for the expected JSON structure, R dependencies, memory and parallelization settings, output objects, and validation checks.

## Block 04: Variant annotation and descriptive analysis

The final block is a modular R workflow that:

- prepares and validates the parsed variant data;
- assigns Known, Absent, and Novel status;
- calculates exact genotype-based allele counts and frequency classes;
- summarizes cumulative variant discovery;
- consolidates ClinVar germline assertions;
- evaluates ClinPGx relationships and priority pharmacogenes;
- identifies strict MANE Select pLOF variants and REVEL-supported pD variants;
- reproduces manuscript figures and tables;
- runs cross-module quality-control checks.

See the [Block 04 README](04_aProVAR-variant-analysis/README.md) for installation, raw-data requirements, path overrides, analytical definitions, validation, and module-level execution.

## Companion repository

Code used to develop the aPRoVAR web interface is available at:

[https://github.com/omatheuspimenta/aPRoVAR](https://github.com/omatheuspimenta/aPRoVAR)

## License

See [`LICENSE`](LICENSE) for the terms governing use of the code in this repository.