# BAM processing and small-variant calling

This directory contains the DRAGEN-based workflow used to process the APROVAR-1010-WES FASTQ files. The script performs read mapping, coordinate sorting, duplicate marking, BAM indexing, HLA typing, and small-variant calling for each sample listed in a comma-separated sample sheet.

The workflow was developed for DRAGEN Bio-IT Platform v4.3.6 and the GRCh38 reference genome.

## Workflow

| Order | Script | Main purpose |
|---:|---|---|
| 1 | `01_DRAGEN_bam_processing_and_variant_call.sh` | Read the sample sheet, process each unique `RGSM`, align reads to GRCh38, generate an indexed BAM, and call SNVs/INDELs in VCF and compact GVCF formats. |

## Required inputs

1. A DRAGEN-compatible GRCh38 hash table:
   `hg38.alt_masked.cnv.graph.hla.rna-10-r4.0-1/`
2. A comma-separated FASTQ sample sheet, named `fastq_list.csv` by default.
3. FASTQ files referenced by the sample sheet.
4. A target-region BED file for the Twist Illumina Exome 2.0 Plus panel.
5. The DRAGEN machine-learning recalibration models distributed with DRAGEN v4.3.6.

## Software requirements

- DRAGEN Bio-IT Platform v4.3.6;
- Bash 4 or later, because the script uses an associative array;
- access to the reference hash table, target BED, temporary directory, and DRAGEN model files.

## Configuration

Edit the variables at the beginning of the script before execution:

```bash
DRAGEN_HASH_TABLE="/path/to/hg38.alt_masked.cnv.graph.hla.rna-10-r4.0-1/"
TEMPD="/path/to/tmp/"
MLMODEL="/opt/dragen/4.3.6/resources/ml_model/"

FASTQ_LIST="/path/to/fastq_list.csv"
OUPTUT_FOLDER="/path/to/results"
VCTARGETBED="/path/to/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed"
```

## Running the workflow

From this directory:

```bash
bash 01_DRAGEN_bam_processing_and_variant_call.sh
```

Samples are processed sequentially. An associative array tracks previously observed `RGSM` values, preventing duplicate processing when the same sample appears in more than one FASTQ-list row.

## DRAGEN configuration

The script enables:

- read mapping and alignment to GRCh38;
- coordinate sorting;
- duplicate marking;
- BAM output and BAM indexing;
- HLA typing, including class II loci;
- sex-aware processing of sex chromosomes;
- small-variant calling restricted to the exome target BED;
- compact GVCF output and standard VCF output;
- machine-learning-based variant recalibration.

Secondary and supplementary alignments are suppressed using:

```text
--Aligner.sec-aligns 0
--Aligner.supp-aligns 0
```

Supplementary alignments, if present, are treated as secondary. CNV and SV calling are disabled in this workflow.