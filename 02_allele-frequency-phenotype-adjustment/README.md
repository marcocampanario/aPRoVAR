# Phenotype-adjusted allele-frequency estimation

This directory contains the workflow used to minimize phenotype-driven bias in aPRoVAR allele-frequency estimates. Because the APROVAR-1010-WES cohort includes individuals recruited through studies of COVID-19, sepsis, and breast cancer, frequencies in genes associated with each phenotype are recalculated after excluding participants recruited because of the corresponding condition.

The workflow preserves the full-cohort allele frequency in `INFO/AF` and adds the phenotype-excluded estimate as `INFO/AF_EXCL`. For variants in phenotype-associated genes, `AF_EXCL` is the frequency reported by default in downstream aPRoVAR analyses.

## Workflow

| Order | Script | Main purpose |
|---:|---|---|
| 1 | `01-anotar_informacao_missingness_campoINFOandAF.sh` | Add full-cohort `F_MISSING` and `AF` INFO tags to the multisample VCF. |
| 2 | `02-escrever-BED-a-partir-de-lista-de-genes.R` | Query Ensembl BioMart for gene coordinates and create a BED file for each phenotype-associated gene list. |
| 3 | `03-recalcular_AF_vies_fenotipo.sh` | Restrict the VCF to each phenotype-specific BED, exclude the matching sample set, and recalculate `AC`, `AN`, and `AF`. |
| 4 | `04-anotar_vcf_original_com_AFs_recalculadas.sh` | Add the recalculated phenotype-excluded frequencies to the original VCF as `INFO/AF_EXCL`. |

For each phenotype, the procedure is:

1. calculate allele frequencies in the complete cohort;
2. identify variants located within the corresponding phenotype-associated gene intervals;
3. exclude all individuals recruited because of that phenotype;
4. recalculate `AC`, `AN`, and `AF` among the remaining individuals;
5. annotate the recalculated frequency back into the complete multisample VCF as `AF_EXCL`.

Variants outside the selected phenotype-associated regions retain frequencies calculated from all available individuals.

## Phenotype-associated gene sets

Distinct evidence sources were used for the three phenotype groups:

| Phenotype | Number of genes | Selection basis |
|---|---:|---|
| Breast cancer | 26 | Predisposition genes reported in the WHO Classification of Breast Tumours (2019) and genes included in commercial hereditary-cancer panels frequently used in Paraná. |
| Critical or hospitalized COVID-19 | 44 | Loci reported by Pathak *et al.* (2022) and the Brazilian BRACOVID study by Pereira *et al.* (2022). |
| Sepsis susceptibility | 8 | Genes selected from the genomic endotype study by Scicluna *et al.* (2017). |

The exact gene lists used in the study are reported in Additional file 2, Table S2 of Campanário & Janke *et al.* (2026).

## Required inputs

1. A bgzip-compressed, indexed multisample VCF.
2. Plain-text gene list per phenotype, with one HGNC symbol per line.
3. Text file per phenotype containing the sample identifiers to exclude, in the format accepted by `bcftools view -S`.
4. A tab-delimited `config.txt` connecting each phenotype to its BED and exclusion list.

Example `config.txt`:

```text
COVID19	covid19_genes.bed	covid19_samples.txt
SEPSIS	sepsis_genes.bed	sepsis_samples.txt
BREAST_CANCER	breast_cancer_genes.bed	breast_cancer_samples.txt
```

Do not include a header unless the scripts are adapted to skip it. The phenotype label is also used to construct output filenames.

## Software requirements

- Bash;
- `bcftools`, including the `+fill-tags` plugin;
- `bgzip` and `tabix`;
- R;
- R packages `biomaRt`, `dplyr`, and `readr`;
- internet access during the BioMart coordinate query.

The scripts assume that the input VCF and all intermediate VCFs are bgzip-compressed and tabix-indexed.

## Interpretation of `AF` and `AF_EXCL`

- `AF`: allele frequency estimated from the complete APROVAR-1010-WES cohort.
- `AF_EXCL`: allele frequency estimated after removing individuals recruited because of the phenotype associated with that genomic region.

## References for gene selection

1. WHO Classification of Tumours Editorial Board. *Breast Tumours*. 5th ed. International Agency for Research on Cancer; 2019.
2. Pathak GA, Karjalainen J, Stevens C, *et al.* A first update on mapping the human genetic architecture of COVID-19. *Nature*. 2022;608:E1–E10. [https://doi.org/10.1038/s41586-022-04826-7](https://doi.org/10.1038/s41586-022-04826-7)
3. Pereira AC, Bes TM, Velho M, *et al.* Genetic risk factors and COVID-19 severity in Brazil: results from BRACOVID study. *Human Molecular Genetics*. 2022;31:3021–3031. [https://doi.org/10.1093/hmg/ddac045](https://doi.org/10.1093/hmg/ddac045)
4. Scicluna BP, van Vught LA, Zwinderman AH, *et al.* Classification of patients with sepsis according to blood genomic endotype: a prospective cohort study. *Lancet Respiratory Medicine*. 2017;5:816–826. [https://doi.org/10.1016/S2213-2600(17)30294-1](https://doi.org/10.1016/S2213-2600(17)30294-1)

