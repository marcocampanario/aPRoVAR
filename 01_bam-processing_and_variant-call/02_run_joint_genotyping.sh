##### 02_run_joint_genotyping.sh #####
# Tiago Minuzzi Freire da Fontoura Gomes @ Fiocruz Paraná, 2026 Jul.
# Generation of multisample VCF via joint genotyping on DRAGEN,
# based on the list of individual gVCFs.

#!/bin/bash

DRAGEN="/opt/dragen/4.3.13/bin/dragen"
REFERENCE="/staging/human/reference/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa"
OUTDIR="/path/to/output-dir/"
PREFIX="output_file_prefix"
BED="/CSPP-DRAGEN-01/minuzzi/labioinfo/BEDs/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed"
LIST="/path/to/individual-GVCFs-list.txt"

$DRAGEN --enable-gvcf-genotyper-iterative true \
	--ht-reference $REFERENCE \
	--output-directory  $OUTDIR \
	--output-file-prefix $PREFIX \
	--gg-regions-bed $BED \
	--variant-list $LIST \
	--gg-discard-ac-zero true \
	--gg-vc-filter true \
	--gg-remove-nonref
