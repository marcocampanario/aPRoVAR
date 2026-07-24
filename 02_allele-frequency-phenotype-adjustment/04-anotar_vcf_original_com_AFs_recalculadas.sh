#!/bin/bash

########################################
# Uso:
# bash anotar_AF_vies_fenotipo.sh input.vcf.gz config.txt
########################################

VCF="$1"
CONFIG="$2"

OUT="multi_gvcf_FINAL_FMISSING_AF_contextual.vcf.gz"

if [ -z "$VCF" ] || [ -z "$CONFIG" ]; then
  echo "Uso: $0 <input.vcf.gz> <config.txt>"
  exit 1
fi

if [ ! -f "$VCF" ]; then
  echo "Erro: VCF não encontrado"
  exit 1
fi

########################################
# Pipeline
########################################

PIPE="$VCF"

while IFS=$'\t' read -r PHENO BED SAMPLES; do

  echo "-----------------------------"
  echo "Anotando: $PHENO"
  echo "-----------------------------"

  AF_TXT="${PHENO}_AF_recalc.txt"
  AF_GZ="${AF_TXT}.gz"
  HEADER="header_${PHENO}.txt"
  TAG="AF_EXCL"

  ########################################
  # Checagem
  ########################################

  if [ ! -f "$AF_TXT" ]; then
    echo "Erro: $AF_TXT não encontrado"
    exit 1
  fi

  ########################################
  # Compressão + indexação
  ########################################

  bgzip -c "$AF_TXT" > "$AF_GZ"
  tabix -s 1 -b 2 -e 2 "$AF_GZ"

  ########################################
  # Header
  ########################################

  echo "##INFO=<ID=${TAG},Number=A,Type=Float,Description=\"AF recalculada excluindo ${PHENO}\">" > "$HEADER"

  ########################################
  # Annotate
  ########################################

  bcftools annotate \
    -a "$AF_GZ" \
    -c CHROM,POS,REF,ALT,INFO/${TAG} \
    -h "$HEADER" \
    "$PIPE" \
    -Oz -o tmp_${PHENO}.vcf.gz

  PIPE="tmp_${PHENO}.vcf.gz"

done < "$CONFIG"

########################################
# Final
########################################

mv "$PIPE" "$OUT"
bcftools index -t "$OUT"

echo "============================="
echo "✔ VCF final:"
echo "$OUT"
