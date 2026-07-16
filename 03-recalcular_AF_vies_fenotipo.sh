#!/bin/bash

########################################
# Uso:
# ./recalculate_af.sh input.vcf.gz config.txt
########################################

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Uso: $0 <input.vcf.gz> <config.txt>"
  exit 1
fi

VCF="$1"
CONFIG="$2"

########################################
# Checagem
########################################

if [ ! -f "$VCF" ]; then
  echo "Erro: VCF não encontrado"
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "Erro: config.txt não encontrado"
  exit 1
fi

########################################
# Loop pelos fenótipos
########################################

while IFS=$'\t' read -r PHENO BED SAMPLES; do

  echo "-----------------------------"
  echo "Processando: $PHENO"
  echo "BED: $BED"
  echo "Samples a excluir: $SAMPLES"
  echo "-----------------------------"

  OUT_VCF="${PHENO}_AF_recalc.vcf.gz"
  OUT_TXT="${PHENO}_AF_recalc.txt"

  ########################################
  # Recalcular AF
  ########################################

  bcftools view \
    -R "$BED" \
    -S ^"$SAMPLES" \
    -Ou "$VCF" | \
  bcftools +fill-tags \
    -Oz -o "$OUT_VCF" \
    -- -t AC,AN,AF

  ########################################
  # Indexar
  ########################################

  bcftools index -t "$OUT_VCF"

  ########################################
  # Extrair tabela (opcional)
  ########################################

  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' \
    "$OUT_VCF" > "$OUT_TXT"

  echo "✔ Gerado:"
  echo "  - $OUT_VCF"
  echo "  - $OUT_TXT"

done < "$CONFIG"

echo "============================="
echo "✔ Todos os fenótipos processados!"
