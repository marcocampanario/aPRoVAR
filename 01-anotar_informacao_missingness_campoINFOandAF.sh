#!/bin/bash

### Precisamos primeiro ativar o env com os plugins do bcftools
# conda activate bcftools_env

########################################
# Objetivo:
# Adicionar INFO/F_MISSING e INFO/AF ao VCF
########################################

# Checagem de argumentos
if [ -z "$1" ]; then
  echo "Uso: $0 <arquivo.vcf.gz>"
  echo "Exemplo: $0 input.vcf.gz"
  exit 1
fi

INPUT="$1"

# Extrai nome base
BASENAME=$(basename "$INPUT")
BASENAME=${BASENAME%.vcf.gz}
BASENAME=${BASENAME%.vcf}

OUTPUT="${BASENAME}_withMissingness-AF.vcf.gz"

echo "-----------------------------"
echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo "-----------------------------"

########################################
# Adiciona INFO/F_MISSING e INFO/AF
########################################

bcftools +fill-tags \
  "$INPUT" \
  -Oz -o "$OUTPUT" \
  -- -t F_MISSING,AF

########################################
# Indexação
########################################

bcftools index -t "$OUTPUT"

echo "✔ Finalizado com sucesso!"
echo "VCF anotado: $OUTPUT"
