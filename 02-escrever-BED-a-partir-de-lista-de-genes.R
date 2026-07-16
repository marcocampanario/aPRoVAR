#!/usr/bin/env Rscript

###################################### PACOTES ######################################

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
})

###################################### ARGUMENTOS ######################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Uso:\n")
  cat("Rscript genes_to_bed.R genes.txt output.bed\n")
  quit(status = 1)
}

input_file  <- args[1]
output_file <- args[2]

###################################### LEITURA ######################################

genes <- read_lines(input_file)
genes <- unique(genes)
genes <- genes[genes != ""]

cat("Genes lidos:", length(genes), "\n")

###################################### BIOMART ######################################

cat("Conectando ao Ensembl...\n")

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

###################################### QUERY ######################################

cat("Buscando coordenadas...\n")

df <- getBM(
  attributes = c(
    "chromosome_name",
    "start_position",
    "end_position",
    "hgnc_symbol"
  ),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

###################################### LIMPEZA ######################################

# Remove cromossomos não padrão
df <- df %>%
  filter(chromosome_name %in% c(1:22, "X", "Y"))

# Remove duplicados
df <- df %>%
  distinct(chromosome_name, start_position, end_position, hgnc_symbol)

# Ajuste BED (0-based)
bed <- df %>%
  mutate(
    chromosome_name = paste0("chr", chromosome_name),
    start = start_position - 1
  ) %>%
  select(chromosome_name, start, end_position, hgnc_symbol)

###################################### GENES NÃO ENCONTRADOS ######################################

genes_found <- unique(df$hgnc_symbol)
genes_missing <- setdiff(genes, genes_found)

cat("Genes encontrados:", length(genes_found), "\n")
cat("Genes NÃO encontrados:", length(genes_missing), "\n")

if (length(genes_missing) > 0) {
  write_lines(genes_missing, "genes_not_found.txt")
  cat("Lista salva em genes_not_found.txt\n")
}

###################################### OUTPUT ######################################

write.table(
  bed,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Arquivo BED salvo em:", output_file, "\n")
cat("✔ Finalizado!\n")
