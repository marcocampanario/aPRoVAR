##### 02_read_json.R #####
# Alysson Henrique Urbanski & Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Parallel parsing and processing of Nirvana JSON annotations.
# Parses position- and gene-level annotations, converts them into
# structured R objects, and saves the results using parallel             
# compression.             

# Required libraries
library(jsonlite)
library(tidyverse)
library(data.table)
library(purrr) # Useful for removing NULL values at the end
library(furrr)
library(future)

# Version 2.0, now robust to files without the "genes" section
read_anno_json_ultra_fast_v2 <- function(file.path) {
  
  print("Reading all file lines...")
  all_lines <- readLines(file.path)
  
  print("Extracting sections...")
  # 1. Attempt to locate the separator line
  gene_section_idx <- which(trimws(all_lines) == '],"genes":[')
  
  # 2. CONDITIONAL LOGIC: Check whether the "genes" section was found
  if (length(gene_section_idx) > 0) {
    # CASE A: The 'genes' section EXISTS (previous behavior)
    print("Section 'genes' found. Splitting the file...")
    position_lines <- all_lines[2:(gene_section_idx - 1)]
    gene_lines <- all_lines[(gene_section_idx + 1):(length(all_lines) - 1)]
    
  } else {
    # CASE B: The 'genes' section DOES NOT EXIST
    print("WARNING: Section 'genes' not found. Assuming file contains just positions.")
    # Positions span from line 2 to the penultimate line of the file.
    position_lines <- all_lines[2:(length(all_lines) - 1)]
    # An empty gene vector is created.
    gene_lines <- character(0)
  }
  
  # From this point onward, the remaining code works for both cases
  
  header_line <- all_lines[1]
  header <- substr(header_line, 11, nchar(header_line) - 14)
  
  # Free the memory used by the original large vector
  rm(all_lines)
  gc()
  
  print("Cleaning strings with data.table...")
  
  # Process 'positions'
  dt_pos <- data.table(lines = position_lines)
  dt_pos[endsWith(lines, ","), lines := substr(lines, 1, nchar(lines) - 1)]
  positions <- dt_pos$lines
  rm(dt_pos, position_lines)
  gc()
  
  # Process 'genes' (if 'gene_lines' is empty, 'genes' will also remain empty, as expected)
  dt_gen <- data.table(lines = gene_lines)
  dt_gen[endsWith(lines, ","), lines := substr(lines, 1, nchar(lines) - 1)]
  genes <- dt_gen$lines
  rm(dt_gen, gene_lines)
  gc()
  
  print("Processing concluded!")
  print(paste('Number of positions:', length(positions)))
  print(paste('Number of genes:', length(genes)))
  
  return(list(header = header, positions = positions, genes = genes))
}

# --- Usage ---
json.file <- "multi_gvcf_FINAL_FMISSING_AF_contextual.json.gz"

tempo_execucao_v2 <- system.time({
  resultados_ultra_v2 <- read_anno_json_ultra_fast_v2(json.file)
}, gcFirst = T)

cat("Reading of json file (version ultra_fast_v2) completed. Execution time:", tempo_execucao_v2)


# -------------------------------------------------------------------------

# Define the output filename
output_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes.rds.gz"

# Create a pipe connection to pigz
# The command runs pigz with 16 threads and redirects the output to output_file
# 'wb' = write binary mode
con <- pipe(sprintf("pigz -p 16 > %s", output_file), "wb")

# Use saveRDS with this connection and DISABLE R internal compression
# R writes the raw data to the pipe as quickly as possible, while pigz handles the compression.
print("Saving object with parallel compression via pigz...")
saveRDS(object = resultados_ultra_v2, file = con, compress = FALSE)

# Closing the connection is ESSENTIAL to ensure that all data are written correctly
close(con)
print("Save successful!")

input_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes.rds.gz"

# The command decompresses the file with pigz (-d) and sends the data to R (-c)
# 'rb' = read binary mode
con <- pipe(sprintf("pigz -dc %s", input_file), "rb")

# Read the object from the connection
print("Reading object with parallel decompression via pigz...")
resultados_lidos <- readRDS(file = con)

# Close the connection
close(con)
print("Read successful!")



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# SCRIPT FOR PARALLEL JSON PARSING WITH future_map
# -------------------------------------------------------------------------
# OBJECTIVE:
# Convert a character vector containing JSON strings ('positions_vec')
# into a list of R objects, using parallel processing to speed up
# the task.
# -------------------------------------------------------------------------


# =========================================================================
# STEP 1: SETUP
# =========================================================================

# -- 1.1: Load the required packages
# furrr for parallel processing, future for configuration,
# jsonlite for parsing, and progressr for the progress bar.
library(furrr)
library(future)
library(jsonlite)
library(progressr)


# =========================================================================
# STEP 2: PARALLEL ENVIRONMENT CONFIGURATION
# =========================================================================

# -- 2.1: Define the parallel execution plan
# 'multisession' starts 35 clean R sessions in the background,
# providing robust and memory-safe isolation between processes.
plan(multisession, workers = 35)
cat(sprintf("Parallel environment configured with %d workers.\n", nbrOfWorkers()))

# -- 2.2: Configure worker options
# Ensure that each of the 35 workers loads the 'jsonlite' package before starting.
# 'seed = TRUE' ensures reproducibility in processes that use random numbers.
options_furrr <- furrr_options(packages = "jsonlite", seed = TRUE)


# =========================================================================
# STEP 3: PARALLEL PARSING
# =========================================================================

# 'positions_vec' should be the character vector already available in memory

cat("Initializing the parsing of", length(resultados_lidos), "records in parellel...\n")

# -- Measure execution time
tempo_execucao_parallel <- system.time({
  
  # Enable progress-bar monitoring
  with_progress({
    
    # Apply fromJSON to each element of 'positions_vec' in parallel
    final_results_list <- future_map(
      .x = resultados_lidos$positions,         # Vector over which to iterate
      .f = fromJSON,              # Function applied to each element
      .options = options_furrr,   # Worker options
      .progress = TRUE            # Enable the progress bar
    )
    
  }) # End of with_progress
  
}) # End of system.time


# =========================================================================
# STEP 4: CLEANUP AND REPORTING
# =========================================================================

# -- 4.1: Release the workers
# Return to sequential mode to close the background sessions.
plan(sequential)
cat("Parallel processing concluded. Workers released.\n")

# -- 4.2: Report execution time and results
cat(sprintf("Total execution time: %.2f seconds\n", tempo_execucao_parallel['elapsed']))
cat(sprintf("Process completed. 'final_results_list' contains %d elements.\n", length(final_results_list)))


# --- Result verification ---
glimpse(final_results_list[[1]])


# -------------------------------------------------------------------------

# Define the output filename
output_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes_positions_list.rds.gz"

# Create a pipe connection to pigz
# The command runs pigz with 16 threads and redirects the output to output_file
# 'wb' = write binary mode
con <- pipe(sprintf("pigz -p 16 > %s", output_file), "wb")

# Use saveRDS with this connection and DISABLE R internal compression
# R writes the raw data to the pipe as quickly as possible, while pigz handles the compression.
print("Saving object with parallel compression via pigz...")
saveRDS(object = final_results_list, file = con, compress = FALSE)

# Closing the connection is ESSENTIAL to ensure that all data are written correctly
close(con)
print("Save successful!")


### GENES

cat("Initializing parsing of", length(resultados_lidos), "records in parallel...\n")

# -- Measure execution time
tempo_execucao_parallel <- system.time({
  
  # Enable progress-bar monitoring
  with_progress({
    
    # Apply fromJSON to each element of 'positions_vec' in parallel
    final_results_list <- future_map(
      .x = resultados_lidos$genes,         # Vector over which to iterate
      .f = fromJSON,              # Function applied to each element
      .options = options_furrr,   # Worker options
      .progress = TRUE            # Enable the progress bar
    )
    
  }) # End of with_progress
  
}) # End of system.time


# =========================================================================
# STEP 4: CLEANUP AND REPORTING
# =========================================================================

# -- 4.1: Release the workers
# Return to sequential mode to close the background sessions.
plan(sequential)
cat("Parallel processing completed. Workers released.\n")

# -- 4.2: Report execution time and results
cat(sprintf("Total execution time: %.2f seconds\n", tempo_execucao_parallel['elapsed']))
cat(sprintf("Process completed. 'final_results_list' contains %d elements.\n", length(final_results_list)))


# --- Result verification ---
glimpse(final_results_list[[1]])


# -------------------------------------------------------------------------

# Define the output filename
output_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_genes_list.rds.gz"

# Create a pipe connection to pigz
# The command runs pigz with 16 threads and redirects the output to output_file
# 'wb' = write binary mode
con <- pipe(sprintf("pigz -p 16 > %s", output_file), "wb")

# Use saveRDS with this connection and DISABLE R internal compression
# R writes the raw data to the pipe as quickly as possible, while pigz handles the compression.
print("Saving object with parallel compression via pigz...")
saveRDS(object = final_results_list, file = con, compress = FALSE)

# Closing the connection is ESSENTIAL to ensure that all data are written correctly
close(con)
print("Save successful!")


header <- resultados_lidos$header %>% fromJSON() 
saveRDS(object = header, file = "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_header.rds")


# -------------------------------------------------------------------------

header <- readRDS("multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_header.rds")
input_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_com_genes_positions_list.rds.gz"

# The command decompresses the file with pigz (-d) and sends the data to R (-c)
# 'rb' = read binary mode
con <- pipe(sprintf("pigz -dc %s", input_file), "rb")

# Read the object from the connection
print("Reading object with parallel decompression via pigz...")
positions_list <- readRDS(file = con)

# Close the connection
close(con)
print("Read successful!")

input_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_genes_list.rds.gz"

# The command decompresses the file with pigz (-d) and sends the data to R (-c)
# 'rb' = read binary mode
con <- pipe(sprintf("pigz -dc %s", input_file), "rb")

# Read the object from the connection
print("Reading object with parallel decompression via pigz...")
genes_list <- readRDS(file = con)

# Close the connection
close(con)
print("Read successful!")

# -------------------------------------------------------------------------

library(purrr)
library(dplyr)
library(tidyr)

# 1. Use purrr::map to extract the names (names()) of each list element
all_names <- map(positions_list, names)
all_names_genes <- map(genes_list, names)

# 2. 'all_names' is now a list of character vectors.
#    Combine them into a single vector with unlist() and retrieve the unique values.
unique_names <- unique(unlist(all_names))
unique_names_genes <- unique(unlist(all_names_genes))


# 3. Display the result
print("Unique fields found in list:")
print(unique_names)

# 3. Display the result
print("Unique fields found in list:")
print(unique_names_genes)

tbl_inicial <- tibble(
  
  chromosome = map_chr(positions_list, "chromosome", .default = NA_character_),
  position = map_int(positions_list, "position", .default = NA_integer_),
  refAllele = map_chr(positions_list, "refAllele", .default = NA_character_),
  altAlleles = map(positions_list, "altAlleles", .default = list(NULL)),
  vcfInfo = map(positions_list, "vcfInfo"),
  quality = map_dbl(positions_list, "quality", .default = NA_real_),
  filters = map_chr(positions_list, "filters", .default = NA_character_),
  cytogeneticBand = map_chr(positions_list, "cytogeneticBand", .default = NA_character_),
  samples = map(positions_list, "samples"),
  variants = map(positions_list, "variants")
)

tbl_genes<- tibble(
  name = map_chr(genes_list, "name", .default = NA_character_),
  ensemblGeneId = map_chr(genes_list, "ensemblGeneId", .default = NA_character_),
  hgncId = map_int(genes_list, "hgncId", .default = NA_integer_),
  ncbiGeneId = map_chr(genes_list, "ncbiGeneId", .default = NA_character_),
  gnomAD = map(genes_list, "gnomAD", .default = list(NULL)),
  omim = map(genes_list, "omim", .default = list(NULL)),
  clingenGeneValidity = map(genes_list, "clingenGeneValidity", .default = list(NULL)),
  clingenDosageSensitivityMap = map(genes_list, "clingenDosageSensitivityMap", .default = list(NULL)),
  cosmic = map(genes_list, "cosmic", .default = list(NULL))
)


# -------------------------------------------------------------------------

# Define the output filename
output_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_inicial.rds.gz"

# Create a pipe connection to pigz
# The command runs pigz with 16 threads and redirects the output to output_file
# 'wb' = write binary mode
con <- pipe(sprintf("pigz -p 16 > %s", output_file), "wb")

# Use saveRDS with this connection and DISABLE R internal compression
# R writes the raw data to the pipe as quickly as possible, while pigz handles the compression.
print("Saving object with parallel compression via pigz...")
saveRDS(object = tbl_inicial, file = con, compress = FALSE)

# Closing the connection is ESSENTIAL to ensure that all data are written correctly
close(con)
print("Save successful!")


# Define the output filename
output_file <- "multi_gvcf_FINAL_FMISSING_AF_contextual_aPRoVAR_22052026_tbl_genes.rds.gz"

# Create a pipe connection to pigz
# The command runs pigz with 16 threads and redirects the output to output_file
# 'wb' = write binary mode
con <- pipe(sprintf("pigz -p 16 > %s", output_file), "wb")

# Use saveRDS with this connection and DISABLE R internal compression
# R writes the raw data to the pipe as quickly as possible, while pigz handles the compression.
print("Saving object with parallel compression via pigz...")
saveRDS(object = tbl_genes, file = con, compress = FALSE)

# Closing the connection is ESSENTIAL to ensure that all data are written correctly
close(con)
print("Save successful!")