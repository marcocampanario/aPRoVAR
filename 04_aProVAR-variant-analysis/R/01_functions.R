##### 01_functions.R #####
# Marco Antonio Campanário @ Fiocruz Paraná, 2026 Jul.
# Provide reusable input/output, nested-annotation, genotype, classification,
# plotting, and statistical helper functions for the aPRoVAR workflow.

# General validation and input/output -----------------------------------------

require_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this module.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

create_output_directories <- function() {
  directories <- c(
    RAW_DATA_DIR,
    PROCESSED_DATA_DIR,
    TABLES_DIR,
    FIGURES_DIR,
    REPORTS_DIR
  )

  invisible(
    vapply(
      directories,
      dir.create,
      logical(1),
      recursive = TRUE,
      showWarnings = FALSE
    )
  )
}

assert_file_exists <- function(path, description = "Input file") {
  if (!file.exists(path)) {
    stop(description, " was not found: ", path, call. = FALSE)
  }

  invisible(path)
}

assert_columns <- function(data, required_columns, object_name = "data") {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

read_rds_parallel <- function(path) {
  assert_file_exists(path, "RDS file")
  pigz <- Sys.which("pigz")

  if (!nzchar(pigz) || !grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }

  command <- sprintf("%s -dc %s", shQuote(pigz), shQuote(path))
  connection <- pipe(command, open = "rb")
  on.exit(close(connection), add = TRUE)
  readRDS(connection)
}

write_rds_parallel <- function(object, path, threads = PIGZ_THREADS) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  pigz <- Sys.which("pigz")

  if (!nzchar(pigz) || !grepl("\\.gz$", path, ignore.case = TRUE)) {
    saveRDS(object, path, compress = grepl("\\.gz$", path, ignore.case = TRUE))
    return(invisible(path))
  }

  command <- sprintf(
    "%s -p %d > %s",
    shQuote(pigz),
    as.integer(threads),
    shQuote(path)
  )

  connection <- pipe(command, open = "wb")
  saveRDS(object, connection, compress = FALSE)
  close(connection)
  invisible(path)
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(path)
}

write_report <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(lines), con = path)
  invisible(path)
}

save_plot_formats <- function(
    plot,
    stem,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    dpi = FIGURE_DPI,
    save_tiff = TRUE,
    save_png = TRUE) {
  require_packages("ggplot2")
  dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(save_tiff)) {
    ggplot2::ggsave(
      filename = file.path(FIGURES_DIR, paste0(stem, ".tiff")),
      plot = plot,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
  }

  if (isTRUE(save_png)) {
    ggplot2::ggsave(
      filename = file.path(FIGURES_DIR, paste0(stem, ".png")),
      plot = plot,
      device = "png",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }

  invisible(plot)
}

# Generic nested-object helpers -----------------------------------------------

is_empty_value <- function(value) {
  is.null(value) || length(value) == 0L || all(is.na(value))
}

first_non_missing_numeric <- function(value) {
  if (is_empty_value(value)) {
    return(NA_real_)
  }

  numeric_value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  numeric_value <- numeric_value[is.finite(numeric_value)]

  if (length(numeric_value) == 0L) NA_real_ else numeric_value[[1L]]
}

# Parse ANNOVAR numeric fields consistently whether fread() inferred the
# source column as character, factor, integer, or double.
parse_annovar_numeric <- function(value) {
  value <- trimws(as.character(value))
  value[is.na(value) | value %in% c("", ".")] <- NA_character_
  suppressWarnings(as.numeric(value))
}

collapse_unique <- function(value, separator = ";") {
  if (is_empty_value(value)) {
    return(NA_character_)
  }

  value <- trimws(as.character(unlist(value, use.names = FALSE)))
  value <- unique(value[!is.na(value) & nzchar(value)])

  if (length(value) == 0L) NA_character_ else paste(value, collapse = separator)
}

unwrap_data_frame <- function(value) {
  if (is.data.frame(value)) {
    return(value)
  }

  if (is.list(value) && length(value) == 1L && is.data.frame(value[[1L]])) {
    return(value[[1L]])
  }

  NULL
}

normalize_variant_annotations <- function(value) {
  if (is.null(value)) {
    return(list(NULL))
  }

  if (is.data.frame(value)) {
    return(list(value))
  }

  if (is.list(value) && !is.data.frame(value)) {
    return(value)
  }

  list(value)
}

match_variant_annotations_to_alts <- function(variant, number_of_alts) {
  if (number_of_alts == 0L) {
    return(list())
  }

  if (is.null(variant)) {
    return(rep(list(NULL), number_of_alts))
  }

  if (is.data.frame(variant)) {
    return(rep(list(variant), number_of_alts))
  }

  # An unnamed list with one annotation object per ALT allele is the standard
  # multiallelic representation in the parsed Nirvana object.
  if (
    is.list(variant) &&
      is.null(names(variant)) &&
      length(variant) == number_of_alts
  ) {
    return(variant)
  }

  # A named list generally represents fields of one annotation object and must
  # remain intact, even if its length happens to equal the number of ALT alleles.
  rep(list(variant), number_of_alts)
}

extract_named_field <- function(value, field) {
  if (is.null(value)) {
    return(NULL)
  }

  if (is.data.frame(value) && field %in% names(value)) {
    return(value[[field]])
  }

  if (is.list(value) && field %in% names(value)) {
    return(value[[field]])
  }

  if (is.list(value)) {
    nested_values <- lapply(value, extract_named_field, field = field)
    nested_values <- nested_values[
      !vapply(nested_values, is_empty_value, logical(1))
    ]

    if (length(nested_values) > 0L) {
      return(unlist(nested_values, recursive = TRUE, use.names = FALSE))
    }
  }

  NULL
}

has_named_annotation <- function(variant, annotation_name) {
  if (is.null(variant)) {
    return(FALSE)
  }

  if (!is.null(names(variant)) && annotation_name %in% names(variant)) {
    return(TRUE)
  }

  if (is.list(variant)) {
    return(any(vapply(
      variant,
      function(item) {
        !is.null(names(item)) && annotation_name %in% names(item)
      },
      logical(1)
    )))
  }

  FALSE
}

extract_annotation_block <- function(variant, annotation_name) {
  if (is.null(variant)) {
    return(NULL)
  }

  if (!is.null(names(variant)) && annotation_name %in% names(variant)) {
    return(variant[[annotation_name]])
  }

  if (is.list(variant)) {
    for (item in variant) {
      if (!is.null(names(item)) && annotation_name %in% names(item)) {
        return(item[[annotation_name]])
      }
    }
  }

  NULL
}

extract_population_all_af <- function(variant, annotation_name) {
  block <- extract_annotation_block(variant, annotation_name)
  data <- unwrap_data_frame(block)

  if (is.null(data) || !"allAf" %in% names(data)) {
    return(NA_real_)
  }

  values <- suppressWarnings(as.numeric(data$allAf))
  values <- values[is.finite(values)]

  if (length(values) == 0L) NA_real_ else mean(values)
}

# Keep both gnomAD resources explicit. For a single reporting frequency,
# prefer gnomAD Genomes and use gnomAD Exomes only when Genomes is unavailable.
select_preferred_gnomad_af <- function(genome_af, exome_af) {
  genome_af <- suppressWarnings(as.numeric(genome_af))
  exome_af <- suppressWarnings(as.numeric(exome_af))

  if (length(genome_af) > 0L && is.finite(genome_af[[1L]])) {
    return(genome_af[[1L]])
  }

  if (length(exome_af) > 0L && is.finite(exome_af[[1L]])) {
    return(exome_af[[1L]])
  }

  NA_real_
}

is_positive_af <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  is.finite(value) & value > 0
}

is_known_from_population_af <- function(
    gnomad_genome_af,
    gnomad_exome_af,
    one_kg_af,
    abraom_af) {
  is_positive_af(gnomad_genome_af) |
    is_positive_af(gnomad_exome_af) |
    is_positive_af(one_kg_af) |
    is_positive_af(abraom_af)
}

extract_dbsnp_id <- function(variant) {
  collapse_unique(extract_named_field(variant, "dbsnp"), separator = ",")
}

extract_variant_type <- function(variant) {
  collapse_unique(extract_named_field(variant, "variantType"), separator = ";")
}

standardize_variant_type <- function(variant_type) {
  if (is.na(variant_type) || !nzchar(variant_type)) {
    return("other")
  }

  if (grepl("SNV", variant_type, ignore.case = TRUE)) {
    return("SNV")
  }

  if (grepl("insertion", variant_type, ignore.case = TRUE)) {
    return("insertion")
  }

  if (grepl("deletion", variant_type, ignore.case = TRUE)) {
    return("deletion")
  }

  "other"
}

extract_transcript_table <- function(variant) {
  transcripts <- extract_annotation_block(variant, "transcripts")

  if (is.null(transcripts)) {
    return(NULL)
  }

  if (is.data.frame(transcripts)) {
    return(transcripts)
  }

  if (is.list(transcripts)) {
    transcript_tables <- transcripts[vapply(transcripts, is.data.frame, logical(1))]

    if (length(transcript_tables) > 0L) {
      return(dplyr::bind_rows(transcript_tables))
    }
  }

  NULL
}

extract_genes <- function(variant) {
  transcripts <- extract_transcript_table(variant)

  if (is.null(transcripts) || !"hgnc" %in% names(transcripts)) {
    return(NA_character_)
  }

  collapse_unique(transcripts$hgnc, separator = ";")
}

# VCF and genotype helpers -----------------------------------------------------

extract_vcf_info <- function(vcf_info, field) {
  value <- if (is.null(vcf_info)) NULL else vcf_info[[field]]

  if (is_empty_value(value)) {
    return(NA_character_)
  }

  collapse_unique(value, separator = ",")
}

parse_af_vector <- function(value) {
  if (is_empty_value(value)) {
    return(NA_real_)
  }

  if (length(value) == 1L && is.character(value)) {
    value <- strsplit(value, ",", fixed = TRUE)[[1L]]
  }

  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) == 0L) NA_real_ else parsed
}

choose_analysis_af <- function(af, af_excl) {
  if (!is_empty_value(af_excl)) parse_af_vector(af_excl) else parse_af_vector(af)
}

make_variant_id <- function(chromosome, position, ref, alt) {
  alt_text <- vapply(
    alt,
    function(value) paste(as.character(unlist(value)), collapse = ","),
    character(1)
  )

  paste(chromosome, position, ref, alt_text, sep = "-")
}

# ANNOVAR row-alignment helpers -----------------------------------------------

normalize_chromosome_label <- function(chromosome) {
  chromosome <- trimws(as.character(chromosome))
  sub("^chr", "", chromosome, ignore.case = TRUE)
}

make_variant_alignment_key <- function(chromosome, position, ref, alt) {
  paste(
    normalize_chromosome_label(chromosome),
    suppressWarnings(as.integer(position)),
    toupper(trimws(as.character(ref))),
    toupper(trimws(as.character(alt))),
    sep = "|"
  )
}

# This standard conversion handles ordinary VCF insertions and deletions. It is
# used only to validate row alignment. ANNOVAR can additionally simplify or
# reposition complex and repetitive alleles, so this key must not be used as a
# complete join key for every variant class.
make_standard_annovar_key_from_vcf <- function(
    chromosome,
    position,
    ref,
    alt) {
  chromosome <- normalize_chromosome_label(chromosome)
  position <- suppressWarnings(as.integer(position))
  ref <- toupper(trimws(as.character(ref)))
  alt <- toupper(trimws(as.character(alt)))

  ref_length <- nchar(ref)
  alt_length <- nchar(alt)
  valid <- !is.na(ref) & !is.na(alt) & !is.na(position)

  ref_is_alt_prefix <- valid &
    alt_length > ref_length &
    substr(alt, 1L, ref_length) == ref

  alt_is_ref_prefix <- valid &
    ref_length > alt_length &
    substr(ref, 1L, alt_length) == alt

  annovar_start <- position
  annovar_ref <- ref
  annovar_alt <- alt

  annovar_start[ref_is_alt_prefix] <-
    position[ref_is_alt_prefix] + ref_length[ref_is_alt_prefix] - 1L
  annovar_ref[ref_is_alt_prefix] <- "-"
  annovar_alt[ref_is_alt_prefix] <- substring(
    alt[ref_is_alt_prefix],
    ref_length[ref_is_alt_prefix] + 1L
  )

  annovar_start[alt_is_ref_prefix] <-
    position[alt_is_ref_prefix] + alt_length[alt_is_ref_prefix]
  annovar_ref[alt_is_ref_prefix] <- substring(
    ref[alt_is_ref_prefix],
    alt_length[alt_is_ref_prefix] + 1L
  )
  annovar_alt[alt_is_ref_prefix] <- "-"

  paste(
    chromosome,
    annovar_start,
    annovar_ref,
    annovar_alt,
    sep = "|"
  )
}

validate_annovar_row_alignment <- function(
    allele_table,
    annovar_table,
    minimum_exact_match = 0.90,
    minimum_standard_match = 0.98) {
  assert_columns(
    allele_table,
    c("chromosome", "position", "refAllele", "altAllele"),
    "ALT-allele table"
  )
  assert_columns(
    annovar_table,
    c("Chr", "Start", "Ref", "Alt"),
    "ANNOVAR table"
  )

  allele_records <- nrow(allele_table)
  annovar_records <- nrow(annovar_table)

  if (allele_records != annovar_records) {
    stop(
      "The ALT-allele and ANNOVAR tables have different row counts (",
      allele_records,
      " versus ",
      annovar_records,
      "). Validated row-index alignment cannot be used.",
      call. = FALSE
    )
  }

  proportion_true <- function(value) {
    if (length(value) == 0L) return(NA_real_)
    sum(value %in% TRUE) / length(value)
  }

  same_chromosome <-
    normalize_chromosome_label(allele_table$chromosome) ==
    normalize_chromosome_label(annovar_table$Chr)

  exact_key_match <-
    make_variant_alignment_key(
      allele_table$chromosome,
      allele_table$position,
      allele_table$refAllele,
      allele_table$altAllele
    ) ==
    make_variant_alignment_key(
      annovar_table$Chr,
      annovar_table$Start,
      annovar_table$Ref,
      annovar_table$Alt
    )

  standard_key_match <-
    make_standard_annovar_key_from_vcf(
      allele_table$chromosome,
      allele_table$position,
      allele_table$refAllele,
      allele_table$altAllele
    ) ==
    make_variant_alignment_key(
      annovar_table$Chr,
      annovar_table$Start,
      annovar_table$Ref,
      annovar_table$Alt
    )

  alignment <- list(
    allele_records = allele_records,
    annovar_records = annovar_records,
    same_chromosome_count = sum(same_chromosome %in% TRUE),
    same_chromosome_rate = proportion_true(same_chromosome),
    exact_key_match_count = sum(exact_key_match %in% TRUE),
    exact_key_match_rate = proportion_true(exact_key_match),
    standard_key_match_count = sum(standard_key_match %in% TRUE),
    standard_key_match_rate = proportion_true(standard_key_match),
    minimum_exact_match = minimum_exact_match,
    minimum_standard_match = minimum_standard_match
  )

  if (!isTRUE(alignment$same_chromosome_rate == 1)) {
    stop(
      "Chromosome labels are not aligned on every corresponding row of the ",
      "ALT-allele and ANNOVAR tables.",
      call. = FALSE
    )
  }

  if (
    !is.finite(alignment$exact_key_match_rate) ||
      alignment$exact_key_match_rate < minimum_exact_match
  ) {
    stop(
      "The same-row exact variant-key match rate is ",
      round(100 * alignment$exact_key_match_rate, 3),
      "%; at least ",
      round(100 * minimum_exact_match, 3),
      "% is required.",
      call. = FALSE
    )
  }

  if (
    !is.finite(alignment$standard_key_match_rate) ||
      alignment$standard_key_match_rate < minimum_standard_match
  ) {
    stop(
      "The same-row standard-normalized variant-key match rate is ",
      round(100 * alignment$standard_key_match_rate, 3),
      "%; at least ",
      round(100 * minimum_standard_match, 3),
      "% is required.",
      call. = FALSE
    )
  }

  alignment
}

extract_sample_names <- function(header) {
  sample_names <- header$samples

  if (is.null(sample_names) || length(sample_names) == 0L) {
    stop("No sample names were found in header$samples.", call. = FALSE)
  }

  as.character(sample_names)
}

extract_genotypes <- function(sample_object) {
  if (is.null(sample_object)) {
    return(character(0))
  }

  if (is.data.frame(sample_object)) {
    genotype_column <- intersect(
      c("genotype", "GT", "gt", "genotype_call"),
      names(sample_object)
    )

    if (length(genotype_column) > 0L) {
      return(as.character(sample_object[[genotype_column[[1L]]]]))
    }
  }

  if (is.list(sample_object) && !is.null(sample_object$genotype)) {
    return(as.character(sample_object$genotype))
  }

  character(0)
}

genotype_has_alt <- function(genotype) {
  if (length(genotype) == 0L || is.na(genotype)) {
    return(FALSE)
  }

  genotype <- trimws(as.character(genotype))

  if (genotype %in% c("", ".", "./.", ".|.", ".:.")) {
    return(FALSE)
  }

  alleles <- strsplit(genotype, "[/|]")[[1L]]
  any(!alleles %in% c("0", "."))
}

make_presence_matrix <- function(sample_list, sample_names, row_ids) {
  expected_length <- length(sample_names)

  presence_list <- lapply(sample_list, function(sample_object) {
    genotypes <- extract_genotypes(sample_object)

    if (length(genotypes) != expected_length) {
      stop(
        "A sample-level genotype object contained ",
        length(genotypes),
        " calls; expected ",
        expected_length,
        ".",
        call. = FALSE
      )
    }

    vapply(genotypes, genotype_has_alt, logical(1))
  })

  presence_matrix <- do.call(rbind, presence_list)
  rownames(presence_matrix) <- row_ids
  colnames(presence_matrix) <- sample_names
  presence_matrix
}

count_site_alleles <- function(sample_object) {
  genotypes <- extract_genotypes(sample_object)

  if (length(genotypes) == 0L) {
    return(list(
      allele_number = NA_integer_,
      alt_counts = integer(0),
      total_alt_count = NA_integer_,
      is_doubleton_hom = FALSE
    ))
  }

  split_alleles <- lapply(genotypes, function(genotype) {
    if (is.na(genotype) || genotype %in% c("", ".", "./.", ".|.")) {
      return(character(0))
    }

    alleles <- strsplit(as.character(genotype), "[/|]")[[1L]]
    alleles[alleles != "."]
  })

  allele_number <- sum(lengths(split_alleles))
  all_alt_alleles <- unlist(split_alleles, use.names = FALSE)
  all_alt_alleles <- all_alt_alleles[all_alt_alleles != "0"]

  if (length(all_alt_alleles) == 0L) {
    alt_counts <- integer(0)
  } else {
    alt_counts <- table(all_alt_alleles)
    alt_counts <- as.integer(alt_counts[order(as.integer(names(alt_counts)))])
  }

  hom_alt <- vapply(split_alleles, function(alleles) {
    length(alleles) == 2L &&
      length(unique(alleles)) == 1L &&
      alleles[[1L]] != "0"
  }, logical(1))

  list(
    allele_number = as.integer(allele_number),
    alt_counts = alt_counts,
    total_alt_count = as.integer(length(all_alt_alleles)),
    is_doubleton_hom = length(all_alt_alleles) == 2L && sum(hom_alt) == 1L
  )
}

classify_frequency_from_counts <- function(
    alt_counts,
    allele_number,
    is_doubleton_hom = FALSE,
    polymorphic_threshold = POLYMORPHIC_AF_THRESHOLD) {
  alt_counts <- suppressWarnings(as.integer(alt_counts))
  alt_counts <- alt_counts[!is.na(alt_counts) & alt_counts > 0L]

  if (length(alt_counts) == 0L || is.na(allele_number) || allele_number <= 0L) {
    return(NA_character_)
  }

  allele_frequencies <- alt_counts / allele_number

  if (any(allele_frequencies >= polymorphic_threshold)) {
    return("polymorphic")
  }

  if (any(alt_counts > 2L)) {
    return("rare")
  }

  if (any(alt_counts == 2L)) {
    if (isTRUE(is_doubleton_hom) && sum(alt_counts) == 2L) {
      return("doubleton_hom")
    }

    return("doubleton")
  }

  if (any(alt_counts == 1L)) {
    return("singleton")
  }

  NA_character_
}

classify_origin <- function(is_known, dbsnp_id) {
  if (isTRUE(is_known)) {
    return("known")
  }

  if (!is.na(dbsnp_id) && nzchar(dbsnp_id)) {
    return("absent")
  }

  "novel"
}

format_origin <- function(origin) {
  factor(
    origin,
    levels = c("known", "absent", "novel"),
    labels = c("Known", "Absent", "Novel")
  )
}

# ClinVar helpers --------------------------------------------------------------

extract_clinvar_classification <- function(variant) {
  clinvar_preview <- extract_annotation_block(variant, "clinvar-preview")

  if (is_empty_value(clinvar_preview)) {
    return(NA_character_)
  }

  classifications <- tryCatch(
    clinvar_preview[[1L]]$classifications$germlineClassification$classification,
    error = function(error) NULL
  )

  collapse_unique(classifications, separator = ";")
}

normalize_clinvar_terms <- function(classification) {
  if (is.na(classification) || !nzchar(classification)) {
    return(NA_character_)
  }

  trimws(strsplit(classification, ";", fixed = TRUE)[[1L]])
}

collapse_clinvar_category <- function(classification) {
  terms <- normalize_clinvar_terms(classification)

  if (length(terms) == 1L && is.na(terms)) {
    return("Not provided")
  }

  pathogenic_terms <- c(
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Likely pathogenic/Pathogenic"
  )

  benign_terms <- c(
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
    "Likely benign/Benign"
  )

  low_penetrance_terms <- c(
    "Pathogenic/Likely pathogenic/Pathogenic, low penetrance",
    "Pathogenic/Pathogenic, low penetrance",
    "Likely risk allele",
    "Uncertain risk allele",
    "Uncertain significance/Uncertain risk allele"
  )

  if (all(terms %in% pathogenic_terms)) {
    return("Likely pathogenic/Pathogenic")
  }

  if (all(terms %in% benign_terms)) {
    return("Likely benign/Benign")
  }

  if (all(terms == "Uncertain significance")) {
    return("Uncertain significance")
  }

  if (all(terms == "Conflicting classifications of pathogenicity")) {
    return("Conflicting classifications of pathogenicity")
  }

  if (all(tolower(terms) == "drug response")) {
    return("Drug response")
  }

  if (all(tolower(terms) == "affects")) {
    return("Affects a non-disease phenotype")
  }

  if (all(tolower(terms) == "protective")) {
    return("Protective")
  }

  if (all(terms %in% low_penetrance_terms)) {
    return("Low penetrance for Mendelian diseases")
  }

  if (all(tolower(terms) %in% c("not provided", "na"))) {
    return("Not provided")
  }

  if (all(tolower(terms) == "association")) {
    return("GWAS hits")
  }

  if (all(tolower(terms) == "risk factor")) {
    return("Risk factor")
  }

  "Other"
}

extract_clinvar_phenotypes <- function(variant) {
  clinvar <- extract_annotation_block(variant, "clinvar")

  if (is_empty_value(clinvar)) {
    return(NA_character_)
  }

  phenotypes <- tryCatch(
    clinvar[[1L]]$phenotypes,
    error = function(error) NULL
  )

  collapse_unique(phenotypes, separator = "; ")
}

# Transcript impact and REVEL helpers -----------------------------------------

extract_selected_transcript_impact <- function(
    variant,
    source = "RefSeq",
    require_canonical = TRUE,
    require_mane = TRUE) {
  transcripts <- extract_transcript_table(variant)

  if (is.null(transcripts) || nrow(transcripts) == 0L) {
    return(NA_character_)
  }

  keep <- rep(TRUE, nrow(transcripts))

  if ("source" %in% names(transcripts)) {
    keep <- keep & transcripts$source == source
  } else {
    keep <- rep(FALSE, nrow(transcripts))
  }

  if (isTRUE(require_canonical)) {
    if ("isCanonical" %in% names(transcripts)) {
      keep <- keep & !is.na(transcripts$isCanonical) & transcripts$isCanonical
    } else {
      keep <- rep(FALSE, nrow(transcripts))
    }
  }

  if (isTRUE(require_mane)) {
    if ("isManeSelect" %in% names(transcripts)) {
      keep <- keep & !is.na(transcripts$isManeSelect) & transcripts$isManeSelect
    } else {
      keep <- rep(FALSE, nrow(transcripts))
    }
  }

  selected <- transcripts[keep, , drop = FALSE]

  if (nrow(selected) == 0L || !"impact" %in% names(selected)) {
    return(NA_character_)
  }

  collapse_unique(tolower(selected$impact), separator = ";")
}

extract_revel_max <- function(variant) {
  revel <- extract_annotation_block(variant, "revel")
  score <- extract_named_field(revel, "score")

  if (is_empty_value(score)) {
    return(NA_real_)
  }

  score <- suppressWarnings(as.numeric(score))
  score <- score[is.finite(score)]

  if (length(score) == 0L) NA_real_ else max(score)
}

# Spreadsheet helpers ---------------------------------------------------------

drop_list_columns <- function(data) {
  data[, !vapply(data, is.list, logical(1)), drop = FALSE]
}

write_workbook <- function(sheets, path) {
  require_packages("openxlsx")
  workbook <- openxlsx::createWorkbook()

  for (sheet_name in names(sheets)) {
    safe_name <- substr(sheet_name, 1L, 31L)
    openxlsx::addWorksheet(workbook, safe_name)
    openxlsx::writeData(
      workbook,
      sheet = safe_name,
      x = sheets[[sheet_name]],
      keepNA = FALSE
    )
  }

  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  invisible(path)
}

# Population-frequency helpers ------------------------------------------------

extract_gnomad_numeric <- function(variant, field) {
  gnomad <- extract_annotation_block(variant, "gnomad")
  first_non_missing_numeric(extract_named_field(gnomad, field))
}

safe_fisher_test <- function(local_ac, local_an, external_ac, external_an) {
  counts <- c(local_ac, local_an, external_ac, external_an)

  if (
    any(!is.finite(counts)) ||
      any(counts < 0) ||
      local_ac > local_an ||
      external_ac > external_an
  ) {
    return(c(p_value = NA_real_, odds_ratio = NA_real_))
  }

  result <- tryCatch(
    stats::fisher.test(
      matrix(
        c(
          local_ac,
          local_an - local_ac,
          external_ac,
          external_an - external_ac
        ),
        nrow = 2L,
        byrow = TRUE
      )
    ),
    error = function(error) NULL
  )

  if (is.null(result)) {
    return(c(p_value = NA_real_, odds_ratio = NA_real_))
  }

  c(
    p_value = unname(result$p.value),
    odds_ratio = unname(as.numeric(result$estimate))
  )
}

x_log_x <- function(count, probability) {
  output <- rep(NA_real_, length(count))
  zero_count <- !is.na(count) & count == 0
  valid <- !is.na(count) & count > 0 & is.finite(probability) & probability > 0
  output[zero_count] <- 0
  output[valid] <- count[valid] * log(probability[valid])
  output
}

safe_lrt <- function(local_ac, local_an, external_ac, external_an) {
  counts <- c(local_ac, local_an, external_ac, external_an)

  if (
    any(!is.finite(counts)) ||
      any(counts < 0) ||
      local_ac > local_an ||
      external_ac > external_an
  ) {
    return(c(statistic = NA_real_, p_value = NA_real_))
  }

  local_ref <- local_an - local_ac
  external_ref <- external_an - external_ac
  pooled_ac <- local_ac + external_ac
  pooled_ref <- local_ref + external_ref
  pooled_an <- pooled_ac + pooled_ref

  local_af <- local_ac / local_an
  external_af <- external_ac / external_an
  pooled_af <- pooled_ac / pooled_an

  alternative_log_likelihood <- sum(c(
    x_log_x(local_ac, local_af),
    x_log_x(local_ref, 1 - local_af),
    x_log_x(external_ac, external_af),
    x_log_x(external_ref, 1 - external_af)
  ))

  null_log_likelihood <- sum(c(
    x_log_x(pooled_ac, pooled_af),
    x_log_x(pooled_ref, 1 - pooled_af)
  ))

  statistic <- 2 * (alternative_log_likelihood - null_log_likelihood)
  statistic <- max(statistic, 0)

  c(
    statistic = statistic,
    p_value = stats::pchisq(statistic, df = 1, lower.tail = FALSE)
  )
}

make_rare_donut <- function(data, origin_name) {
  plot_data <- data |>
    dplyr::filter(Origin == origin_name) |>
    dplyr::mutate(
      rare_component = factor(
        rare_component,
        levels = rare_components
      )
    )
  
  # Generate three shades from the origin-specific color.
  component_colors <- grDevices::colorRampPalette(
    c("#F2F2F2", ORIGIN_COLORS[[origin_name]])
  )(4)[2:4]
  
  names(component_colors) <- rare_components
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = 2,
      y = proportion,
      fill = rare_component
    )
  ) +
    ggplot2::geom_col(
      color = "white",
      linewidth = 0.6,
      width = 1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 3,
        label = label
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3,
      lineheight = 0.9
    ) +
    ggplot2::coord_polar(
      theta = "y",
      clip = "off"
    ) +
    ggplot2::xlim(0.5, 3.2) +
    ggplot2::scale_fill_manual(
      values = component_colors,
      breaks = rare_components,
      drop = FALSE
    ) +
    ggplot2::labs(
      fill = NULL,
      title = paste0(origin_name, " rare variants")
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = ggplot2::margin(
        t = 15,
        r = 35,
        b = 15,
        l = 35
      )
    )
}

make_gene_wordcloud <- function(data, origin_name) {
  plot_data <- data |>
    dplyr::filter(Origin == origin_name) |>
    dplyr::arrange(
      dplyr::desc(n_variants),
      genes
    )
  
  if (nrow(plot_data) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      label = genes,
      size = n_variants,
      color = cor_categoria,
      angle = angle
    )
  ) +
    ggwordcloud::geom_text_wordcloud_area(
      grid_size = 2,
      eccentricity = 1,
      rm_outside = FALSE
    ) +
    ggplot2::scale_size_area(max_size = 20) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = paste0(
        origin_name,
        " genes harboring P/LP variants"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}
