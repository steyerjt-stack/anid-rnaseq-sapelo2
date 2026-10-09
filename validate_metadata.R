#!/usr/bin/env Rscript

###############################################################################
# validate_metadata.R
#
# Validates RNA-seq sample metadata and confirms that each biological sample
# has exactly one merged, compressed FASTQ file.
#
# Usage:
#
#   Rscript validate_metadata.R metadata/sample_key.csv raw_fastq
#
# Exit status:
#   0 = validation passed
#   1 = validation failed
###############################################################################

required_packages <- c("readr", "dplyr", "stringr", "tibble")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
    stop(
        "Missing required R package(s): ",
        paste(missing_packages, collapse = ", "),
        ". Load the correct R module and install the missing package(s) in your user library.",
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tibble)
})

###############################################################################
# HELPERS
###############################################################################

fail <- function(...) {
    message("")
    message("========================================")
    message("VALIDATION FAILED")
    message("========================================")
    message("")
    stop(..., call. = FALSE)
}

pass <- function(message_text) {
    cat("[PASS] ", message_text, "\n", sep = "")
}

warn <- function(message_text) {
    cat("[WARNING] ", message_text, "\n", sep = "")
}

###############################################################################
# COMMAND-LINE ARGUMENTS
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop(
        paste0(
            "Usage:\n",
            "  Rscript validate_metadata.R metadata/sample_key.csv raw_fastq\n"
        ),
        call. = FALSE
    )
}

metadata_file <- args[1]
fastq_directory <- args[2]

###############################################################################
# SETTINGS
###############################################################################

required_columns <- c(
    "fastq_name",
    "sample",
    "genotype",
    "nitrogen",
    "replicate",
    "batch"
)

allowed_genotypes <- c(
    "WT",
    "leuB_del",
    "leuR_del",
    "double_del"
)

expected_comparison_groups <- allowed_genotypes
safe_identifier_pattern <- "^[A-Za-z0-9._-]+$"
allowed_fastq_suffixes <- c(".fastq.gz", ".fq.gz")

###############################################################################
# HEADER
###############################################################################

cat("\n")
cat("========================================\n")
cat("RNA-seq Metadata Validation\n")
cat("========================================\n\n")
cat("Metadata file:\n  ", metadata_file, "\n", sep = "")
cat("FASTQ directory:\n  ", fastq_directory, "\n\n", sep = "")

###############################################################################
# INPUT PATH VALIDATION
###############################################################################

if (!file.exists(metadata_file)) {
    fail("Metadata file not found: ", metadata_file)
}
pass("Metadata file found")

if (!dir.exists(fastq_directory)) {
    fail("FASTQ directory not found: ", fastq_directory)
}
pass("FASTQ directory found")

###############################################################################
# READ METADATA
###############################################################################

metadata <- tryCatch(
    read_csv(
        metadata_file,
        show_col_types = FALSE,
        trim_ws = TRUE,
        na = c("", "NA")
    ),
    error = function(error_condition) {
        fail("Could not read metadata CSV: ", conditionMessage(error_condition))
    }
)

if (nrow(metadata) == 0) {
    fail("The metadata file contains zero samples.")
}
pass(paste("Metadata contains", nrow(metadata), "biological samples"))

###############################################################################
# COLUMN VALIDATION
###############################################################################

missing_columns <- setdiff(required_columns, names(metadata))

if (length(missing_columns) > 0) {
    fail(
        "Missing required metadata column(s): ",
        paste(missing_columns, collapse = ", ")
    )
}
pass("All required metadata columns are present")

extra_columns <- setdiff(names(metadata), required_columns)
if (length(extra_columns) > 0) {
    warn(
        paste0(
            "Additional metadata column(s) will be retained in the input file but ignored by setup: ",
            paste(extra_columns, collapse = ", ")
        )
    )
}

metadata_required <- metadata |>
    select(all_of(required_columns))

###############################################################################
# MISSING AND WHITESPACE VALUES
###############################################################################

character_columns <- c("fastq_name", "sample", "genotype", "nitrogen")

for (column_name in character_columns) {
    metadata_required[[column_name]] <- trimws(
        as.character(metadata_required[[column_name]])
    )
}

missing_rows <- metadata_required |>
    mutate(.metadata_row = row_number() + 1L) |>
    filter(
        if_any(
            all_of(required_columns),
            ~ is.na(.x) || trimws(as.character(.x)) == ""
        )
    )

if (nrow(missing_rows) > 0) {
    print(missing_rows)
    fail(
        nrow(missing_rows),
        " metadata row(s) contain missing or empty required values."
    )
}
pass("No missing or empty required values")

###############################################################################
# SAFE IDENTIFIERS
###############################################################################

invalid_sample_names <- metadata_required |>
    filter(!str_detect(sample, safe_identifier_pattern)) |>
    distinct(sample)

if (nrow(invalid_sample_names) > 0) {
    print(invalid_sample_names)
    fail(
        "Invalid sample name(s). Use only letters, numbers, periods, underscores, and hyphens."
    )
}
pass("Sample names use shell-safe characters")

invalid_fastq_ids <- metadata_required |>
    filter(!str_detect(fastq_name, safe_identifier_pattern)) |>
    distinct(fastq_name)

if (nrow(invalid_fastq_ids) > 0) {
    print(invalid_fastq_ids)
    fail(
        "Invalid FASTQ identifier(s). Enter the FASTQ basename without .fastq.gz or .fq.gz, using only letters, numbers, periods, underscores, and hyphens."
    )
}
pass("FASTQ identifiers use shell-safe characters")

suffix_in_metadata <- vapply(
    metadata_required$fastq_name,
    function(fastq_id) any(endsWith(fastq_id, allowed_fastq_suffixes)),
    logical(1)
)

if (any(suffix_in_metadata)) {
    fail(
        "The fastq_name column must contain basenames without .fastq.gz or .fq.gz. Invalid value(s): ",
        paste(metadata_required$fastq_name[suffix_in_metadata], collapse = ", ")
    )
}
pass("FASTQ identifiers omit compressed FASTQ suffixes")

###############################################################################
# DUPLICATE VALIDATION
###############################################################################

duplicate_samples <- metadata_required |>
    count(sample, name = "occurrences") |>
    filter(occurrences > 1)

if (nrow(duplicate_samples) > 0) {
    print(duplicate_samples)
    fail(
        "Duplicate biological sample names were detected. Technical sequencing-cell FASTQs must be concatenated before analysis, leaving one metadata row per biological sample."
    )
}
pass("Biological sample names are unique")

duplicate_fastq_ids <- metadata_required |>
    count(fastq_name, name = "occurrences") |>
    filter(occurrences > 1)

if (nrow(duplicate_fastq_ids) > 0) {
    print(duplicate_fastq_ids)
    fail("Duplicate FASTQ identifiers were detected.")
}
pass("FASTQ identifiers are unique")

###############################################################################
# GENOTYPE VALIDATION
###############################################################################

invalid_genotypes <- setdiff(
    unique(metadata_required$genotype),
    allowed_genotypes
)

if (length(invalid_genotypes) > 0) {
    fail(
        "Invalid genotype value(s): ",
        paste(invalid_genotypes, collapse = ", "),
        ". Allowed values: ",
        paste(allowed_genotypes, collapse = ", ")
    )
}
pass("Genotype names are valid")

if (!"WT" %in% metadata_required$genotype) {
    fail("No WT samples were found.")
}
pass("WT samples are present")

###############################################################################
# REPLICATE AND BATCH VALIDATION
###############################################################################

parse_positive_integer <- function(values, column_name) {
    character_values <- trimws(as.character(values))

    if (any(!str_detect(character_values, "^[0-9]+$"))) {
        fail(column_name, " must contain whole numbers only.")
    }

    numeric_values <- suppressWarnings(as.integer(character_values))

    if (any(is.na(numeric_values)) || any(numeric_values < 1)) {
        fail(column_name, " must contain positive integers beginning at 1.")
    }

    numeric_values
}

metadata_required$replicate <- parse_positive_integer(
    metadata_required$replicate,
    "replicate"
)
pass("Replicate values are positive integers")

metadata_required$batch <- parse_positive_integer(
    metadata_required$batch,
    "batch"
)
pass("Batch values are positive integers")

replicate_duplicates <- metadata_required |>
    count(genotype, nitrogen, batch, replicate, name = "occurrences") |>
    filter(occurrences > 1)

if (nrow(replicate_duplicates) > 0) {
    print(replicate_duplicates)
    fail(
        "Duplicate replicate identifiers were found within the same genotype, nitrogen condition, and batch."
    )
}
pass("Replicate identifiers are unique within genotype, condition, and batch")

###############################################################################
# FASTQ RESOLUTION
###############################################################################

resolve_fastq <- function(fastq_id) {
    candidates <- file.path(
        fastq_directory,
        paste0(fastq_id, allowed_fastq_suffixes)
    )

    existing_candidates <- candidates[file.exists(candidates)]

    tibble(
        fastq_name = fastq_id,
        matches = length(existing_candidates),
        resolved_fastq = if (length(existing_candidates) == 1) {
            normalizePath(existing_candidates, mustWork = TRUE)
        } else {
            paste(existing_candidates, collapse = ";")
        }
    )
}

fastq_resolution <- bind_rows(
    lapply(metadata_required$fastq_name, resolve_fastq)
)

missing_fastqs <- fastq_resolution |>
    filter(matches == 0)

if (nrow(missing_fastqs) > 0) {
    print(missing_fastqs |>
              select(fastq_name))
    fail(
        nrow(missing_fastqs),
        " metadata FASTQ identifier(s) could not be resolved to .fastq.gz or .fq.gz files."
    )
}
pass("Every metadata row resolves to a FASTQ file")

ambiguous_fastqs <- fastq_resolution |>
    filter(matches > 1)

if (nrow(ambiguous_fastqs) > 0) {
    print(ambiguous_fastqs)
    fail(
        "Both .fastq.gz and .fq.gz files exist for one or more FASTQ identifiers. Keep exactly one compressed FASTQ per biological sample."
    )
}
pass("Each FASTQ identifier resolves unambiguously")

empty_fastqs <- fastq_resolution |>
    mutate(size_bytes = file.info(resolved_fastq)$size) |>
    filter(is.na(size_bytes) | size_bytes <= 0)

if (nrow(empty_fastqs) > 0) {
    print(empty_fastqs |>
              select(fastq_name, resolved_fastq, size_bytes))
    fail("One or more resolved FASTQ files are empty or unreadable.")
}
pass("All resolved FASTQ files are nonempty")

unreadable_fastqs <- fastq_resolution$resolved_fastq[
    file.access(fastq_resolution$resolved_fastq, mode = 4) != 0
]

if (length(unreadable_fastqs) > 0) {
    fail(
        "One or more FASTQ files are not readable: ",
        paste(unreadable_fastqs, collapse = ", ")
    )
}
pass("All resolved FASTQ files are readable")

###############################################################################
# GROUP AND COMPARISON VALIDATION
###############################################################################

group_summary <- metadata_required |>
    count(genotype, nitrogen, name = "biological_replicates") |>
    arrange(nitrogen, genotype)

cat("\nBiological replicate summary:\n\n")
print(group_summary, n = Inf)

if (sum(metadata_required$genotype == "WT") < 2) {
    fail("At least two WT biological replicates are required.")
}
pass("WT replicate count is sufficient for DESeq2")

low_replication_groups <- group_summary |>
    filter(biological_replicates < 2)

if (nrow(low_replication_groups) > 0) {
    print(low_replication_groups)
    fail("Every genotype and nitrogen group must contain at least two biological replicates.")
}
pass("Every observed genotype-condition group has at least two replicates")

missing_comparison_groups <- setdiff(
    expected_comparison_groups,
    unique(metadata_required$genotype)
)

if (length(missing_comparison_groups) > 0) {
    warn(
        paste0(
            "Planned comparison group(s) are absent and will be skipped: ",
            paste(missing_comparison_groups, collapse = ", ")
        )
    )
} else {
    pass("All planned comparison groups are present")
}

###############################################################################
# WRITE VALIDATION OUTPUTS
###############################################################################

manifest_directory <- "manifests"
dir.create(manifest_directory, recursive = TRUE, showWarnings = FALSE)

validation_summary_path <- file.path(
    manifest_directory,
    "validation_summary.txt"
)

validation_table_path <- file.path(
    manifest_directory,
    "validated_samples.tsv"
)

validated_samples <- metadata_required |>
    left_join(
        fastq_resolution |>
            select(fastq_name, resolved_fastq),
        by = "fastq_name"
    ) |>
    select(
        sample,
        fastq_name,
        resolved_fastq,
        genotype,
        nitrogen,
        replicate,
        batch
    )

write_tsv(validated_samples, validation_table_path)

summary_lines <- c(
    "validation_passed=true",
    paste0("metadata_file=", normalizePath(metadata_file, mustWork = TRUE)),
    paste0("fastq_directory=", normalizePath(fastq_directory, mustWork = TRUE)),
    paste0("samples=", nrow(metadata_required)),
    paste0("genotypes=", length(unique(metadata_required$genotype))),
    paste0("nitrogen_conditions=", length(unique(metadata_required$nitrogen))),
    paste0("batches=", length(unique(metadata_required$batch))),
    paste0("wt_samples=", sum(metadata_required$genotype == "WT")),
    paste0("validation_time=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)

writeLines(summary_lines, validation_summary_path)

###############################################################################
# SUCCESS REPORT
###############################################################################

cat("\n")
cat("========================================\n")
cat("VALIDATION PASSED\n")
cat("========================================\n\n")
cat("Biological samples:    ", nrow(metadata_required), "\n", sep = "")
cat("Genotypes:             ", length(unique(metadata_required$genotype)), "\n", sep = "")
cat("Nitrogen conditions:   ", length(unique(metadata_required$nitrogen)), "\n", sep = "")
cat("Batches:               ", length(unique(metadata_required$batch)), "\n", sep = "")
cat("WT samples:            ", sum(metadata_required$genotype == "WT"), "\n", sep = "")
cat("\nValidation summary:\n  ", validation_summary_path, "\n", sep = "")
cat("Validated sample table:\n  ", validation_table_path, "\n\n", sep = "")

quit(save = "no", status = 0)
