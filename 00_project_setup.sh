#!/bin/bash

###############################################################################
# 00_project_setup.sh
#
# Initializes the Aspergillus nidulans RNA-seq project before Slurm jobs run.
#
# Run from the project root:
#
#   bash 00_project_setup.sh
#
# Responsibilities:
#   1. Load config.sh.
#   2. Detect the genome FASTA and GFF/GFF3 annotation.
#   3. Validate metadata and FASTQ files.
#   4. Create the project directory structure.
#   5. Create standardized FASTQ symlinks.
#   6. Generate sample and STAR manifests.
#   7. Generate a STAR-compatible GTF with gffread.
#   8. Perform final setup sanity checks.
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: Setup failed near line ${LINENO}." >&2' ERR

fail() {
    echo >&2
    echo "ERROR: $*" >&2
    echo >&2
    exit 1
}

###############################################################################
# PROJECT ROOT AND CONFIGURATION
###############################################################################

# Resolve paths relative to this script so setup can be launched from another
# directory without silently writing outputs to the wrong location.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}"

[[ -f "config.sh" ]] || fail "config.sh was not found in ${SCRIPT_DIR}."

# shellcheck source=/dev/null
source "config.sh"

###############################################################################
# HEADER
###############################################################################

echo
echo "========================================================="
echo "Aspergillus nidulans RNA-seq Pipeline"
echo "Project Setup"
echo "========================================================="
echo
echo "Project root:"
echo "  ${SCRIPT_DIR}"
echo

###############################################################################
# REQUIRED COMMANDS
###############################################################################

if ! command -v Rscript >/dev/null 2>&1; then
    fail "Rscript was not found. Load the same versioned R module used by the pipeline, then rerun setup. Example: module spider R; module load R/VERSION"
fi

if ! command -v realpath >/dev/null 2>&1; then
    fail "The realpath command is required but was not found."
fi

###############################################################################
# CONFIGURATION VALIDATION
###############################################################################

required_config_variables=(
    PROJECT_NAME
    RAW_FASTQ_DIR
    REFERENCE_DIR
    RESULTS_DIR
    FASTP_DIR
    STAR_DIR
    COUNTS_DIR
    DESEQ2_DIR
    MULTIQC_DIR
    QC_DIR
    LOG_DIR
    SAMPLE_KEY
    MANIFEST_DIR
    SAMPLE_MANIFEST
    STAR_MANIFEST
    SAMPLE_COUNT_FILE
    VALIDATION_SUMMARY
    STAR_INDEX_DIR
    GENERATED_GTF
    GFFREAD_EXECUTABLE
)

for variable_name in "${required_config_variables[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        fail "Required config variable ${variable_name} is empty or undefined."
    fi
done

###############################################################################
# REQUIRED INPUTS
###############################################################################

[[ -d "${RAW_FASTQ_DIR}" ]] || fail "FASTQ directory not found: ${RAW_FASTQ_DIR}"
[[ -d "${REFERENCE_DIR}" ]] || fail "Reference directory not found: ${REFERENCE_DIR}"
[[ -f "${SAMPLE_KEY}" ]] || fail "Sample key not found: ${SAMPLE_KEY}"
[[ -f "validate_metadata.R" ]] || fail "validate_metadata.R was not found in the project root."

echo "Detecting reference files..."
detect_reference_files

[[ -s "${GENOME_FASTA}" ]] || fail "Detected genome FASTA is missing or empty: ${GENOME_FASTA}"
[[ -s "${ANNOTATION_GFF}" ]] || fail "Detected annotation is missing or empty: ${ANNOTATION_GFF}"

echo
echo "Genome FASTA:"
echo "  ${GENOME_FASTA}"
echo
echo "Annotation:"
echo "  ${ANNOTATION_GFF}"
echo

###############################################################################
# DIRECTORY STRUCTURE
###############################################################################

echo "Creating project directories..."

mkdir -p \
    "${RESULTS_DIR}" \
    "${FASTP_DIR}" \
    "${STAR_DIR}" \
    "${COUNTS_DIR}" \
    "${DESEQ2_DIR}" \
    "${MULTIQC_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}" \
    "${MANIFEST_DIR}" \
    "standardized_fastq" \
    "${STAR_INDEX_DIR}"

echo "Project directories are ready."
echo

###############################################################################
# METADATA VALIDATION
###############################################################################

echo "Running metadata validation..."

Rscript --vanilla validate_metadata.R \
    "${SAMPLE_KEY}" \
    "${RAW_FASTQ_DIR}"

echo
echo "Metadata validation completed successfully."
echo

###############################################################################
# CREATE STANDARDIZED FASTQ LINKS AND MANIFESTS
###############################################################################

echo "Creating standardized FASTQ links and manifests..."

# R reads the CSV safely. Bash deliberately does not parse CSV because quoted
# fields and embedded commas can make shell CSV parsing unreliable.
Rscript --vanilla - \
    "${SAMPLE_KEY}" \
    "${RAW_FASTQ_DIR}" \
    "${MANIFEST_DIR}" \
    "${SAMPLE_MANIFEST}" \
    "${STAR_MANIFEST}" \
    "${SAMPLE_COUNT_FILE}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
    stop("Internal setup error: expected six manifest-generation arguments.")
}

sample_key_file <- args[1]
raw_fastq_dir <- args[2]
manifest_dir <- args[3]
sample_manifest_file <- args[4]
star_manifest_file <- args[5]
sample_count_file <- args[6]

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

metadata <- read_csv(
    sample_key_file,
    show_col_types = FALSE,
    trim_ws = TRUE
)

required_columns <- c(
    "fastq_name",
    "sample",
    "genotype",
    "nitrogen",
    "replicate",
    "batch"
)

missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
    stop(
        "Missing sample-key columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

if (nrow(metadata) == 0) {
    stop("The sample key contains no samples.")
}

if (anyDuplicated(metadata$sample)) {
    stop("Sample names are duplicated in the sample key.")
}

if (anyDuplicated(metadata$fastq_name)) {
    stop("FASTQ identifiers are duplicated in the sample key.")
}

dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("standardized_fastq", recursive = TRUE, showWarnings = FALSE)

resolved_fastqs <- character(nrow(metadata))

for (i in seq_len(nrow(metadata))) {
    fastq_id <- metadata$fastq_name[i]
    sample_id <- metadata$sample[i]

    candidates <- c(
        file.path(raw_fastq_dir, paste0(fastq_id, ".fastq.gz")),
        file.path(raw_fastq_dir, paste0(fastq_id, ".fq.gz"))
    )

    existing_candidates <- candidates[file.exists(candidates)]

    if (length(existing_candidates) == 0) {
        stop("No FASTQ found for sample ", sample_id, ": ", fastq_id)
    }

    if (length(existing_candidates) > 1) {
        stop(
            "Multiple FASTQ extensions found for sample ", sample_id,
            ". Keep only one of .fastq.gz or .fq.gz."
        )
    }

    source_file <- normalizePath(existing_candidates[1], mustWork = TRUE)
    target_file <- file.path(
        "standardized_fastq",
        paste0(sample_id, ".fastq.gz")
    )

    if (file.exists(target_file) || Sys.readlink(target_file) != "") {
        unlink(target_file)
    }

    link_created <- file.symlink(source_file, target_file)

    if (!isTRUE(link_created) || !file.exists(target_file)) {
        stop("Could not create FASTQ symlink: ", target_file)
    }

    resolved_fastqs[i] <- target_file
}

sample_manifest <- metadata |>
    select(sample, genotype, nitrogen, replicate, batch)

star_manifest <- tibble(
    sample = metadata$sample,
    fastq = resolved_fastqs
)

write_tsv(sample_manifest, sample_manifest_file)
write_tsv(star_manifest, star_manifest_file)
writeLines(as.character(nrow(metadata)), sample_count_file)

cat("Created ", nrow(metadata), " standardized FASTQ links.\n", sep = "")
cat("Wrote sample manifest: ", sample_manifest_file, "\n", sep = "")
cat("Wrote STAR manifest: ", star_manifest_file, "\n", sep = "")
cat("Wrote sample count: ", sample_count_file, "\n", sep = "")
RSCRIPT

echo
echo "FASTQ links and manifests created successfully."
echo

###############################################################################
# PREPARE STAR-COMPATIBLE GTF
###############################################################################

echo "Preparing STAR-compatible annotation..."

GTF_FILE="${GENERATED_GTF}"

if [[ -s "${GTF_FILE}" ]]; then
    echo "Existing nonempty GTF detected:"
    echo "  ${GTF_FILE}"
else
    if ! command -v "${GFFREAD_EXECUTABLE}" >/dev/null 2>&1; then
        fail "${GFFREAD_EXECUTABLE} was not found. Load the versioned gffread module or create ${GTF_FILE} manually, then rerun setup."
    fi

    echo "Generating GTF with ${GFFREAD_EXECUTABLE}..."

    "${GFFREAD_EXECUTABLE}" \
        "${ANNOTATION_GFF}" \
        -T \
        -o "${GTF_FILE}"

    [[ -s "${GTF_FILE}" ]] || fail "gffread completed but the GTF is missing or empty: ${GTF_FILE}"

    echo "Generated GTF:"
    echo "  ${GTF_FILE}"
fi

echo

###############################################################################
# REFERENCE CONSISTENCY CHECKS
###############################################################################

echo "Checking reference sequence-name compatibility..."

FASTA_NAMES_FILE=$(mktemp)
GTF_NAMES_FILE=$(mktemp)
trap 'rm -f "${FASTA_NAMES_FILE:-}" "${GTF_NAMES_FILE:-}"' EXIT

awk '/^>/ {
    name = substr($0, 2)
    sub(/[[:space:]].*$/, "", name)
    print name
}' "${GENOME_FASTA}" | sort -u > "${FASTA_NAMES_FILE}"

awk -F '\t' '$0 !~ /^#/ && NF >= 9 {print $1}' \
    "${GTF_FILE}" | sort -u > "${GTF_NAMES_FILE}"

[[ -s "${FASTA_NAMES_FILE}" ]] || fail "No sequence names were found in the genome FASTA."
[[ -s "${GTF_NAMES_FILE}" ]] || fail "No sequence names were found in the generated GTF."

MISSING_GTF_SEQUENCES=$(comm -23 "${GTF_NAMES_FILE}" "${FASTA_NAMES_FILE}" || true)

if [[ -n "${MISSING_GTF_SEQUENCES}" ]]; then
    echo >&2
    echo "ERROR: GTF sequence names absent from the genome FASTA:" >&2
    echo "${MISSING_GTF_SEQUENCES}" >&2
    echo >&2
    exit 1
fi

echo "Reference sequence names are compatible."
echo

###############################################################################
# FINAL SETUP VALIDATION
###############################################################################

[[ -s "${SAMPLE_MANIFEST}" ]] || fail "Sample manifest is missing or empty: ${SAMPLE_MANIFEST}"
[[ -s "${STAR_MANIFEST}" ]] || fail "STAR manifest is missing or empty: ${STAR_MANIFEST}"
[[ -s "${SAMPLE_COUNT_FILE}" ]] || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"
[[ -s "${GTF_FILE}" ]] || fail "Generated GTF is missing or empty: ${GTF_FILE}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count recorded in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

LINK_COUNT=$(find standardized_fastq -maxdepth 1 -type l -name "*.fastq.gz" | wc -l)
STAR_MANIFEST_ROWS=$(awk 'END {print NR - 1}' "${STAR_MANIFEST}")
SAMPLE_MANIFEST_ROWS=$(awk 'END {print NR - 1}' "${SAMPLE_MANIFEST}")

if [[ "${LINK_COUNT}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} standardized FASTQ links but found ${LINK_COUNT}."
fi

if [[ "${STAR_MANIFEST_ROWS}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} STAR manifest records but found ${STAR_MANIFEST_ROWS}."
fi

if [[ "${SAMPLE_MANIFEST_ROWS}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} sample manifest records but found ${SAMPLE_MANIFEST_ROWS}."
fi

# Record a setup-completion flag only after all checks pass.
cat > "${MANIFEST_DIR}/setup_complete.flag" <<EOF
setup_complete=true
project=${PROJECT_NAME}
samples=${SAMPLE_COUNT}
genome=${GENOME_FASTA}
annotation_gff=${ANNOTATION_GFF}
annotation_gtf=${GTF_FILE}
completed=$(date --iso-8601=seconds)
EOF

###############################################################################
# FINAL SUMMARY
###############################################################################

echo "========================================================="
echo "SETUP COMPLETE"
echo "========================================================="
echo
echo "Project:"
echo "  ${PROJECT_NAME}"
echo
echo "Biological samples:"
echo "  ${SAMPLE_COUNT}"
echo
echo "Genome FASTA:"
echo "  ${GENOME_FASTA}"
echo
echo "Annotation GFF/GFF3:"
echo "  ${ANNOTATION_GFF}"
echo
echo "STAR-compatible GTF:"
echo "  ${GTF_FILE}"
echo
echo "Sample manifest:"
echo "  ${SAMPLE_MANIFEST}"
echo
echo "STAR manifest:"
echo "  ${STAR_MANIFEST}"
echo
echo "Standardized FASTQ links:"
echo "  standardized_fastq/"
echo
echo "Setup flag:"
echo "  ${MANIFEST_DIR}/setup_complete.flag"
echo
echo "Next step:"
echo "  SAMPLE_COUNT=\$(cat ${SAMPLE_COUNT_FILE})"
echo "  sbatch --array=1-\"\${SAMPLE_COUNT}\" 01_fastp.slurm"
echo
