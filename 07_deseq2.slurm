#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/deseq2_%j.out
#SBATCH --error=logs/deseq2_%j.err

###############################################################################
# 07_deseq2.slurm
#
# Runs the final differential-expression analysis using the featureCounts matrix
# and biological-sample metadata.
#
# Inputs:
#   metadata/sample_key.csv
#   results/counts/gene_counts.txt
#   06_deseq2_analysis.R
#
# Expected output:
#   results/deseq2/analysis_complete.flag
#
# Submit from the project root after featureCounts completes:
#
#   sbatch 07_deseq2.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: DESeq2 job failed near line ${LINENO}." >&2' ERR

fail() {
    echo >&2
    echo "ERROR: $*" >&2
    echo >&2
    exit 1
}

###############################################################################
# PROJECT ROOT AND CONFIGURATION
###############################################################################

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "${SCRIPT_DIR}"

[[ -f "config.sh" ]] || fail "config.sh was not found in ${SCRIPT_DIR}."

# shellcheck source=/dev/null
source "config.sh"

###############################################################################
# MODULES
###############################################################################

module purge

# Replace this with the exact versioned R module under which the required user
# packages were installed. Using the same R version is essential because user
# libraries are normally specific to the R minor version.
module load R

command -v Rscript >/dev/null 2>&1 \
    || fail "Rscript was not found after loading the R module."

###############################################################################
# INPUT PATHS
###############################################################################

METADATA_FILE="${SAMPLE_KEY}"
COUNT_MATRIX="${COUNTS_DIR}/gene_counts.txt"
ANALYSIS_SCRIPT="06_deseq2_analysis.R"
COUNTS_SUMMARY="${COUNT_MATRIX}.summary"

[[ -s "${METADATA_FILE}" ]] \
    || fail "Metadata file is missing or empty: ${METADATA_FILE}"

[[ -s "${COUNT_MATRIX}" ]] \
    || fail "Count matrix is missing or empty: ${COUNT_MATRIX}"

[[ -s "${COUNTS_SUMMARY}" ]] \
    || fail "featureCounts summary is missing or empty: ${COUNTS_SUMMARY}"

[[ -s "${ANALYSIS_SCRIPT}" ]] \
    || fail "DESeq2 analysis script is missing or empty: ${ANALYSIS_SCRIPT}"

[[ -s "${SAMPLE_COUNT_FILE}" ]] \
    || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

###############################################################################
# SCRIPT AND PACKAGE VALIDATION
###############################################################################

echo "Checking DESeq2 analysis script syntax..."

Rscript --vanilla -e \
    "invisible(parse(file='${ANALYSIS_SCRIPT}')); cat('R syntax OK\n')"

# The analysis script imports these packages individually. Check them in the
# exact batch-job R environment before performing any expensive model fitting.
Rscript --vanilla - <<'RSCRIPT'
packages <- c(
    "DESeq2",
    "apeglm",
    "dplyr",
    "tidyr",
    "tibble",
    "ggplot2",
    "purrr",
    "readr",
    "pheatmap",
    "vsn",
    "RColorBrewer",
    "ggrepel"
)

missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
    stop(
        "Missing R package(s) in the batch-job environment: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

cat("All required R packages are available.\n")
RSCRIPT

###############################################################################
# COUNT-MATRIX SAMPLE VALIDATION
###############################################################################

COUNT_SAMPLE_COUNT=$(
    awk -F '\t' '$0 !~ /^#/ {print NF - 6; exit}' "${COUNT_MATRIX}"
)

if [[ ! "${COUNT_SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Could not determine the number of samples in ${COUNT_MATRIX}."
fi

if [[ "${COUNT_SAMPLE_COUNT}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Count matrix contains ${COUNT_SAMPLE_COUNT} samples but setup reports ${SAMPLE_COUNT}."
fi

###############################################################################
# OUTPUT DIRECTORY
###############################################################################

mkdir -p "${DESEQ2_DIR}"

SUCCESS_FLAG="${DESEQ2_DIR}/analysis_complete.flag"
RUN_RECORD="${DESEQ2_DIR}/slurm_run_record.txt"

# Remove success markers before starting. If the job fails, an old completion
# flag cannot be mistaken for success. The R script overwrites its output files.
rm -f "${SUCCESS_FLAG}" "${RUN_RECORD}"

###############################################################################
# RUN INFORMATION
###############################################################################

echo
echo "========================================================="
echo "DESEQ2 ANALYSIS"
echo "========================================================="
echo
echo "Started:"
date
echo
echo "Slurm job ID:"
echo "  ${SLURM_JOB_ID:-unknown}"
echo
echo "R version:"
Rscript --version 2>&1
echo
echo "Biological samples:"
echo "  ${SAMPLE_COUNT}"
echo
echo "Metadata:"
echo "  ${METADATA_FILE}"
echo
echo "Count matrix:"
echo "  ${COUNT_MATRIX}"
echo
echo "Analysis script:"
echo "  ${ANALYSIS_SCRIPT}"
echo
echo "Output directory:"
echo "  ${DESEQ2_DIR}"
echo

###############################################################################
# RUN ANALYSIS
###############################################################################

Rscript --vanilla \
    "${ANALYSIS_SCRIPT}" \
    "${METADATA_FILE}" \
    "${COUNT_MATRIX}"

###############################################################################
# REQUIRED OUTPUT VALIDATION
###############################################################################

[[ -s "${SUCCESS_FLAG}" ]] \
    || fail "DESeq2 analysis completed without creating the success flag: ${SUCCESS_FLAG}"

required_outputs=(
    "${DESEQ2_DIR}/comparison_summary.csv"
    "${DESEQ2_DIR}/analysis_report.csv"
    "${DESEQ2_DIR}/raw_counts.csv"
    "${DESEQ2_DIR}/normalized_counts.csv"
    "${DESEQ2_DIR}/vst_counts.csv"
    "${DESEQ2_DIR}/all_comparisons_combined.csv"
    "${DESEQ2_DIR}/sessionInfo.txt"
    "${DESEQ2_DIR}/software_versions.csv"
    "${DESEQ2_DIR}/qc/qc_pca.pdf"
)

for output_file in "${required_outputs[@]}"; do
    [[ -s "${output_file}" ]] \
        || fail "Required DESeq2 output is missing or empty: ${output_file}"
done

RESULT_TABLE_COUNT=$(
    find "${DESEQ2_DIR}/results" \
        -maxdepth 1 \
        -type f \
        -name "*.csv" \
        -size +0c \
        | wc -l
)

MA_PLOT_COUNT=$(
    find "${DESEQ2_DIR}/plots/MA" \
        -maxdepth 1 \
        -type f \
        -name "*.pdf" \
        -size +0c \
        | wc -l
)

VOLCANO_PLOT_COUNT=$(
    find "${DESEQ2_DIR}/plots/volcano" \
        -maxdepth 1 \
        -type f \
        -name "*.pdf" \
        -size +0c \
        | wc -l
)

# This experiment has three mutant-versus-WT comparisons and two custom
# non-WT contrasts. Treat a different output count as a failed final package.
EXPECTED_COMPARISONS=5

if [[ "${RESULT_TABLE_COUNT}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    fail "Expected ${EXPECTED_COMPARISONS} comparison tables but found ${RESULT_TABLE_COUNT}."
fi

if [[ "${MA_PLOT_COUNT}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    fail "Expected ${EXPECTED_COMPARISONS} MA plots but found ${MA_PLOT_COUNT}."
fi

if [[ "${VOLCANO_PLOT_COUNT}" -ne "${EXPECTED_COMPARISONS}" ]]; then
    fail "Expected ${EXPECTED_COMPARISONS} volcano plots but found ${VOLCANO_PLOT_COUNT}."
fi

###############################################################################
# SLURM RUN RECORD
###############################################################################

{
    echo "deseq2_slurm_complete=true"
    echo "job_id=${SLURM_JOB_ID:-unknown}"
    echo "samples=${SAMPLE_COUNT}"
    echo "metadata=${METADATA_FILE}"
    echo "count_matrix=${COUNT_MATRIX}"
    echo "analysis_script=${ANALYSIS_SCRIPT}"
    echo "results_directory=${DESEQ2_DIR}"
    echo "comparison_tables=${RESULT_TABLE_COUNT}"
    echo "ma_plots=${MA_PLOT_COUNT}"
    echo "volcano_plots=${VOLCANO_PLOT_COUNT}"
    echo "completed=$(date --iso-8601=seconds)"
} > "${RUN_RECORD}"

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "DESeq2 output directory size:"
du -sh "${DESEQ2_DIR}"
echo
echo "Comparison tables:"
echo "  ${RESULT_TABLE_COUNT}"
echo
echo "MA plots:"
echo "  ${MA_PLOT_COUNT}"
echo
echo "Volcano plots:"
echo "  ${VOLCANO_PLOT_COUNT}"
echo
echo "Comparison summary:"
column -t -s ',' "${DESEQ2_DIR}/comparison_summary.csv" \
    || cat "${DESEQ2_DIR}/comparison_summary.csv"
echo
echo "========================================================="
echo "DESEQ2 SUCCESSFUL"
echo "========================================================="
echo
echo "Results:"
echo "  ${DESEQ2_DIR}"
echo
echo "Success flag:"
echo "  ${SUCCESS_FLAG}"
echo
echo "Slurm run record:"
echo "  ${RUN_RECORD}"
echo
echo "Finished:"
date
echo
