#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err

###############################################################################
# 05_multiqc.slurm
#
# Aggregates quality-control reports from fastp, STAR, and featureCounts.
#
# Expected output:
#
#   results/multiqc/multiqc_report.html
#
# Submit from the project root after featureCounts completes:
#
#   sbatch 05_multiqc.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: MultiQC failed near line ${LINENO}." >&2' ERR

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

# Replace this with the full versioned module name if required by Sapelo2.
module load MultiQC

command -v "${MULTIQC_EXECUTABLE}" >/dev/null 2>&1 \
    || fail "MultiQC executable not found after loading the MultiQC module: ${MULTIQC_EXECUTABLE}"

###############################################################################
# INPUT VALIDATION
###############################################################################

[[ -d "${FASTP_DIR}" ]] \
    || fail "fastp results directory is missing: ${FASTP_DIR}"

[[ -d "${STAR_DIR}" ]] \
    || fail "STAR results directory is missing: ${STAR_DIR}"

[[ -d "${COUNTS_DIR}" ]] \
    || fail "featureCounts results directory is missing: ${COUNTS_DIR}"

[[ -s "${SAMPLE_COUNT_FILE}" ]] \
    || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

FASTP_JSON_COUNT=$(
    find "${FASTP_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.fastp.json" \
        -size +0c \
        | wc -l
)

STAR_FINAL_LOG_COUNT=$(
    find "${STAR_DIR}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name "*.Log.final.out" \
        -size +0c \
        | wc -l
)

if [[ "${FASTP_JSON_COUNT}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} nonempty fastp JSON reports but found ${FASTP_JSON_COUNT}."
fi

if [[ "${STAR_FINAL_LOG_COUNT}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} nonempty STAR Log.final.out files but found ${STAR_FINAL_LOG_COUNT}."
fi

COUNTS_SUMMARY="${COUNTS_DIR}/gene_counts.txt.summary"

[[ -s "${COUNTS_SUMMARY}" ]] \
    || fail "featureCounts summary is missing or empty: ${COUNTS_SUMMARY}"

###############################################################################
# OUTPUT PATHS
###############################################################################

mkdir -p "${MULTIQC_DIR}"

REPORT_BASENAME="multiqc_report"
REPORT_HTML="${MULTIQC_DIR}/${REPORT_BASENAME}.html"
DATA_DIR="${MULTIQC_DIR}/${REPORT_BASENAME}_data"
RUN_RECORD="${MULTIQC_DIR}/multiqc_run_record.txt"

# Remove outputs recreated by this run so stale reports cannot be mistaken for
# a successful new MultiQC job.
rm -f "${REPORT_HTML}" "${RUN_RECORD}"

if [[ -d "${DATA_DIR}" ]]; then
    rm -rf "${DATA_DIR}"
fi

###############################################################################
# RUN INFORMATION
###############################################################################

echo
echo "========================================================="
echo "MULTIQC"
echo "========================================================="
echo
echo "Started:"
date
echo
echo "Slurm job ID:"
echo "  ${SLURM_JOB_ID:-unknown}"
echo
echo "Biological samples:"
echo "  ${SAMPLE_COUNT}"
echo
echo "fastp JSON reports:"
echo "  ${FASTP_JSON_COUNT}"
echo
echo "STAR final logs:"
echo "  ${STAR_FINAL_LOG_COUNT}"
echo
echo "featureCounts summary:"
echo "  ${COUNTS_SUMMARY}"
echo
echo "Output directory:"
echo "  ${MULTIQC_DIR}"
echo
echo "MultiQC version:"
"${MULTIQC_EXECUTABLE}" --version
echo

###############################################################################
# RUN MULTIQC
###############################################################################

# Scan only the three tool-output directories needed by this workflow. This
# avoids recursively processing old DESeq2 outputs or archived failed runs.
"${MULTIQC_EXECUTABLE}" \
    "${FASTP_DIR}" \
    "${STAR_DIR}" \
    "${COUNTS_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "${REPORT_BASENAME}.html" \
    --force

###############################################################################
# OUTPUT VALIDATION
###############################################################################

[[ -s "${REPORT_HTML}" ]] \
    || fail "MultiQC HTML report was not generated or is empty: ${REPORT_HTML}"

[[ -d "${DATA_DIR}" ]] \
    || fail "MultiQC data directory was not generated: ${DATA_DIR}"

# MultiQC creates a source-tracking file inside the data directory. Its
# presence confirms that input files were parsed rather than only creating an
# empty report shell.
SOURCES_FILE="${DATA_DIR}/multiqc_sources.txt"

[[ -s "${SOURCES_FILE}" ]] \
    || fail "MultiQC source record is missing or empty: ${SOURCES_FILE}"

SOURCE_COUNT=$(awk 'NF > 0 {count++} END {print count + 0}' "${SOURCES_FILE}")

if [[ "${SOURCE_COUNT}" -lt 1 ]]; then
    fail "MultiQC did not record any parsed source files."
fi

###############################################################################
# RUN RECORD
###############################################################################

{
    echo "multiqc_complete=true"
    echo "job_id=${SLURM_JOB_ID:-unknown}"
    echo "samples=${SAMPLE_COUNT}"
    echo "fastp_json_reports=${FASTP_JSON_COUNT}"
    echo "star_final_logs=${STAR_FINAL_LOG_COUNT}"
    echo "featurecounts_summary=${COUNTS_SUMMARY}"
    echo "report_html=${REPORT_HTML}"
    echo "data_directory=${DATA_DIR}"
    echo "source_records=${SOURCE_COUNT}"
    echo "completed=$(date --iso-8601=seconds)"
} > "${RUN_RECORD}"

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "MultiQC outputs:"
ls -lh "${REPORT_HTML}" "${RUN_RECORD}"
echo
echo "MultiQC data directory:"
du -sh "${DATA_DIR}"
echo
echo "Parsed source records:"
echo "  ${SOURCE_COUNT}"
echo
echo "========================================================="
echo "MULTIQC SUCCESSFUL"
echo "========================================================="
echo
echo "Report:"
echo "  ${REPORT_HTML}"
echo
echo "Finished:"
date
echo
