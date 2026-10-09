#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastp_%A_%a.out
#SBATCH --error=logs/fastp_%A_%a.err

###############################################################################
# 01_fastp.slurm
#
# Trims and quality-checks merged, single-end RNA-seq FASTQ files with fastp.
#
# Submit from the project root:
#
#   SAMPLE_COUNT=$(cat manifests/sample_count.txt)
#   sbatch --array=1-"${SAMPLE_COUNT}" 01_fastp.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: fastp failed near line ${LINENO}." >&2' ERR

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
# SLURM ARRAY VALIDATION
###############################################################################

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] \
    || fail "SLURM_ARRAY_TASK_ID is undefined. Submit this script with --array."

if [[ ! "${SLURM_ARRAY_TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
fi

###############################################################################
# MODULES
###############################################################################

module purge

# Replace this with the full versioned module name if required by Sapelo2.
module load fastp

command -v "${FASTP_EXECUTABLE}" >/dev/null 2>&1 \
    || fail "fastp executable not found after loading the fastp module: ${FASTP_EXECUTABLE}"

###############################################################################
# INPUT VALIDATION
###############################################################################

[[ -s "${STAR_MANIFEST}" ]] \
    || fail "STAR manifest is missing or empty: ${STAR_MANIFEST}"

[[ -s "${SAMPLE_COUNT_FILE}" ]] \
    || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

if (( SLURM_ARRAY_TASK_ID > SAMPLE_COUNT )); then
    fail "Array task ${SLURM_ARRAY_TASK_ID} exceeds the sample count ${SAMPLE_COUNT}."
fi

###############################################################################
# GET SAMPLE FROM MANIFEST
###############################################################################

# The manifest contains one header row. Array task 1 corresponds to row 2.
LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))

IFS=$'\t' read -r SAMPLE INPUT_FASTQ < <(
    awk -F '\t' -v line="${LINE_NUMBER}" 'NR == line {print $1 "\t" $2}' \
        "${STAR_MANIFEST}"
)

[[ -n "${SAMPLE:-}" ]] \
    || fail "No sample was found for array task ${SLURM_ARRAY_TASK_ID}."

[[ -n "${INPUT_FASTQ:-}" ]] \
    || fail "No input FASTQ was found for array task ${SLURM_ARRAY_TASK_ID}."

if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    fail "Unsafe sample name in STAR manifest: ${SAMPLE}"
fi

[[ -s "${INPUT_FASTQ}" ]] \
    || fail "Input FASTQ is missing or empty: ${INPUT_FASTQ}"

###############################################################################
# OUTPUT PATHS
###############################################################################

mkdir -p "${FASTP_DIR}"

OUTPUT_FASTQ="${FASTP_DIR}/${SAMPLE}.trimmed.fastq.gz"
JSON_REPORT="${FASTP_DIR}/${SAMPLE}.fastp.json"
HTML_REPORT="${FASTP_DIR}/${SAMPLE}.fastp.html"
COMPLETION_FLAG="${FASTP_DIR}/${SAMPLE}.fastp_complete.flag"

# Remove only outputs that this array task recreates. This prevents stale files
# from being mistaken for successful outputs after a failed rerun.
rm -f \
    "${OUTPUT_FASTQ}" \
    "${JSON_REPORT}" \
    "${HTML_REPORT}" \
    "${COMPLETION_FLAG}"

###############################################################################
# RUN INFORMATION
###############################################################################

echo
echo "========================================================="
echo "FASTP"
echo "========================================================="
echo
echo "Started:"
date
echo
echo "Slurm job ID:"
echo "  ${SLURM_JOB_ID:-unknown}"
echo
echo "Array task:"
echo "  ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_COUNT}"
echo
echo "Sample:"
echo "  ${SAMPLE}"
echo
echo "Input FASTQ:"
echo "  ${INPUT_FASTQ}"
echo
echo "Output FASTQ:"
echo "  ${OUTPUT_FASTQ}"
echo
echo "Threads:"
echo "  ${SLURM_CPUS_PER_TASK}"
echo
echo "fastp version:"
"${FASTP_EXECUTABLE}" --version
echo

###############################################################################
# RUN FASTP
###############################################################################

# This is single-end data. Do not use paired-end-only options such as
# --detect_adapter_for_pe.
"${FASTP_EXECUTABLE}" \
    --thread "${SLURM_CPUS_PER_TASK}" \
    --in1 "${INPUT_FASTQ}" \
    --out1 "${OUTPUT_FASTQ}" \
    --json "${JSON_REPORT}" \
    --html "${HTML_REPORT}" \
    --qualified_quality_phred "${FASTP_QUALIFIED_QUALITY_PHRED}" \
    --length_required "${FASTP_MINIMUM_LENGTH}" \
    --overrepresentation_analysis

###############################################################################
# OUTPUT VALIDATION
###############################################################################

[[ -s "${OUTPUT_FASTQ}" ]] \
    || fail "Trimmed FASTQ was not generated or is empty: ${OUTPUT_FASTQ}"

[[ -s "${JSON_REPORT}" ]] \
    || fail "fastp JSON report was not generated or is empty: ${JSON_REPORT}"

[[ -s "${HTML_REPORT}" ]] \
    || fail "fastp HTML report was not generated or is empty: ${HTML_REPORT}"

# gzip -t validates all concatenated gzip members and catches truncated output.
gzip -t "${OUTPUT_FASTQ}"

###############################################################################
# COMPLETION FLAG
###############################################################################

cat > "${COMPLETION_FLAG}" <<EOF
fastp_complete=true
sample=${SAMPLE}
job_id=${SLURM_JOB_ID:-unknown}
array_task=${SLURM_ARRAY_TASK_ID}
input_fastq=${INPUT_FASTQ}
output_fastq=${OUTPUT_FASTQ}
json_report=${JSON_REPORT}
html_report=${HTML_REPORT}
completed=$(date --iso-8601=seconds)
EOF

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "Output files:"
ls -lh \
    "${OUTPUT_FASTQ}" \
    "${JSON_REPORT}" \
    "${HTML_REPORT}" \
    "${COMPLETION_FLAG}"
echo
echo "========================================================="
echo "FASTP SUCCESSFUL"
echo "========================================================="
echo
echo "Sample:"
echo "  ${SAMPLE}"
echo
echo "Finished:"
date
echo
