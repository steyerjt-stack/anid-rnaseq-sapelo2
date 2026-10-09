#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/star_align_%A_%a.out
#SBATCH --error=logs/star_align_%A_%a.err

###############################################################################
# 03_star_align_array.slurm
#
# Aligns fastp-trimmed, single-end RNA-seq reads with STAR.
# Produces coordinate-sorted BAM files and BAM indexes.
#
# Submit from the project root:
#
#   SAMPLE_COUNT=$(cat manifests/sample_count.txt)
#   sbatch --array=1-"${SAMPLE_COUNT}" 03_star_align_array.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: STAR alignment failed near line ${LINENO}." >&2' ERR

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

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || fail "SLURM_ARRAY_TASK_ID is undefined. Submit this script with --array."

if [[ ! "${SLURM_ARRAY_TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
fi

###############################################################################
# MODULES
###############################################################################

module purge

# Replace these with full versioned module names if required by Sapelo2.
module load STAR
module load SAMtools

command -v "${STAR_EXECUTABLE}" >/dev/null 2>&1 \
    || fail "STAR executable not found after loading the STAR module: ${STAR_EXECUTABLE}"

command -v samtools >/dev/null 2>&1 \
    || fail "samtools executable not found after loading the SAMtools module."

###############################################################################
# INPUT VALIDATION
###############################################################################

[[ -s "${STAR_MANIFEST}" ]] || fail "STAR manifest is missing or empty: ${STAR_MANIFEST}"
[[ -s "${SAMPLE_COUNT_FILE}" ]] || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"
[[ -d "${STAR_INDEX_DIR}" ]] || fail "STAR index directory is missing: ${STAR_INDEX_DIR}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

if (( SLURM_ARRAY_TASK_ID > SAMPLE_COUNT )); then
    fail "Array task ${SLURM_ARRAY_TASK_ID} exceeds the sample count ${SAMPLE_COUNT}."
fi

required_index_files=(
    Genome
    SA
    SAindex
    chrName.txt
    genomeParameters.txt
)

for index_file in "${required_index_files[@]}"; do
    [[ -s "${STAR_INDEX_DIR}/${index_file}" ]] \
        || fail "Required STAR index component is missing or empty: ${STAR_INDEX_DIR}/${index_file}"
done

###############################################################################
# GET SAMPLE FROM MANIFEST
###############################################################################

# The STAR manifest contains one header row. Array task 1 corresponds to row 2.
LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))

SAMPLE=$(awk -F '\t' -v line="${LINE_NUMBER}" 'NR == line {print $1}' "${STAR_MANIFEST}")

[[ -n "${SAMPLE}" ]] || fail "No sample was found for array task ${SLURM_ARRAY_TASK_ID}."

if [[ ! "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]]; then
    fail "Unsafe sample name in STAR manifest: ${SAMPLE}"
fi

# Align the fastp-trimmed read file, not the original standardized FASTQ.
FASTQ="${FASTP_DIR}/${SAMPLE}.trimmed.fastq.gz"

[[ -s "${FASTQ}" ]] || fail "Trimmed FASTQ is missing or empty: ${FASTQ}"

###############################################################################
# OUTPUT PATHS
###############################################################################

SAMPLE_DIR="${STAR_DIR}/${SAMPLE}"
PREFIX="${SAMPLE_DIR}/${SAMPLE}."
BAM_FILE="${PREFIX}Aligned.sortedByCoord.out.bam"
BAM_INDEX="${BAM_FILE}.bai"
FINAL_LOG="${PREFIX}Log.final.out"
PROGRESS_LOG="${PREFIX}Log.progress.out"

mkdir -p "${SAMPLE_DIR}"

###############################################################################
# REMOVE INCOMPLETE OUTPUTS FOR THIS SAMPLE
###############################################################################

# This script may be safely resubmitted for one failed array task. Remove only
# outputs belonging to the current sample so stale partial files are not reused.
rm -f \
    "${BAM_FILE}" \
    "${BAM_INDEX}" \
    "${PREFIX}Log.out" \
    "${PREFIX}Log.progress.out" \
    "${PREFIX}Log.final.out" \
    "${PREFIX}SJ.out.tab"

if [[ -d "${PREFIX}_STARtmp" ]]; then
    rm -rf "${PREFIX}_STARtmp"
fi

###############################################################################
# RUN INFORMATION
###############################################################################

echo
echo "========================================================="
echo "STAR ALIGNMENT"
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
echo "Input trimmed FASTQ:"
echo "  ${FASTQ}"
echo
echo "STAR index:"
echo "  ${STAR_INDEX_DIR}"
echo
echo "Output directory:"
echo "  ${SAMPLE_DIR}"
echo
echo "Threads:"
echo "  ${SLURM_CPUS_PER_TASK}"
echo
echo "STAR version:"
"${STAR_EXECUTABLE}" --version
echo
echo "SAMtools version:"
samtools --version | head -n 1
echo

###############################################################################
# RUN STAR
###############################################################################

# STAR_EXTRA_ARGS is intentionally expanded without quotes because config.sh
# stores it as multiple command-line arguments.
# shellcheck disable=SC2086
"${STAR_EXECUTABLE}" \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --readFilesIn "${FASTQ}" \
    --readFilesCommand zcat \
    --twopassMode "${STAR_TWO_PASS_MODE}" \
    --outFileNamePrefix "${PREFIX}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM XS \
    ${STAR_EXTRA_ARGS}

###############################################################################
# VALIDATE AND INDEX BAM
###############################################################################

[[ -s "${BAM_FILE}" ]] || fail "STAR did not produce a nonempty BAM file: ${BAM_FILE}"
[[ -s "${FINAL_LOG}" ]] || fail "STAR final log is missing or empty: ${FINAL_LOG}"

# quickcheck exits nonzero for a truncated or malformed BAM.
samtools quickcheck -v "${BAM_FILE}"

samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM_FILE}"

[[ -s "${BAM_INDEX}" ]] || fail "BAM index was not produced: ${BAM_INDEX}"

###############################################################################
# COMPLETION FLAG
###############################################################################

cat > "${SAMPLE_DIR}/${SAMPLE}.alignment_complete.flag" <<EOF
alignment_complete=true
sample=${SAMPLE}
job_id=${SLURM_JOB_ID:-unknown}
array_task=${SLURM_ARRAY_TASK_ID}
fastq=${FASTQ}
bam=${BAM_FILE}
completed=$(date --iso-8601=seconds)
EOF

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "BAM output:"
ls -lh "${BAM_FILE}" "${BAM_INDEX}"
echo
echo "STAR final alignment summary:"
cat "${FINAL_LOG}"
echo
echo "========================================================="
echo "STAR ALIGNMENT SUCCESSFUL"
echo "========================================================="
echo
echo "Sample:"
echo "  ${SAMPLE}"
echo
echo "BAM:"
echo "  ${BAM_FILE}"
echo
echo "BAM index:"
echo "  ${BAM_INDEX}"
echo
echo "Progress log:"
echo "  ${PROGRESS_LOG}"
echo
echo "Finished:"
date
echo
