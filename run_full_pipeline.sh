#!/bin/bash

###############################################################################
# run_full_pipeline.sh
#
# Submits the complete A. nidulans RNA-seq workflow to Sapelo2 with Slurm job
# dependencies.
#
# Normal use for a new or clean project:
#
#   bash run_full_pipeline.sh
#
# Intentional rerun when completed or partial outputs already exist:
#
#   bash run_full_pipeline.sh --force
#
# The --force option permits submission but does not delete or archive existing
# results. Individual pipeline stages control which stage-specific outputs are
# replaced. For the safest rerun, create a new project directory or archive the
# existing results directory first.
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: Pipeline submission failed near line ${LINENO}." >&2' ERR

fail() {
    echo >&2
    echo "ERROR: $*" >&2
    echo >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  bash run_full_pipeline.sh
  bash run_full_pipeline.sh --force
  bash run_full_pipeline.sh --help

Options:
  --force  Allow submission when prior pipeline outputs are present.
  --help   Display this help message.
EOF
}

###############################################################################
# ARGUMENTS
###############################################################################

FORCE_RERUN=false

if [[ $# -gt 1 ]]; then
    usage
    fail "Only one optional argument is supported."
fi

if [[ $# -eq 1 ]]; then
    case "$1" in
        --force)
            FORCE_RERUN=true
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            usage
            fail "Unknown argument: $1"
            ;;
    esac
fi

###############################################################################
# PROJECT ROOT AND CONFIGURATION
###############################################################################

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
echo "Sapelo2 Workflow Submission"
echo "========================================================="
echo
echo "Project root:"
echo "  ${SCRIPT_DIR}"
echo
echo "Force rerun:"
echo "  ${FORCE_RERUN}"
echo

###############################################################################
# REQUIRED COMMANDS
###############################################################################

command -v sbatch >/dev/null 2>&1 \
    || fail "sbatch was not found. Run this launcher on a Sapelo2 login node."

command -v Rscript >/dev/null 2>&1 \
    || fail "Rscript was not found. Load the versioned R module used by the pipeline before running this launcher."

command -v sha256sum >/dev/null 2>&1 \
    || fail "sha256sum was not found."

###############################################################################
# REQUIRED FILES
###############################################################################

required_files=(
    config.sh
    validate_metadata.R
    00_project_setup.sh
    01_fastp.slurm
    02_star_index.slurm
    03_star_align_array.slurm
    04_featureCounts.slurm
    05_multiqc.slurm
    06_deseq2_analysis.R
    07_deseq2.slurm
    run_full_pipeline.sh
)

for required_file in "${required_files[@]}"; do
    [[ -s "${required_file}" ]] \
        || fail "Required pipeline file is missing or empty: ${required_file}"
done

###############################################################################
# STATIC SCRIPT VALIDATION
###############################################################################

echo "Checking Bash and Slurm script syntax..."

bash_scripts=(
    config.sh
    00_project_setup.sh
    01_fastp.slurm
    02_star_index.slurm
    03_star_align_array.slurm
    04_featureCounts.slurm
    05_multiqc.slurm
    07_deseq2.slurm
    run_full_pipeline.sh
)

for script_file in "${bash_scripts[@]}"; do
    bash -n "${script_file}"
    echo "  [OK] ${script_file}"
done

echo
echo "Checking R script syntax..."

Rscript --vanilla -e \
    "invisible(parse(file='validate_metadata.R')); cat('  [OK] validate_metadata.R\n')"

Rscript --vanilla -e \
    "invisible(parse(file='06_deseq2_analysis.R')); cat('  [OK] 06_deseq2_analysis.R\n')"

###############################################################################
# BROWSER-COPY CORRUPTION CHECK
###############################################################################

echo
echo "Checking scripts for browser HTML artifacts..."

# Build the entity names from adjacent shell strings so this detector does not
# mistake its own source code for a corrupted HTML entity.
HTML_PATTERN='<strong|</strong>|&g''t;|&l''t;|&a''mp;|<br'

for script_file in "${required_files[@]}"; do
    if grep -nE "${HTML_PATTERN}" "${script_file}" >/dev/null 2>&1; then
        echo >&2
        echo "HTML-like browser artifacts were found in ${script_file}:" >&2
        grep -nE "${HTML_PATTERN}" "${script_file}" >&2 || true
        echo >&2
        exit 1
    fi
done

echo "  [OK] No browser HTML artifacts detected."
echo

###############################################################################
# EXISTING OUTPUT PROTECTION
###############################################################################

existing_outputs=()

output_markers=(
    "${FASTP_DIR}"
    "${STAR_INDEX_DIR}"
    "${STAR_DIR}"
    "${COUNTS_DIR}/gene_counts.txt"
    "${MULTIQC_DIR}/multiqc_report.html"
    "${DESEQ2_DIR}/analysis_complete.flag"
)

for output_marker in "${output_markers[@]}"; do
    if [[ -e "${output_marker}" ]]; then
        # Ignore empty directories created during an unused setup attempt.
        if [[ -d "${output_marker}" ]]; then
            if find "${output_marker}" -mindepth 1 -print -quit | grep -q .; then
                existing_outputs+=("${output_marker}")
            fi
        else
            existing_outputs+=("${output_marker}")
        fi
    fi
done

if [[ ${#existing_outputs[@]} -gt 0 && "${FORCE_RERUN}" != true ]]; then
    echo >&2
    echo "Existing pipeline outputs were detected:" >&2
    printf '  %s\n' "${existing_outputs[@]}" >&2
    echo >&2
    echo "Submission stopped to protect the existing run." >&2
    echo >&2
    echo "Recommended options:" >&2
    echo "  1. Create a new project directory for the rerun." >&2
    echo "  2. Archive and move the existing results, then rerun." >&2
    echo "  3. Rerun intentionally with: bash run_full_pipeline.sh --force" >&2
    echo >&2
    exit 1
fi

if [[ ${#existing_outputs[@]} -gt 0 && "${FORCE_RERUN}" == true ]]; then
    echo "WARNING: Existing outputs are present and --force was supplied."
    printf '  %s\n' "${existing_outputs[@]}"
    echo
    echo "The launcher will submit a rerun. Existing outputs are not archived"
    echo "automatically and may be replaced by individual pipeline stages."
    echo
fi

###############################################################################
# PROJECT SETUP
###############################################################################

echo "Running project setup..."

bash 00_project_setup.sh

[[ -s "${MANIFEST_DIR}/setup_complete.flag" ]] \
    || fail "Setup completed without creating ${MANIFEST_DIR}/setup_complete.flag."

[[ -s "${SAMPLE_COUNT_FILE}" ]] \
    || fail "Sample-count file is missing or empty after setup: ${SAMPLE_COUNT_FILE}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count after setup: ${SAMPLE_COUNT}"
fi

[[ -s "${STAR_MANIFEST}" ]] \
    || fail "STAR manifest is missing or empty after setup: ${STAR_MANIFEST}"

STAR_MANIFEST_ROWS=$(awk 'END {print NR - 1}' "${STAR_MANIFEST}")

if [[ "${STAR_MANIFEST_ROWS}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "STAR manifest contains ${STAR_MANIFEST_ROWS} samples but sample_count.txt reports ${SAMPLE_COUNT}."
fi

echo
echo "Setup validated."
echo "Biological samples: ${SAMPLE_COUNT}"
echo

###############################################################################
# SUBMISSION RECORD DIRECTORY
###############################################################################

mkdir -p "${MANIFEST_DIR}" "${LOG_DIR}"

RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
RUN_ID="${PROJECT_NAME}_${RUN_TIMESTAMP}"
RUN_MANIFEST="${MANIFEST_DIR}/pipeline_run_${RUN_TIMESTAMP}.txt"
CHECKSUM_FILE="${MANIFEST_DIR}/pipeline_scripts_${RUN_TIMESTAMP}.sha256"
LATEST_JOB_FILE="${MANIFEST_DIR}/job_ids.txt"

sha256sum "${required_files[@]}" > "${CHECKSUM_FILE}"

###############################################################################
# JOB SUBMISSION HELPER
###############################################################################

submit_job() {
    local submission_output
    local job_id

    submission_output=$(sbatch --parsable "$@")

    # On some Slurm installations --parsable returns JOBID;CLUSTER.
    job_id=${submission_output%%;*}

    if [[ ! "${job_id}" =~ ^[0-9]+(_[0-9]+)?$ ]]; then
        fail "Could not parse a Slurm job ID from: ${submission_output}"
    fi

    printf '%s\n' "${job_id}"
}

###############################################################################
# SUBMIT FASTP AND STAR INDEX
###############################################################################

echo "Submitting fastp array..."
FASTP_JOB=$(submit_job \
    --array="1-${SAMPLE_COUNT}" \
    01_fastp.slurm)
echo "  fastp job: ${FASTP_JOB}"

echo
echo "Submitting STAR index..."
STAR_INDEX_JOB=$(submit_job 02_star_index.slurm)
echo "  STAR index job: ${STAR_INDEX_JOB}"

###############################################################################
# SUBMIT STAR ALIGNMENT AFTER BOTH PREREQUISITES
###############################################################################

echo
echo "Submitting STAR alignment array..."
STAR_ALIGN_JOB=$(submit_job \
    --dependency="afterok:${FASTP_JOB}:${STAR_INDEX_JOB}" \
    --array="1-${SAMPLE_COUNT}" \
    03_star_align_array.slurm)
echo "  STAR alignment job: ${STAR_ALIGN_JOB}"

###############################################################################
# SUBMIT FEATURECOUNTS AFTER ALIGNMENT ARRAY
###############################################################################

echo
echo "Submitting featureCounts..."
FEATURECOUNTS_JOB=$(submit_job \
    --dependency="afterok:${STAR_ALIGN_JOB}" \
    04_featureCounts.slurm)
echo "  featureCounts job: ${FEATURECOUNTS_JOB}"

###############################################################################
# SUBMIT MULTIQC AND DESEQ2 AFTER FEATURECOUNTS
###############################################################################

echo
echo "Submitting MultiQC..."
MULTIQC_JOB=$(submit_job \
    --dependency="afterok:${FEATURECOUNTS_JOB}" \
    05_multiqc.slurm)
echo "  MultiQC job: ${MULTIQC_JOB}"

echo
echo "Submitting DESeq2..."
DESEQ2_JOB=$(submit_job \
    --dependency="afterok:${FEATURECOUNTS_JOB}" \
    07_deseq2.slurm)
echo "  DESeq2 job: ${DESEQ2_JOB}"

###############################################################################
# WRITE RUN MANIFESTS
###############################################################################

cat > "${RUN_MANIFEST}" <<EOF
pipeline_run_id=${RUN_ID}
project=${PROJECT_NAME}
project_root=${SCRIPT_DIR}
user=${USER:-unknown}
host=$(hostname)
submission_time=$(date --iso-8601=seconds)
force_rerun=${FORCE_RERUN}
samples=${SAMPLE_COUNT}
setup_flag=${MANIFEST_DIR}/setup_complete.flag
script_checksums=${CHECKSUM_FILE}
fastp_job=${FASTP_JOB}
star_index_job=${STAR_INDEX_JOB}
star_align_job=${STAR_ALIGN_JOB}
featurecounts_job=${FEATURECOUNTS_JOB}
multiqc_job=${MULTIQC_JOB}
deseq2_job=${DESEQ2_JOB}
EOF

cat > "${LATEST_JOB_FILE}" <<EOF
RUN_ID=${RUN_ID}
FASTP=${FASTP_JOB}
STAR_INDEX=${STAR_INDEX_JOB}
STAR_ALIGN=${STAR_ALIGN_JOB}
FEATURECOUNTS=${FEATURECOUNTS_JOB}
MULTIQC=${MULTIQC_JOB}
DESEQ2=${DESEQ2_JOB}
RUN_MANIFEST=${RUN_MANIFEST}
SUBMISSION_TIME=$(date --iso-8601=seconds)
EOF

###############################################################################
# FINAL SUMMARY
###############################################################################

echo
echo "========================================================="
echo "PIPELINE SUBMITTED"
echo "========================================================="
echo
echo "Run ID:"
echo "  ${RUN_ID}"
echo
echo "Samples:"
echo "  ${SAMPLE_COUNT}"
echo
echo "Job IDs:"
echo "  fastp:          ${FASTP_JOB}"
echo "  STAR index:     ${STAR_INDEX_JOB}"
echo "  STAR alignment: ${STAR_ALIGN_JOB}"
echo "  featureCounts:  ${FEATURECOUNTS_JOB}"
echo "  MultiQC:        ${MULTIQC_JOB}"
echo "  DESeq2:         ${DESEQ2_JOB}"
echo
echo "Monitor current jobs:"
echo "  squeue -u \"\$USER\""
echo "  sq --me"
echo
echo "Inspect all submitted job records:"
echo "  sacct -j ${FASTP_JOB},${STAR_INDEX_JOB},${STAR_ALIGN_JOB},${FEATURECOUNTS_JOB},${MULTIQC_JOB},${DESEQ2_JOB} --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem"
echo
echo "Cancel this workflow:"
echo "  scancel ${FASTP_JOB} ${STAR_INDEX_JOB} ${STAR_ALIGN_JOB} ${FEATURECOUNTS_JOB} ${MULTIQC_JOB} ${DESEQ2_JOB}"
echo
echo "Run manifest:"
echo "  ${RUN_MANIFEST}"
echo
echo "Latest job IDs:"
echo "  ${LATEST_JOB_FILE}"
echo
echo "Script checksums:"
echo "  ${CHECKSUM_FILE}"
echo
