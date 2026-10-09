#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/featureCounts_%j.out
#SBATCH --error=logs/featureCounts_%j.err

###############################################################################
# 04_featureCounts.slurm
#
# Generates a reverse-stranded gene-level count matrix from STAR BAM files.
#
# Expected settings for the supplied A. nidulans GFF3 annotation:
#
#   feature type: exon
#   gene attribute: gene_id
#   strandedness: 2 (reverse stranded)
#
# Submit from the project root:
#
#   sbatch 04_featureCounts.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: featureCounts failed near line ${LINENO}." >&2' ERR

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
module load Subread

command -v "${FEATURECOUNTS_EXECUTABLE}" >/dev/null 2>&1 \
    || fail "featureCounts executable not found after loading Subread: ${FEATURECOUNTS_EXECUTABLE}"

###############################################################################
# REFERENCE AND INPUT VALIDATION
###############################################################################

detect_reference_files

[[ -s "${ANNOTATION_GFF}" ]] \
    || fail "Annotation file is missing or empty: ${ANNOTATION_GFF}"

[[ -s "${STAR_MANIFEST}" ]] \
    || fail "STAR manifest is missing or empty: ${STAR_MANIFEST}"

[[ -s "${SAMPLE_COUNT_FILE}" ]] \
    || fail "Sample-count file is missing or empty: ${SAMPLE_COUNT_FILE}"

[[ -d "${STAR_DIR}" ]] \
    || fail "STAR results directory is missing: ${STAR_DIR}"

SAMPLE_COUNT=$(tr -d '[:space:]' < "${SAMPLE_COUNT_FILE}")

if [[ ! "${SAMPLE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "Invalid sample count in ${SAMPLE_COUNT_FILE}: ${SAMPLE_COUNT}"
fi

###############################################################################
# ANNOTATION VALIDATION
###############################################################################

echo
echo "========================================================="
echo "FEATURECOUNTS"
echo "========================================================="
echo
echo "Started:"
date
echo
echo "featureCounts version:"
"${FEATURECOUNTS_EXECUTABLE}" -v
echo
echo "Validating annotation structure..."

# Count all records matching the configured feature type. This avoids
# grep/head pipelines, which can produce broken-pipe failures under pipefail.
FEATURE_COUNT=$(
    awk -F '\t' \
        -v feature_type="${FEATURECOUNTS_FEATURE_TYPE}" '
        BEGIN {
            count = 0
        }

        $0 !~ /^#/ && NF >= 9 && $3 == feature_type {
            count++
        }

        END {
            print count
        }
        ' "${ANNOTATION_GFF}"
)

if [[ ! "${FEATURE_COUNT}" =~ ^[0-9]+$ ]]; then
    fail "Could not determine the number of ${FEATURECOUNTS_FEATURE_TYPE} records in ${ANNOTATION_GFF}."
fi

if [[ "${FEATURE_COUNT}" -eq 0 ]]; then
    fail "No ${FEATURECOUNTS_FEATURE_TYPE} records were found in ${ANNOTATION_GFF}. Check FEATURECOUNTS_FEATURE_TYPE."
fi

# Count configured features whose ninth column contains the gene attribute.
# This recognizes both GFF3 syntax (gene_id=AN0005) and GTF syntax
# (gene_id \"AN0005\";), although this workflow counts directly from the GFF3.
ATTRIBUTE_COUNT=$(
    awk -F '\t' \
        -v feature_type="${FEATURECOUNTS_FEATURE_TYPE}" \
        -v attribute="${FEATURECOUNTS_GENE_ATTRIBUTE}" '
        BEGIN {
            count = 0
        }

        $0 !~ /^#/ && NF >= 9 && $3 == feature_type {
            gff3_pattern = "(^|;)" attribute "="
            gtf_pattern = "(^|;[[:space:]]*)" attribute "[[:space:]]+\""

            if ($9 ~ gff3_pattern || $9 ~ gtf_pattern) {
                count++
            }
        }

        END {
            print count
        }
        ' "${ANNOTATION_GFF}"
)

if [[ ! "${ATTRIBUTE_COUNT}" =~ ^[0-9]+$ ]]; then
    fail "Could not determine ${FEATURECOUNTS_GENE_ATTRIBUTE} attribute counts in ${ANNOTATION_GFF}."
fi

if [[ "${ATTRIBUTE_COUNT}" -eq 0 ]]; then
    fail "The ${FEATURECOUNTS_GENE_ATTRIBUTE} attribute was not found on any ${FEATURECOUNTS_FEATURE_TYPE} records."
fi

if [[ "${ATTRIBUTE_COUNT}" -lt "${FEATURE_COUNT}" ]]; then
    MISSING_ATTRIBUTE_COUNT=$((FEATURE_COUNT - ATTRIBUTE_COUNT))

    echo
echo "WARNING:"
echo "  ${MISSING_ATTRIBUTE_COUNT} of ${FEATURE_COUNT} ${FEATURECOUNTS_FEATURE_TYPE} records"
echo "  do not contain ${FEATURECOUNTS_GENE_ATTRIBUTE}."
echo "  featureCounts will not assign those records using the configured attribute."
echo
fi

echo "Annotation validation passed."
echo "  Feature records: ${FEATURE_COUNT}"
echo "  Records with ${FEATURECOUNTS_GENE_ATTRIBUTE}: ${ATTRIBUTE_COUNT}"
echo

###############################################################################
# COLLECT BAM FILES IN MANIFEST ORDER
###############################################################################

echo "Collecting BAM files in sample-manifest order..."

BAM_FILES=()
SAMPLE_NAMES=()

# Skip the manifest header. The first column contains the standardized sample ID.
while IFS=$'\t' read -r sample_name ignored_fastq; do
    [[ -n "${sample_name}" ]] || continue

    if [[ ! "${sample_name}" =~ ^[A-Za-z0-9._-]+$ ]]; then
        fail "Unsafe sample name in STAR manifest: ${sample_name}"
    fi

    bam_file="${STAR_DIR}/${sample_name}/${sample_name}.Aligned.sortedByCoord.out.bam"
    bam_index="${bam_file}.bai"
    completion_flag="${STAR_DIR}/${sample_name}/${sample_name}.alignment_complete.flag"

    [[ -s "${bam_file}" ]] \
        || fail "BAM file is missing or empty for ${sample_name}: ${bam_file}"

    [[ -s "${bam_index}" ]] \
        || fail "BAM index is missing or empty for ${sample_name}: ${bam_index}"

    # Completion flags are created by the revised STAR script. To support BAMs
    # from the already-successful earlier run, absence is a warning, not fatal.
    if [[ ! -s "${completion_flag}" ]]; then
        echo "WARNING: Alignment completion flag not found for ${sample_name}."
    fi

    SAMPLE_NAMES+=("${sample_name}")
    BAM_FILES+=("${bam_file}")
done < <(tail -n +2 "${STAR_MANIFEST}")

if [[ "${#BAM_FILES[@]}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} BAM files from the manifest but collected ${#BAM_FILES[@]}."
fi

if [[ "${#SAMPLE_NAMES[@]}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} sample names but collected ${#SAMPLE_NAMES[@]}."
fi

echo "Collected ${#BAM_FILES[@]} BAM files."
echo

###############################################################################
# OUTPUT PATHS
###############################################################################

mkdir -p "${COUNTS_DIR}"

COUNTS_FILE="${COUNTS_DIR}/gene_counts.txt"
SUMMARY_FILE="${COUNTS_FILE}.summary"
RUN_RECORD="${COUNTS_DIR}/featureCounts_run_record.txt"

# Remove only outputs that this job recreates. This prevents stale outputs from
# being mistaken for a successful new run.
rm -f "${COUNTS_FILE}" "${SUMMARY_FILE}" "${RUN_RECORD}"

###############################################################################
# RUN INFORMATION
###############################################################################

echo "Annotation:"
echo "  ${ANNOTATION_GFF}"
echo
echo "BAM files:"
echo "  ${#BAM_FILES[@]}"
echo
echo "Feature type:"
echo "  ${FEATURECOUNTS_FEATURE_TYPE}"
echo
echo "Gene attribute:"
echo "  ${FEATURECOUNTS_GENE_ATTRIBUTE}"
echo
echo "Strandedness:"
echo "  ${FEATURECOUNTS_STRAND}"
echo
echo "Threads:"
echo "  ${SLURM_CPUS_PER_TASK}"
echo

###############################################################################
# RUN FEATURECOUNTS
###############################################################################

# FEATURECOUNTS_EXTRA_ARGS is intentionally expanded without quotes because
# config.sh may store multiple command-line arguments in this variable.
# shellcheck disable=SC2086
"${FEATURECOUNTS_EXECUTABLE}" \
    -T "${SLURM_CPUS_PER_TASK}" \
    -s "${FEATURECOUNTS_STRAND}" \
    -t "${FEATURECOUNTS_FEATURE_TYPE}" \
    -g "${FEATURECOUNTS_GENE_ATTRIBUTE}" \
    -a "${ANNOTATION_GFF}" \
    -o "${COUNTS_FILE}" \
    ${FEATURECOUNTS_EXTRA_ARGS} \
    "${BAM_FILES[@]}"

###############################################################################
# OUTPUT VALIDATION
###############################################################################

[[ -s "${COUNTS_FILE}" ]] \
    || fail "Count matrix was not generated or is empty: ${COUNTS_FILE}"

[[ -s "${SUMMARY_FILE}" ]] \
    || fail "featureCounts summary was not generated or is empty: ${SUMMARY_FILE}"

# featureCounts output has six annotation columns followed by one column per BAM.
HEADER_COLUMN_COUNT=$(
    awk -F '\t' '$0 !~ /^#/ {print NF; exit}' "${COUNTS_FILE}"
)

EXPECTED_COLUMN_COUNT=$((SAMPLE_COUNT + 6))

if [[ "${HEADER_COLUMN_COUNT}" -ne "${EXPECTED_COLUMN_COUNT}" ]]; then
    fail "Expected ${EXPECTED_COLUMN_COUNT} count-matrix columns but found ${HEADER_COLUMN_COUNT}."
fi

# The first non-comment line is the header. Count only subsequent gene rows.
GENE_COUNT=$(
    awk '$0 !~ /^#/ {
        if (seen_header) {
            count++
        } else {
            seen_header = 1
        }
    }
    END {
        print count + 0
    }' "${COUNTS_FILE}"
)

if [[ ! "${GENE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
    fail "No gene rows were detected in ${COUNTS_FILE}."
fi

SUMMARY_SAMPLE_COUNT=$(
    awk -F '\t' 'NR == 1 {print NF - 1; exit}' "${SUMMARY_FILE}"
)

if [[ "${SUMMARY_SAMPLE_COUNT}" -ne "${SAMPLE_COUNT}" ]]; then
    fail "Expected ${SAMPLE_COUNT} samples in the featureCounts summary but found ${SUMMARY_SAMPLE_COUNT}."
fi

###############################################################################
# RUN RECORD
###############################################################################

{
    echo "featurecounts_complete=true"
    echo "job_id=${SLURM_JOB_ID:-unknown}"
    echo "samples=${SAMPLE_COUNT}"
    echo "gene_rows=${GENE_COUNT}"
    echo "annotation=${ANNOTATION_GFF}"
    echo "feature_type=${FEATURECOUNTS_FEATURE_TYPE}"
    echo "gene_attribute=${FEATURECOUNTS_GENE_ATTRIBUTE}"
    echo "strandedness=${FEATURECOUNTS_STRAND}"
    echo "counts_file=${COUNTS_FILE}"
    echo "summary_file=${SUMMARY_FILE}"
    echo "completed=$(date --iso-8601=seconds)"
} > "${RUN_RECORD}"

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "Output files:"
ls -lh "${COUNTS_FILE}" "${SUMMARY_FILE}" "${RUN_RECORD}"
echo
echo "Gene rows:"
echo "  ${GENE_COUNT}"
echo
echo "Sample columns:"
echo "  ${SAMPLE_COUNT}"
echo
echo "featureCounts assignment summary:"
column -t -s $'\t' "${SUMMARY_FILE}" || cat "${SUMMARY_FILE}"
echo
echo "========================================================="
echo "FEATURECOUNTS SUCCESSFUL"
echo "========================================================="
echo
echo "Finished:"
date
echo
