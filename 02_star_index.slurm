#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/star_index_%j.out
#SBATCH --error=logs/star_index_%j.err

###############################################################################
# 02_star_index.slurm
#
# Generates the STAR genome index for Aspergillus nidulans RNA-seq alignment.
#
# Run after:
#
#   bash 00_project_setup.sh
#
# Submit from the project root:
#
#   sbatch 02_star_index.slurm
###############################################################################

set -euo pipefail

###############################################################################
# ERROR HANDLING
###############################################################################

trap 'echo; echo "ERROR: STAR index generation failed near line ${LINENO}." >&2' ERR

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
module load STAR

command -v "${STAR_EXECUTABLE}" >/dev/null 2>&1 \
    || fail "STAR executable not found after loading the STAR module: ${STAR_EXECUTABLE}"

###############################################################################
# DETECT AND VALIDATE REFERENCES
###############################################################################

detect_reference_files

[[ -s "${GENOME_FASTA}" ]] \
    || fail "Genome FASTA is missing or empty: ${GENOME_FASTA}"

[[ -s "${ANNOTATION_GFF}" ]] \
    || fail "Annotation GFF/GFF3 is missing or empty: ${ANNOTATION_GFF}"

[[ -s "${GENERATED_GTF}" ]] \
    || fail "STAR-compatible GTF is missing or empty: ${GENERATED_GTF}. Run bash 00_project_setup.sh first."

###############################################################################
# CONFIGURATION VALIDATION
###############################################################################

if [[ ! "${STAR_GENOME_SA_INDEX_NBASES}" =~ ^[1-9][0-9]*$ ]]; then
    fail "STAR_GENOME_SA_INDEX_NBASES must be a positive integer: ${STAR_GENOME_SA_INDEX_NBASES}"
fi

if [[ ! "${STAR_SJDB_OVERHANG}" =~ ^[1-9][0-9]*$ ]]; then
    fail "STAR_SJDB_OVERHANG must be a positive integer: ${STAR_SJDB_OVERHANG}"
fi

if [[ ! "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
    fail "SLURM_CPUS_PER_TASK is missing or invalid."
fi

###############################################################################
# REFERENCE SEQUENCE-NAME CHECK
###############################################################################

echo
echo "Checking FASTA and GTF sequence-name compatibility..."

FASTA_NAMES_FILE=$(mktemp)
GTF_NAMES_FILE=$(mktemp)

cleanup_temp_files() {
    rm -f "${FASTA_NAMES_FILE:-}" "${GTF_NAMES_FILE:-}"
}

trap cleanup_temp_files EXIT

awk '/^>/ {
    name = substr($0, 2)
    sub(/[[:space:]].*$/, "", name)
    print name
}' "${GENOME_FASTA}" | sort -u > "${FASTA_NAMES_FILE}"

awk -F '\t' '$0 !~ /^#/ && NF >= 9 {print $1}' \
    "${GENERATED_GTF}" | sort -u > "${GTF_NAMES_FILE}"

[[ -s "${FASTA_NAMES_FILE}" ]] \
    || fail "No sequence names were found in the genome FASTA."

[[ -s "${GTF_NAMES_FILE}" ]] \
    || fail "No sequence names were found in the GTF."

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
# SAFE INDEX-DIRECTORY CLEANUP
###############################################################################

EXPECTED_INDEX_DIR="${REFERENCE_DIR}/star_index"

if [[ "${STAR_INDEX_DIR}" != "${EXPECTED_INDEX_DIR}" ]]; then
    fail "Unexpected STAR_INDEX_DIR: ${STAR_INDEX_DIR}. Expected: ${EXPECTED_INDEX_DIR}"
fi

if [[ -z "${STAR_INDEX_DIR}" || "${STAR_INDEX_DIR}" == "/" || "${STAR_INDEX_DIR}" == "." ]]; then
    fail "Refusing to remove unsafe STAR index path: ${STAR_INDEX_DIR}"
fi

case "${STAR_INDEX_DIR}" in
    "${REFERENCE_DIR}"/*)
        ;;
    *)
        fail "STAR index directory is not located inside ${REFERENCE_DIR}: ${STAR_INDEX_DIR}"
        ;;
esac

echo "Recreating STAR index directory..."
rm -rf -- "${STAR_INDEX_DIR}"
mkdir -p "${STAR_INDEX_DIR}"
echo

###############################################################################
# RUN INFORMATION
###############################################################################

echo "========================================================="
echo "STAR INDEX GENERATION"
echo "========================================================="
echo
echo "Started:"
date
echo
echo "Slurm job ID:"
echo "  ${SLURM_JOB_ID:-unknown}"
echo
echo "Genome FASTA:"
echo "  ${GENOME_FASTA}"
echo
echo "Annotation GFF/GFF3:"
echo "  ${ANNOTATION_GFF}"
echo
echo "STAR-compatible GTF:"
echo "  ${GENERATED_GTF}"
echo
echo "STAR index directory:"
echo "  ${STAR_INDEX_DIR}"
echo
echo "Threads:"
echo "  ${SLURM_CPUS_PER_TASK}"
echo
echo "genomeSAindexNbases:"
echo "  ${STAR_GENOME_SA_INDEX_NBASES}"
echo
echo "sjdbOverhang:"
echo "  ${STAR_SJDB_OVERHANG}"
echo
echo "STAR version:"
"${STAR_EXECUTABLE}" --version
echo

###############################################################################
# GENERATE STAR INDEX
###############################################################################

"${STAR_EXECUTABLE}" \
    --runMode genomeGenerate \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FASTA}" \
    --sjdbGTFfile "${GENERATED_GTF}" \
    --sjdbOverhang "${STAR_SJDB_OVERHANG}" \
    --genomeSAindexNbases "${STAR_GENOME_SA_INDEX_NBASES}"

###############################################################################
# INDEX VALIDATION
###############################################################################

required_index_files=(
    Genome
    SA
    SAindex
    chrLength.txt
    chrName.txt
    chrNameLength.txt
    chrStart.txt
    genomeParameters.txt
)

for index_file in "${required_index_files[@]}"; do
    [[ -s "${STAR_INDEX_DIR}/${index_file}" ]] \
        || fail "Required STAR index component is missing or empty: ${STAR_INDEX_DIR}/${index_file}"
done

# STAR writes this message only after successful genome generation.
STAR_LOG="${STAR_INDEX_DIR}/Log.out"

if [[ -s "${STAR_LOG}" ]]; then
    if ! grep -q "finished successfully" "${STAR_LOG}"; then
        fail "STAR index components exist, but Log.out does not report successful completion."
    fi
fi

###############################################################################
# COMPLETION RECORD
###############################################################################

INDEX_SIZE=$(du -sh "${STAR_INDEX_DIR}" | awk '{print $1}')
COMPLETION_FLAG="${STAR_INDEX_DIR}/star_index_complete.flag"

cat > "${COMPLETION_FLAG}" <<EOF
star_index_complete=true
job_id=${SLURM_JOB_ID:-unknown}
genome_fasta=${GENOME_FASTA}
annotation_gtf=${GENERATED_GTF}
index_directory=${STAR_INDEX_DIR}
genome_sa_index_nbases=${STAR_GENOME_SA_INDEX_NBASES}
sjdb_overhang=${STAR_SJDB_OVERHANG}
threads=${SLURM_CPUS_PER_TASK}
index_size=${INDEX_SIZE}
completed=$(date --iso-8601=seconds)
EOF

###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "STAR index size:"
echo "  ${INDEX_SIZE}"
echo
echo "Validated index components:"
for index_file in "${required_index_files[@]}"; do
    echo "  [OK] ${STAR_INDEX_DIR}/${index_file}"
done

echo
echo "Completion flag:"
echo "  ${COMPLETION_FLAG}"
echo
echo "========================================================="
echo "STAR INDEX SUCCESSFUL"
echo "========================================================="
echo
echo "Finished:"
date
echo
