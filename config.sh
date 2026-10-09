#!/bin/bash

###############################################################################
# config.sh
#
# Central configuration for the Aspergillus nidulans RNA-seq pipeline on
# Sapelo2. All shell and Slurm scripts source this file.
###############################################################################

###############################################################################
# PROJECT
###############################################################################

export PROJECT_NAME="A_nidulans_RNAseq"

###############################################################################
# TOOL EXECUTABLES
###############################################################################

# Module loading occurs in the applicable Slurm scripts. These variables permit
# use of nonstandard executable names or full executable paths when necessary.
export GFFREAD_EXECUTABLE="gffread"
export STAR_EXECUTABLE="STAR"
export FASTP_EXECUTABLE="fastp"
export FEATURECOUNTS_EXECUTABLE="featureCounts"
export MULTIQC_EXECUTABLE="multiqc"

###############################################################################
# INPUT DIRECTORIES
###############################################################################

# raw_fastq must contain one merged FASTQ per biological sample.
export RAW_FASTQ_DIR="raw_fastq"
export REFERENCE_DIR="reference"
export METADATA_DIR="metadata"

###############################################################################
# OUTPUT DIRECTORIES
###############################################################################

export RESULTS_DIR="results"
export FASTP_DIR="${RESULTS_DIR}/fastp"
export STAR_DIR="${RESULTS_DIR}/star"
export COUNTS_DIR="${RESULTS_DIR}/counts"
export DESEQ2_DIR="${RESULTS_DIR}/deseq2"
export MULTIQC_DIR="${RESULTS_DIR}/multiqc"
export QC_DIR="${RESULTS_DIR}/qc"
export LOG_DIR="logs"

###############################################################################
# METADATA AND MANIFEST FILES
###############################################################################

export SAMPLE_KEY="${METADATA_DIR}/sample_key.csv"
export MANIFEST_DIR="manifests"
export SAMPLE_MANIFEST="${MANIFEST_DIR}/sample_manifest.tsv"
export STAR_MANIFEST="${MANIFEST_DIR}/star_samples.tsv"
export SAMPLE_COUNT_FILE="${MANIFEST_DIR}/sample_count.txt"
export VALIDATION_SUMMARY="${MANIFEST_DIR}/validation_summary.txt"

###############################################################################
# REFERENCE OUTPUTS
###############################################################################

export GENERATED_GTF="${REFERENCE_DIR}/annotation.gtf"
export STAR_INDEX_DIR="${REFERENCE_DIR}/star_index"

###############################################################################
# FASTQ SETTINGS
###############################################################################

# Used for documentation. Setup currently accepts .fastq.gz and .fq.gz.
export ALLOWED_FASTQ_EXTENSIONS=(
    "fastq.gz"
    "fq.gz"
)

###############################################################################
# FASTP SETTINGS
###############################################################################

export FASTP_QUALIFIED_QUALITY_PHRED=20
export FASTP_MINIMUM_LENGTH=30

###############################################################################
# STAR SETTINGS
###############################################################################

# STAR selected 11 for the successful A. nidulans index build. Using 11 avoids
# the warning produced when 12 was supplied for this approximately 31 Mb genome.
export STAR_GENOME_SA_INDEX_NBASES=11

# Reads are 75 bp; STAR recommends read length minus one for sjdbOverhang.
export STAR_SJDB_OVERHANG=74

export STAR_TWO_PASS_MODE="Basic"

# Expanded intentionally as separate command-line arguments by the STAR script.
export STAR_EXTRA_ARGS="\
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04"

###############################################################################
# FEATURECOUNTS SETTINGS
###############################################################################

# NEBNext Stranded mRNA single-end libraries are reverse stranded.
export FEATURECOUNTS_STRAND=2
export FEATURECOUNTS_FEATURE_TYPE="exon"
export FEATURECOUNTS_GENE_ATTRIBUTE="gene_id"

# Optional extra featureCounts arguments. Leave empty unless intentionally set.
export FEATURECOUNTS_EXTRA_ARGS=""

###############################################################################
# DESEQ2 SETTINGS
###############################################################################

# These values are documented here. The final R script currently uses the same
# values directly so the completed analysis remains reproducible.
export ALPHA=0.05
export LOG2FC_THRESHOLD=1
export MIN_COUNT=10
export MIN_SAMPLES=3
export WT_REFERENCE="WT"

# Bash arrays cannot be exported to child processes. This array is available to
# scripts that source config.sh directly.
CUSTOM_COMPARISONS=(
    "double_del:leuB_del"
    "double_del:leuR_del"
)

###############################################################################
# RESOURCE DOCUMENTATION
###############################################################################

# Slurm resources are intentionally declared in each .slurm header because
# #SBATCH directives cannot use shell variables from this configuration file.
export FASTP_CPUS=4
export FASTP_MEM="8G"
export FASTP_TIME="02:00:00"

export INDEX_CPUS=8
export INDEX_MEM="32G"
export INDEX_TIME="04:00:00"

export ALIGN_CPUS=8
export ALIGN_MEM="16G"
export ALIGN_TIME="04:00:00"

export FC_CPUS=8
export FC_MEM="8G"
export FC_TIME="02:00:00"

export MULTIQC_CPUS=2
export MULTIQC_MEM="4G"
export MULTIQC_TIME="01:00:00"

export DESEQ_CPUS=8
export DESEQ_MEM="32G"
export DESEQ_TIME="04:00:00"

###############################################################################
# REFERENCE FILE AUTO-DETECTION
###############################################################################

detect_reference_files() {
    local fasta_files=()
    local gff_files=()

    if [[ ! -d "${REFERENCE_DIR}" ]]; then
        echo >&2
        echo "ERROR: Reference directory does not exist: ${REFERENCE_DIR}" >&2
        echo >&2
        return 1
    fi

    mapfile -t fasta_files < <(
        find "${REFERENCE_DIR}" \
            -maxdepth 1 \
            -type f \
            \( \
                -iname "*.fa" -o \
                -iname "*.fasta" -o \
                -iname "*.fna" \
            \) \
            -print \
            | sort
    )

    mapfile -t gff_files < <(
        find "${REFERENCE_DIR}" \
            -maxdepth 1 \
            -type f \
            \( \
                -iname "*.gff" -o \
                -iname "*.gff3" \
            \) \
            -print \
            | sort
    )

    if [[ ${#fasta_files[@]} -ne 1 ]]; then
        echo >&2
        echo "ERROR: Expected exactly one FASTA file in ${REFERENCE_DIR}." >&2
        echo "Found ${#fasta_files[@]} candidate(s):" >&2

        if [[ ${#fasta_files[@]} -gt 0 ]]; then
            printf '  %s\n' "${fasta_files[@]}" >&2
        fi

        echo >&2
        return 1
    fi

    if [[ ${#gff_files[@]} -ne 1 ]]; then
        echo >&2
        echo "ERROR: Expected exactly one GFF or GFF3 file in ${REFERENCE_DIR}." >&2
        echo "Found ${#gff_files[@]} candidate(s):" >&2

        if [[ ${#gff_files[@]} -gt 0 ]]; then
            printf '  %s\n' "${gff_files[@]}" >&2
        fi

        echo >&2
        return 1
    fi

    export GENOME_FASTA="${fasta_files[0]}"
    export ANNOTATION_GFF="${gff_files[0]}"
}

###############################################################################
# CONFIGURATION REPORT
###############################################################################

print_config() {
    echo
    echo "========================================="
    echo "Pipeline Configuration"
    echo "========================================="
    echo
    echo "Project:"
    echo "  ${PROJECT_NAME}"
    echo
    echo "Sample key:"
    echo "  ${SAMPLE_KEY}"
    echo
    echo "Raw FASTQ directory:"
    echo "  ${RAW_FASTQ_DIR}"
    echo
    echo "Reference directory:"
    echo "  ${REFERENCE_DIR}"
    echo
    echo "Results directory:"
    echo "  ${RESULTS_DIR}"
    echo
    echo "STAR genomeSAindexNbases:"
    echo "  ${STAR_GENOME_SA_INDEX_NBASES}"
    echo
    echo "featureCounts strandedness:"
    echo "  ${FEATURECOUNTS_STRAND}"
    echo
}
