# Aspergillus nidulans RNA-seq Pipeline on Sapelo2

A production RNA-seq workflow using **fastp**, **STAR**, **featureCounts**, **MultiQC**, and **DESeq2** on the University of Georgia GACRC Sapelo2 cluster.

## Experiment Summary

This pipeline was developed for:

- Organism: *Aspergillus nidulans* FGSC A4, genome version 4
- Genome size: approximately 31 Mb
- Sequencing: Illumina NextSeq HO, single-end 75 bp
- Library preparation: NEBNext Stranded mRNA
- Strandedness: reverse stranded
- Biological samples: 15
  - 6 WT
  - 3 `leuB_del`
  - 3 `leuR_del`
  - 3 `double_del`
- Two sequencing-cell FASTQ files per biological sample were concatenated before trimming and alignment

## Workflow

```text
Technical FASTQ files from two sequencing cells
                     |
                     v
       Concatenate by biological sample
                     |
                     v
          Metadata and input validation
                     |
                     v
                   fastp
                     |
                     v
          STAR genome index and alignment
                     |
                     v
               featureCounts
                     |
             +-------+-------+
             |               |
             v               v
          MultiQC          DESeq2
                             |
             +---------------+----------------+
             |               |                |
             v               v                v
            PCA          Heatmaps       DE tables and plots
```

## Pipeline Files

```text
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
README.md
```

The working copies on Sapelo2 are authoritative because they include fixes made during testing.

## Important Final Code Changes

The following files changed during testing and should be preserved from the successful Sapelo2 run.

### `config.sh`

Use:

```bash
export STAR_GENOME_SA_INDEX_NBASES=11
```

STAR automatically selected 11 rather than 12 for this approximately 31 Mb reference genome.

### `00_project_setup.sh`

The final working version:

- accepts `GFFREAD_EXECUTABLE` from `config.sh`
- validates metadata and FASTQ files
- creates standardized sample links
- creates manifests
- converts GFF/GFF3 to GTF for STAR
- contains a valid nested `if`, `else`, and `fi` block for GTF generation

### `03_star_align_array.slurm`

The final working version aligns the fastp-trimmed reads:

```bash
FASTQ="${FASTP_DIR}/${SAMPLE}.trimmed.fastq.gz"
```

It does not align the original untrimmed symlinks.

### `04_featureCounts.slurm`

The final working version:

- validates exon records with `awk`
- validates `gene_id` on exon features
- avoids `grep | head` pipelines under `set -o pipefail`
- counts reverse-stranded libraries with `-s 2`
- uses `-t exon -g gene_id`

### `06_deseq2_analysis.R`

The final clean version:

- imports the featureCounts matrix
- constructs the DESeq2 dataset correctly
- uses `~ genotype` for the present single-condition, single-batch experiment
- filters genes with at least 10 counts in at least 3 samples
- uses explicit nonparallel DESeq2 execution
- reports raw and shrunken log2 fold changes
- uses `apeglm` for genotype-versus-WT coefficient shrinkage
- uses normal shrinkage for non-WT custom contrasts
- classifies genes using FDR and shrunken fold change
- creates PCA, distance, correlation, clustering, and variable-gene heatmaps
- creates MA and labeled volcano plots
- exports top 25 upregulated and top 25 downregulated genes
- records R package and session information
- writes `analysis_complete.flag` after successful completion

## Project Directory Structure

```text
Anid_RNA/
├── config.sh
├── validate_metadata.R
├── 00_project_setup.sh
├── 01_fastp.slurm
├── 02_star_index.slurm
├── 03_star_align_array.slurm
├── 04_featureCounts.slurm
├── 05_multiqc.slurm
├── 06_deseq2_analysis.R
├── 07_deseq2.slurm
├── run_full_pipeline.sh
├── README.md
├── metadata/
│   └── sample_key.csv
├── raw_fastq/
├── reference/
│   ├── genome.fa
│   ├── annotation.gff
│   ├── annotation.gtf
│   └── star_index/
├── standardized_fastq/
├── manifests/
├── logs/
└── results/
    ├── fastp/
    ├── star/
    ├── counts/
    ├── multiqc/
    └── deseq2/
```

## Connect to Sapelo2

From a local terminal:

```bash
ssh YourMyID@sapelo2.gacrc.uga.edu
```

When connecting from off campus, connect to the UGA VPN first.

Use the file-transfer node for uploads and downloads:

```text
xfer.gacrc.uga.edu
```

## Software Modules

Search for current modules before running the workflow:

```bash
module spider fastp
module spider STAR
module spider SAMtools
module spider Subread
module spider MultiQC
module spider gffread
module spider R
```

Use versioned module names when possible. Example:

```bash
module purge
module load R/VERSION
```

Each Slurm script should load the modules needed for that step.

## Required R Packages

The final DESeq2 script requires:

```text
DESeq2
apeglm
dplyr
tidyr
tibble
ggplot2
purrr
readr
pheatmap
vsn
RColorBrewer
ggrepel
```

Check them under the same R module used by the Slurm job:

```bash
module purge
module load R/VERSION

Rscript -e "
packages <- c(
  'DESeq2', 'apeglm', 'dplyr', 'tidyr', 'tibble', 'ggplot2',
  'purrr', 'readr', 'pheatmap', 'vsn', 'RColorBrewer', 'ggrepel'
)
for (package in packages) {
  cat(sprintf('%-15s %s\\n', package,
              requireNamespace(package, quietly = TRUE)))
}
"
```

Install Bioconductor packages without trying to update the cluster-managed library:

```r
BiocManager::install(
  "DESeq2",
  ask = FALSE,
  update = FALSE
)
```

## Reference Files

Place exactly one genome FASTA and one GFF/GFF3 annotation in `reference/`.

Accepted FASTA extensions:

```text
.fa
.fasta
.fna
```

Accepted annotation extensions:

```text
.gff
.gff3
```

The setup script generates:

```text
reference/annotation.gtf
```

using `gffread`.

If needed, generate the GTF manually:

```bash
gffread \
    reference/annotation.gff \
    -T \
    -o reference/annotation.gtf
```

The supplied annotation contains exon records with `gene_id`, allowing featureCounts to use:

```bash
-t exon -g gene_id
```

## Merge Technical FASTQ Files

The two FASTQ files from separate sequencing cells represent the same biological sample and should not be entered as separate DESeq2 samples.

Compressed FASTQ files can be concatenated directly:

```bash
cat cell1_sample.fastq.gz cell2_sample.fastq.gz \
    > merged_fastq/WT_R1.fastq.gz
```

After merging, there should be one FASTQ per biological sample.

Verify the expected count:

```bash
find raw_fastq -maxdepth 1 -name "*.fastq.gz" | wc -l
```

Expected for this experiment:

```text
15
```

## Final Metadata Design

`metadata/sample_key.csv` should contain one row per biological sample.

```csv
fastq_name,sample,genotype,nitrogen,replicate,batch
WT_R1,WT_R1,WT,N1,1,1
WT_R2,WT_R2,WT,N1,2,1
WT_R3,WT_R3,WT,N1,3,1
WT_R4,WT_R4,WT,N1,4,1
WT_R5,WT_R5,WT,N1,5,1
WT_R6,WT_R6,WT,N1,6,1
leuB_R1,leuB_R1,leuB_del,N1,1,1
leuB_R2,leuB_R2,leuB_del,N1,2,1
leuB_R3,leuB_R3,leuB_del,N1,3,1
leuR_R1,leuR_R1,leuR_del,N1,1,1
leuR_R2,leuR_R2,leuR_del,N1,2,1
leuR_R3,leuR_R3,leuR_del,N1,3,1
double_R1,double_R1,double_del,N1,1,1
double_R2,double_R2,double_del,N1,2,1
double_R3,double_R3,double_del,N1,3,1
```

Allowed genotypes:

```text
WT
leuB_del
leuR_del
double_del
```

Sample names may contain letters, numbers, periods, underscores, and hyphens.

## Validate Scripts Before Running

Check Bash and Slurm scripts:

```bash
bash -n 00_project_setup.sh
bash -n 01_fastp.slurm
bash -n 02_star_index.slurm
bash -n 03_star_align_array.slurm
bash -n 04_featureCounts.slurm
bash -n 05_multiqc.slurm
bash -n 07_deseq2.slurm
```

Check the R script:

```bash
Rscript -e "invisible(parse(file='06_deseq2_analysis.R')); cat('R syntax OK\\n')"
```

Check for browser HTML accidentally copied into code:

```bash
grep -nE '<strong|</strong>|&gt;|&lt;|&amp;|<br' \
    06_deseq2_analysis.R
```

The HTML check should return no output.

## Project Setup

Run setup from the project root:

```bash
bash 00_project_setup.sh
```

Successful setup should:

- validate 15 samples
- identify the genome and annotation
- create standardized FASTQ links
- create manifests
- create or detect `reference/annotation.gtf`

Verify:

```bash
cat manifests/sample_count.txt
```

Expected:

```text
15
```

## Educational Step-by-Step Workflow

### 1. fastp

```bash
SAMPLE_COUNT=$(cat manifests/sample_count.txt)

sbatch --array=1-"${SAMPLE_COUNT}" \
    01_fastp.slurm
```

Verify 15 outputs:

```bash
find results/fastp \
    -maxdepth 1 \
    -name "*.trimmed.fastq.gz" \
    | wc -l
```

### 2. STAR index

```bash
sbatch 02_star_index.slurm
```

A successful index contains:

```text
Genome
SA
SAindex
chrName.txt
genomeParameters.txt
```

Verify:

```bash
for file in Genome SA SAindex chrName.txt genomeParameters.txt; do
    if [[ -s "reference/star_index/${file}" ]]; then
        echo "[OK] ${file}"
    else
        echo "[MISSING] ${file}"
    fi
done
```

### 3. STAR alignment

```bash
SAMPLE_COUNT=$(cat manifests/sample_count.txt)

sbatch --array=1-"${SAMPLE_COUNT}" \
    03_star_align_array.slurm
```

Verify 15 BAM files and indexes:

```bash
find results/star \
    -name "*.Aligned.sortedByCoord.out.bam" \
    -size +0c \
    | wc -l
```

```bash
find results/star \
    -name "*.Aligned.sortedByCoord.out.bam.bai" \
    -size +0c \
    | wc -l
```

### 4. featureCounts

```bash
sbatch 04_featureCounts.slurm
```

Expected files:

```text
results/counts/gene_counts.txt
results/counts/gene_counts.txt.summary
```

Verify 15 sample columns:

```bash
awk 'BEGIN {FS="\t"} !/^#/ {
    print "Total columns:", NF
    print "Sample columns:", NF - 6
    exit
}' results/counts/gene_counts.txt
```

### 5. MultiQC

```bash
sbatch 05_multiqc.slurm
```

Expected report:

```text
results/multiqc/multiqc_report.html
```

### 6. DESeq2

```bash
sbatch 07_deseq2.slurm
```

Confirm completion:

```bash
cat results/deseq2/analysis_complete.flag
```

## Automated Workflow

After independently testing the steps, submit the dependency-managed workflow with:

```bash
bash run_full_pipeline.sh
```

The launcher writes job IDs to:

```text
manifests/job_ids.txt
```

## Monitor Jobs

```bash
squeue -u "$USER"
```

Sapelo2's formatted wrapper:

```bash
sq --me
```

Inspect a finished job:

```bash
sacct -j JOBID \
    --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem
```

A successful job should normally report:

```text
State     COMPLETED
ExitCode  0:0
```

Inspect logs:

```bash
cat logs/job_JOBID.err
```

```bash
tail -n 100 logs/job_JOBID.out
```

## DESeq2 Model and Filtering

For the completed experiment, the design resolved to:

```r
~ genotype
```

because only one nitrogen condition and one batch were observed.

Genes were retained when they had:

```text
at least 10 counts in at least 3 samples
```

The successful run loaded 10,988 genes and retained 9,480 after filtering.

## Differential-Expression Comparisons

Automatically generated genotype-versus-WT comparisons:

```text
double_del vs WT
leuB_del vs WT
leuR_del vs WT
```

Custom comparisons:

```text
double_del vs leuB_del
double_del vs leuR_del
```

## Significance Criteria

A gene is classified as significant when:

```text
adjusted p-value < 0.05
```

and:

```text
absolute shrunken log2 fold change >= 1
```

Classification values:

```text
UP
DOWN
NOT_SIGNIFICANT
```

## Fold-Change Reporting

Each full comparison table contains:

```text
gene_id
baseMean
raw_log2FC
shrunken_log2FC
lfc_difference
shrinkage_method
lfcSE
stat
pvalue
padj
neg_log10_fdr
classification
```

Shrinkage methods:

```text
Mutant vs WT: apeglm
Custom non-WT contrasts: normal
```

Use the shrunken fold change for ranking and plots. Preserve the raw fold change for historical comparisons with the older HISAT2 workflow.

## Main DESeq2 Outputs

```text
results/deseq2/raw_counts.csv
results/deseq2/normalized_counts.csv
results/deseq2/vst_counts.csv
results/deseq2/comparison_summary.csv
results/deseq2/all_comparisons_combined.csv
results/deseq2/results/
results/deseq2/significant/
results/deseq2/plots/MA/
results/deseq2/plots/volcano/
results/deseq2/qc/
```

### QC plots

```text
qc_pca.pdf
qc_sample_distance_heatmap.pdf
qc_sample_correlation_heatmap.pdf
qc_sample_clustering.pdf
qc_top100_variable_genes_heatmap.pdf
```

### Per-comparison outputs

- full results table
- significant-only table
- top 25 upregulated genes
- top 25 downregulated genes
- MA plot
- volcano plot with up to 15 upregulated and 15 downregulated labels

## Recommended Scientific Review Order

Review the following before interpreting individual genes:

1. `results/multiqc/multiqc_report.html`
2. STAR mapping summaries
3. featureCounts assignment summary
4. `results/deseq2/qc/qc_pca.pdf`
5. sample-distance and correlation heatmaps
6. `results/deseq2/dispersion_plot.pdf`
7. `results/deseq2/comparison_summary.csv`
8. individual comparison tables and volcano plots

Do not exclude a sample using PCA alone. Cross-check alignment, assignment, library size, correlation, and experimental notes.

## Reproducibility Outputs

The DESeq2 step writes:

```text
analysis_complete.flag
analysis_report.csv
analysis_summary.csv
filter_summary.csv
results_names.csv
sessionInfo.txt
software_versions.csv
run_completion.csv
```

Preserve these with the final results.

# Exporting and Archiving Results

The working Sapelo2 copies should be exported after successful completion.

## Create a Pipeline Code Package

```bash
mkdir -p final_pipeline_files/metadata
mkdir -p final_pipeline_files/reference_record
```

```bash
cp \
    config.sh \
    validate_metadata.R \
    00_project_setup.sh \
    01_fastp.slurm \
    02_star_index.slurm \
    03_star_align_array.slurm \
    04_featureCounts.slurm \
    05_multiqc.slurm \
    06_deseq2_analysis.R \
    07_deseq2.slurm \
    run_full_pipeline.sh \
    README.md \
    final_pipeline_files/
```

```bash
cp metadata/sample_key.csv \
    final_pipeline_files/metadata/
```

Record reference filenames:

```bash
find reference \
    -maxdepth 1 \
    -type f \
    \( -name "*.fa" -o \
       -name "*.fasta" -o \
       -name "*.fna" -o \
       -name "*.gff" -o \
       -name "*.gff3" -o \
       -name "*.gtf" \) \
    -print \
    > final_pipeline_files/reference_record/reference_files.txt
```

Record checksums:

```bash
md5sum reference/* \
    > final_pipeline_files/reference_record/reference_md5sums.txt
```

Record R reproducibility files:

```bash
cp \
    results/deseq2/sessionInfo.txt \
    results/deseq2/software_versions.csv \
    final_pipeline_files/
```

Create the code archive:

```bash
tar -czf Anid_RNA_final_pipeline.tar.gz \
    final_pipeline_files
```

## Create a Main Results Package

```bash
tar -czf Anid_RNA_main_results.tar.gz \
    results/deseq2 \
    results/multiqc \
    results/counts/gene_counts.txt \
    results/counts/gene_counts.txt.summary
```

## Create a Reproducibility Package

This archive includes logs and tool-specific QC but excludes large BAM files:

```bash
tar -czf Anid_RNA_reproducibility_package.tar.gz \
    final_pipeline_files \
    results/deseq2 \
    results/multiqc \
    results/counts \
    results/fastp/*.fastp.json \
    results/fastp/*.fastp.html \
    results/star/*/*.Log.final.out \
    metadata/sample_key.csv \
    logs
```

## Verify Archives

```bash
ls -lh Anid_RNA_*.tar.gz
```

```bash
tar -tzf Anid_RNA_final_pipeline.tar.gz | head -n 30
```

## Download Through the Transfer Node

Run these commands from a terminal on the local computer:

```bash
scp \
    YourMyID@xfer.gacrc.uga.edu:/full/path/to/Anid_RNA/Anid_RNA_final_pipeline.tar.gz \
    .
```

```bash
scp \
    YourMyID@xfer.gacrc.uga.edu:/full/path/to/Anid_RNA/Anid_RNA_main_results.tar.gz \
    .
```

```bash
scp \
    YourMyID@xfer.gacrc.uga.edu:/full/path/to/Anid_RNA/Anid_RNA_reproducibility_package.tar.gz \
    .
```

Find the project path on Sapelo2 with:

```bash
pwd
```

## Recommended Files to Preserve

At minimum, retain:

```text
Anid_RNA_final_pipeline.tar.gz
Anid_RNA_main_results.tar.gz
Anid_RNA_reproducibility_package.tar.gz
```

Also retain the original FASTQ files in appropriate long-term storage. BAM files are useful for browser visualization and re-counting, but are not required to reproduce the analysis if the original FASTQs, reference files, and working pipeline are preserved.

## Storage Notes

Use:

- `/home` for scripts and stable files
- `/scratch` for active temporary analyses
- `/work` for reusable project data
- `/project` for retained project data accessible through transfer nodes

Move important results out of scratch because scratch is not backed up and may be subject to purge policies.

## Troubleshooting Checklist

If a job disappears from `squeue`:

```bash
sacct -j JOBID \
    --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem
```

Then inspect:

```bash
cat logs/job_JOBID.err
```

```bash
tail -n 100 logs/job_JOBID.out
```

If a Bash script fails immediately:

```bash
bash -n script.slurm
```

If an R script fails immediately:

```bash
Rscript -e "invisible(parse(file='script.R')); cat('R syntax OK\\n')"
```

If R reports a missing package, test the package under the exact R module used by the Slurm job.

If code was copied from a formatted browser window, check for HTML corruption:

```bash
grep -nE '<strong|</strong>|&gt;|&lt;|&amp;|<br' script.R
```

## Companion Command Reference

See:

```text
Sapelo2_Command_Cheat_Sheet.md
```

for commonly used Sapelo2, Slurm, module, transfer, validation, and troubleshooting commands.
