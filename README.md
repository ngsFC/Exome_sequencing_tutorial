# Exome Sequencing Tutorial — From Reads to Variants, and Somatic Mutations in Wilms Tumor

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ngsFC/Exome_sequencing_tutorial/HEAD?urlpath=rstudio)

> An easy, hands-on introduction to whole-exome sequencing (WES): first the real command-line pipeline from raw reads to an annotated, benchmarked VCF; then the analysis of somatic mutations in a real pediatric cancer cohort (TARGET Wilms tumor).

**Author:** Francesco Cecere (Senior Bioinformatic Researcher) · AORN Santobono-Pausilipon IRCCS

## Modules

| File | Content |
|------|---------|
| `01_from_reads_to_variants.Rmd` | FASTQ → `fastp` → `bwa mem` → `samtools` (sort, duplicates, coverage) → `bcftools` (calling, filtering, normalization, `csq` annotation) → benchmark against the GIAB truth set |
| `02_somatic_mutations_wilms.Rmd` | Somatic MAF with `maftools`: one sample per patient, cohort summary, tumor mutational burden vs TCGA, oncoplots, TP53 lollipop, Favorable vs Unfavorable histology |

Module 1 runs real shell commands (`bash` chunks) inside the Binder container, so it works the same from any computer (Windows, macOS, Linux): only a web browser is needed. Its outputs go to `work/` (not committed). Pre-rendered reports are in `html/`.

## Data

**Module 1 — `data/germline/`** (GRCh37, chromosome 11:1–3,000,000, the 11p15.5 region)
- `NA12878_exome_R1/R2.fastq.gz` — paired-end reads of GIAB sample NA12878 (HG001), Illumina TruSeq exome, from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) Nebraska exome BAM.
- `reference_chr11p15.fa` — reference sequence of the region (UCSC hg19 = GRCh37).
- `exome_targets.bed` — TruSeq exome capture targets.
- `genes_Ensembl87.gff3` — Ensembl GRCh37 release 87 gene annotation, for `bcftools csq`.
- `giab_truth.vcf.gz`, `giab_confident.bed` — GIAB HG001 v4.2.1 benchmark variants and high-confidence regions.

**Module 2 — `data/somatic/`**
- `TARGET_WT_somatic_mutations.maf.gz` — all open-access WXS masked somatic mutations of TARGET-WT from the [GDC Data Portal](https://portal.gdc.cancer.gov/projects/TARGET-WT) (42 tumor samples, GRCh38).
- `TARGET_WT_clinical.tsv` — histology, age, sex and vital status for the same samples.

All built once by `data-prep/build_germline_data.sh` and `data-prep/build_somatic_data.R` (not run by the tutorial — included for audit/rebuild).

## Environment

Built by [repo2docker](https://repo2docker.readthedocs.io/) from `.binder/Dockerfile` (`FROM rocker/binder:4.4.1`), adding the Ubuntu packages `bwa`, `samtools`, `bcftools`, `tabix`, `fastp`, `bedtools`, and the R packages `tidyverse`, `jsonlite`, `DT`, `knitr`, `rmarkdown`, `R.utils` (needed by `maftools::tcgaCompare()` to read its bundled TCGA files) and Bioconductor's `maftools`. No network access is needed at runtime.
