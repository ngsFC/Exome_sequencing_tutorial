#!/usr/bin/env bash
# Builds every file under data/germline/ (Module 1) from public sources:
#   - reads:     GIAB NA12878 (HG001) TruSeq exome, Nebraska, aligned BAM (GRCh37)
#   - reference: hg19 chr11:1-3,000,000 from the UCSC API (identical to GRCh37 "11")
#   - targets:   TruSeq exome capture regions from the same GIAB folder
#   - genes:     Ensembl GRCh37 release 87 GFF3, chromosome 11
#   - truth:     GIAB HG001 GRCh37 v4.2.1 benchmark VCF + confident regions
#
# This script is NOT run during the tutorial -- it was run once to produce the
# files committed under data/germline/, so that Binder sessions need no large
# download. Only the reads overlapping the region are fetched (samtools reads
# the remote BAM through its index).
#
# Requires: samtools, bcftools, tabix, curl, python3.
set -euo pipefail

REGION_CHR=11
REGION_END=3000000
OUT=data/germline
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
BAM=$GIAB/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam
TRUTH=$GIAB/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"

# ---- 1. Reads: every read pair aligned to the region, turned back into FASTQ ----
# The BAM is the only public form of these reads small enough to slice; turning
# it back into FASTQ lets the tutorial start from raw reads, as a real project does.
curl -sS --retry 5 --retry-all-errors -o "$TMP/region.bam.bai" "$BAM.bai"
samtools view -b -X -o "$TMP/region.bam" "$BAM" "$TMP/region.bam.bai" "$REGION_CHR:1-$REGION_END"
samtools collate -u -O "$TMP/region.bam" |
  samtools fastq -n -1 "$OUT/NA12878_exome_R1.fastq.gz" -2 "$OUT/NA12878_exome_R2.fastq.gz" -0 /dev/null -s /dev/null -

# ---- 2. Reference: the same 3 Mb, as a FASTA with a single contig "11" ----
curl -sS --retry 5 --retry-all-errors "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr$REGION_CHR;start=0;end=$REGION_END" -o "$TMP/seq.json"
python3 -I -c '
import json, sys, textwrap
seq = json.load(open(sys.argv[1]))["dna"].upper()
open(sys.argv[2], "w").write(">" + sys.argv[3] + "\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n")
' "$TMP/seq.json" "$OUT/reference_chr11p15.fa" "$REGION_CHR"

# ---- 3. Exome capture targets inside the region ----
curl -sS --retry 5 --retry-all-errors "$GIAB/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed" |
  awk -v c="$REGION_CHR" -v e="$REGION_END" 'BEGIN{OFS="\t"} {sub(/^chr/, "", $1)} $1 == c && $3 <= e {print $1, $2, $3}' \
  > "$OUT/exome_targets.bed"

# ---- 4. Gene annotation for bcftools csq: genes entirely inside the region ----
curl -sS --retry 5 --retry-all-errors -o "$TMP/genes.gff3.gz" \
  https://ftp.ebi.ac.uk/ensemblorg/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chromosome.$REGION_CHR.gff3.gz
python3 -I data-prep/subset_gff.py "$TMP/genes.gff3.gz" "$OUT/genes_Ensembl87.gff3" "$REGION_END"

# ---- 5. GIAB truth set: variants and high-confidence regions inside the region ----
curl -sS --retry 5 --retry-all-errors -o "$TMP/truth.vcf.gz.tbi" "$TRUTH.vcf.gz.tbi"
bcftools view -r "$REGION_CHR:1-$REGION_END" -Oz -o "$OUT/giab_truth.vcf.gz" "$TRUTH.vcf.gz##idx##$TMP/truth.vcf.gz.tbi"
tabix -p vcf "$OUT/giab_truth.vcf.gz"
curl -sS --retry 5 --retry-all-errors "$TRUTH.bed" | awk -v c="$REGION_CHR" -v e="$REGION_END" 'BEGIN{OFS="\t"} $1 == c && $2 < e {if ($3 > e) $3 = e; print $1, $2, $3}' \
  > "$OUT/giab_confident.bed"

ls -la "$OUT"
