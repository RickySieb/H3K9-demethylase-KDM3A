#!/usr/bin/env bash
# bedtools intersect - common peaks + peaks unique to each condition
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
BED_A="/path/to/conditionA_IDR_peaks.bed"
BED_B="/path/to/conditionB_IDR_peaks.bed"
NAMEA="conditionA"
NAMEB="conditionB"
OUTDIR="./bedtools_output"

mkdir -p "$OUTDIR"

module load bedtools

# peaks present in both
bedtools intersect -a "$BED_A" -b "$BED_B" > "$OUTDIR/${NAMEA}_${NAMEB}_common.bed"

# peaks only in A
bedtools intersect -a "$BED_A" -b "$BED_B" -v > "$OUTDIR/${NAMEA}_unique.bed"

# peaks only in B
bedtools intersect -a "$BED_B" -b "$BED_A" -v > "$OUTDIR/${NAMEB}_unique.bed"

echo "done, output in $OUTDIR"
