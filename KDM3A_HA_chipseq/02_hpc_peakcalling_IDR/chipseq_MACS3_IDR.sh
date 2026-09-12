#!/usr/bin/env bash
# MACS3 peak calling + IDR for two ChIP-seq replicates
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
BAMA="/path/to/replicate_A_filtered.bam"
BAMB="/path/to/replicate_B_filtered.bam"
PVAL=0.01
OUTDIR="./IDR_output/my_sample"
IDR_THRESHOLD=0.05

CONDA_ENV_MACS3="/path/to/conda_envs/macs3_env"
CONDA_ENV_IDR="/path/to/conda_envs/idr_env"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

BASENAMEA=$(basename "$BAMA" .bam)
BASENAMEB=$(basename "$BAMB" .bam)

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# peak calling
conda activate "$CONDA_ENV_MACS3"

macs3 callpeak -t "$BAMA" --outdir . -n "${BASENAMEA}_narrow.MACS3.p${PVAL}" -p "$PVAL" --verbose 0
macs3 callpeak -t "$BAMB" --outdir . -n "${BASENAMEB}_narrow.MACS3.p${PVAL}" -p "$PVAL" --verbose 0

conda deactivate

# sort peaks, needed for IDR
sort -k8,8nr "${BASENAMEA}_narrow.MACS3.p${PVAL}_peaks.narrowPeak" > "${BASENAMEA}_sortedPeaks.narrowPeak"
sort -k8,8nr "${BASENAMEB}_narrow.MACS3.p${PVAL}_peaks.narrowPeak" > "${BASENAMEB}_sortedPeaks.narrowPeak"

# IDR between replicates
conda activate "$CONDA_ENV_IDR"

idr --samples "${BASENAMEA}_sortedPeaks.narrowPeak" "${BASENAMEB}_sortedPeaks.narrowPeak" \
    --input-file-type narrowPeak \
    --output-file "${BASENAMEA}_${BASENAMEB}_idr.txt" \
    --plot \
    --log-output-file "${BASENAMEA}_${BASENAMEB}_idr.log"

conda deactivate

# col 5 is -log10(local IDR)*100, so 540 = IDR 0.05
awk '($5 >= 540)' "${BASENAMEA}_${BASENAMEB}_idr.txt" > "${BASENAMEA}_${BASENAMEB}_IDR${IDR_THRESHOLD}Peaks.bed"

echo "done, reproducible peaks: ${BASENAMEA}_${BASENAMEB}_IDR${IDR_THRESHOLD}Peaks.bed"
