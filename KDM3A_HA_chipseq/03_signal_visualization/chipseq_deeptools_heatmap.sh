#!/bin/bash
#SBATCH --job-name=deept
#SBATCH --cpus-per-task=15
#SBATCH --time=12:00:00
#SBATCH --output=deept_%j.output
#SBATCH --error=deept_%j.errors

# deeptools computeMatrix + plotHeatmap for KDM3A ChIP-seq
# CTRL/LPS-specific and shared peaks, centered on peak summit
# replicate BigWigs are averaged with bigwigAverage before plotting
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
CTRL_REPLICATES=(
    "/path/to/CTRL_rep1.bw"
    "/path/to/CTRL_rep2.bw"
)
LPS_REPLICATES=(
    "/path/to/LPS_rep1.bw"
    "/path/to/LPS_rep2.bw"
)
LABELS=(CTRL LPS)

BEDFILES=(
    "/path/to/KDM3A_CTRL_specific.bed"
    "/path/to/KDM3A_LPS_specific.bed"
    "/path/to/intersect_KDM3A_CTRL_LPS.bed"
)

OUTDIR="./deeptools_output"
OUTNAME="KDM3A_heatmap"
CONDA_ENV_DEEPTOOLS="/path/to/conda_envs/deeptools_env"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_DEEPTOOLS"

# average replicate bigwigs per condition
CTRL_AVG="${OUTDIR}/CTRL_averaged.bw"
LPS_AVG="${OUTDIR}/LPS_averaged.bw"

bigwigAverage -b "${CTRL_REPLICATES[@]}" -o "$CTRL_AVG"
bigwigAverage -b "${LPS_REPLICATES[@]}" -o "$LPS_AVG"

BIGWIGS=("$CTRL_AVG" "$LPS_AVG")

MATRIX_OUT="${OUTNAME}_matrix.gz"

computeMatrix reference-point \
    --referencePoint center \
    -S "${BIGWIGS[@]}" \
    -R "${BEDFILES[@]}" \
    -a 5000 -b 5000 \
    --sortRegions descend \
    -o "$MATRIX_OUT"

plotHeatmap -m "$MATRIX_OUT" \
    -out "${OUTNAME}_heatmap.png" \
    --whatToShow "heatmap and colorbar" \
    --samplesLabel "${LABELS[@]}"

conda deactivate

echo "done: ${OUTDIR}/${OUTNAME}_heatmap.png"
