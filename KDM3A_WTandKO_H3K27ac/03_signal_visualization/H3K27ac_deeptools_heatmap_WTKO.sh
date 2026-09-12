#!/bin/bash
#SBATCH --job-name=deept
#SBATCH --cpus-per-task=15
#SBATCH --time=12:00:00
#SBATCH --output=deept_%j.output
#SBATCH --error=deept_%j.errors

# deeptools computeMatrix + plotHeatmap
# WT/KO, CTRL/LPS H3K27ac signal, centered on peak summit
# replicate BigWigs are averaged per condition with bigwigAverage before plotting
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
WT_CTRL_REPLICATES=(
    "/path/to/H3K27ac_WT-trimmed.sorted.dups_marked.filtered.bw"
    "/path/to/H3K27ac_WT-trimmed.sorted.dups_marked.filtered.bw"
)
KO_CTRL_REPLICATES=(
    "/path/to/H3K27ac_KO-trimmed.sorted.dups_marked.filtered.bw"
    "/path/to/H3K27ac_KO-trimmed.sorted.dups_marked.filtered.bw"
)
WT_LPS_REPLICATES=(
    "/path/to/H3K27ac_LPSWT-trimmed.sorted.dups_marked.filtered.bw"
    "/path/to/H3K27ac_LPSWT-trimmed.sorted.dups_marked.filtered.bw"
)
KO_LPS_REPLICATES=(
    "/path/to/H3K27ac_LPSKO-trimmed.sorted.dups_marked.filtered.bw"
    "/path/to/H3K27ac_LPSKO-trimmed.sorted.dups_marked.filtered.bw"
)
LABELS=(WT_CTRL KO_CTRL WT_LPS KO_LPS)

BEDFILES=(
    "/path/to/KDM3A_CTRL_specific.bed"
    "/path/to/KDM3A_LPS_specific.bed"
    "/path/to/intersect_KDM3A_CTRL_LPS.bed"
)

OUTDIR="./deeptools_output"
OUTNAME="H3K27ac_WTKO_heatmap"
CONDA_ENV_DEEPTOOLS="/path/to/conda_envs/deeptools_env"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_DEEPTOOLS"

# average replicate bigwigs per condition
WT_CTRL_AVG="${OUTDIR}/WT_CTRL_averaged.bw"
KO_CTRL_AVG="${OUTDIR}/KO_CTRL_averaged.bw"
WT_LPS_AVG="${OUTDIR}/WT_LPS_averaged.bw"
KO_LPS_AVG="${OUTDIR}/KO_LPS_averaged.bw"

bigwigAverage -b "${WT_CTRL_REPLICATES[@]}" -o "$WT_CTRL_AVG"
bigwigAverage -b "${KO_CTRL_REPLICATES[@]}" -o "$KO_CTRL_AVG"
bigwigAverage -b "${WT_LPS_REPLICATES[@]}"  -o "$WT_LPS_AVG"
bigwigAverage -b "${KO_LPS_REPLICATES[@]}"  -o "$KO_LPS_AVG"

BIGWIGS=("$WT_CTRL_AVG" "$KO_CTRL_AVG" "$WT_LPS_AVG" "$KO_LPS_AVG")

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
