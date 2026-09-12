#!/bin/bash
#SBATCH --job-name=deept
#SBATCH --cpus-per-task=15
#SBATCH --time=12:00:00
#SBATCH --output=deept_%j.output
#SBATCH --error=deept_%j.errors

# merge RNA Pol II ChIP-seq replicate BigWigs per condition using deeptools bigwigAverage
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
WT_CTRL_REPLICATES=(
    "/path/to/RNApolII_WT_rep1.bw"
    "/path/to/RNApolII_WT_rep2.bw"
)
KO_CTRL_REPLICATES=(
    "/path/to/RNApolII_KO_rep1.bw"
    "/path/to/RNApolII_KO_rep2.bw"
)
WT_LPS_REPLICATES=(
    "/path/to/RNApolII_WTLPS_rep1.bw"
    "/path/to/RNApolII_WTLPS_rep2.bw"
)
KO_LPS_REPLICATES=(
    "/path/to/RNApolII_KOLPS_rep1.bw"
    "/path/to/RNApolII_KOLPS_rep2.bw"
)

OUTDIR="./avgBigWig"
CONDA_ENV_DEEPTOOLS="/path/to/conda_envs/deeptools_env"

mkdir -p "$OUTDIR"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_DEEPTOOLS"

bigwigAverage -b "${WT_CTRL_REPLICATES[@]}" -o "$OUTDIR/RNApolII_WT_averaged.bw"
bigwigAverage -b "${KO_CTRL_REPLICATES[@]}" -o "$OUTDIR/RNApolII_KO_averaged.bw"
bigwigAverage -b "${WT_LPS_REPLICATES[@]}"  -o "$OUTDIR/RNApolII_WT_LPS_averaged.bw"
bigwigAverage -b "${KO_LPS_REPLICATES[@]}"  -o "$OUTDIR/RNApolII_KO_LPS_averaged.bw"

conda deactivate

echo "done, averaged bigwigs in $OUTDIR"
