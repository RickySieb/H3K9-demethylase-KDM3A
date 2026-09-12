#!/bin/bash
#SBATCH --job-name=HOMER
#SBATCH --cpus-per-task=63
#SBATCH --time=01:00:00
#SBATCH --output=HOMER_%j.output
#SBATCH --error=HOMER_%j.errors

# HOMER peak annotation + GO enrichment for the three KDM3A peak sets
# (CTRL-specific, LPS-specific, constitutive), mm10
# R. Siebeler

set -euo pipefail

# --- settings, edit before running ---
GENOME="mm10"
HOMER_BIN="$HOME/bin/homer/bin/annotatePeaks.pl"

BEDFILES=(
    "/path/to/bedtools/KDM3A_CTRL_specific.bed"
    "/path/to/bedtools/KDM3A_LPS_specific.bed"
    "/path/to/bedtools/intersect_KDM3A_CTRL_LPS.bed"
)

OUTDIR="./Peak_annotation"
mkdir -p "$OUTDIR"

for BED in "${BEDFILES[@]}"; do
    NAME=$(basename "$BED" .bed)
    echo "annotating $NAME"

    GO_DIR="$OUTDIR/${NAME}_GO"
    mkdir -p "$GO_DIR"

    "$HOMER_BIN" "$BED" "$GENOME" -go "$GO_DIR" > "$OUTDIR/${NAME}_annotated.txt"
done

echo "done, output in $OUTDIR"
