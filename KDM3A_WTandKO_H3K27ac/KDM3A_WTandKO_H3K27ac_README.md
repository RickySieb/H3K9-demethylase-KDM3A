# H3K27ac ChIP-seq: WT vs KO, ± LPS

H3K27ac ChIP-seq in WT and KDM3A-KO cells, ± LPS stimulation. Signal is quantified
near genes downregulated upon KDM3A loss, to assess whether reduced H3K27ac
(active enhancer/promoter mark) accompanies their downregulation.

## Pipeline order

1. **`01_hpc_mapping/`** — Trimming, alignment, and BAM filtering via the omnomnomics
   pipeline (see top-level README for version/commit). Exact command used is in this
   folder.
2. **`02_signal_analysis/`** — Quantify H3K27ac signal near downregulated genes
   across all 8 samples (WT/KO, ±LPS, 2 replicates each), then compare groups
   statistically.
   - `H3K27ac_analysis.R` — coverage extraction per region, box plots (log2 signal),
     Wilcoxon test with BH correction
3. **`03_signal_visualization/`** — *(add once available)* deeptools heatmap and/or
   karyoploteR trackshot of H3K27ac signal, following the same pattern as the ATAC
   and KDM3A ChIP-seq visualization scripts.

## Figure map

| Figure | Script | Notes |
|---|---|---|
| Fig Xa (box plots + stats) | `02_signal_analysis/H3K27ac_analysis.R` | H3K27ac near downregulated genes, WT vs KO |

## Known issue to resolve before upload

`H3K27ac_analysis.R` currently references `df1`/`df2` (the BED/region set defining
"genes downregulated in KO" and "genes downregulated in KO+LPS") that are not loaded
anywhere in the script itself — see the `TODO` comment in the script. Add the actual
loading code (or the BED file paths) for these two gene sets before this script can
run standalone. Until resolved, this script will only work if `df1`/`df2` already
exist in the R session from a prior analysis step.

## Notes

- Raw BAMs and BigWigs are not included; see top-level README for GEO accession.
