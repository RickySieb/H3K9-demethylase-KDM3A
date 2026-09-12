# ATAC-seq: WT vs KO, ± LPS

ATAC-seq in WT and KDM3A-KO cells, ± LPS stimulation. Signal is quantified at
KDM3A ChIP-seq peaks (CTRL-specific, LPS-specific, constitutive) to assess whether
KDM3A loss affects chromatin accessibility at its own binding sites.

## Pipeline order

1. **`01_hpc_mapping/`** — Trimming, alignment, and BAM filtering via the omnomnomics
   pipeline (see top-level README for version/commit). Exact command used is in this
   folder.
2. **`02_signal_analysis/`** — Quantify ATAC signal at the three KDM3A peak sets
   (from `KDM3A_HA_chipseq/03_bedtools_peaks/`) across all 8 samples (WT/KO, ±LPS,
   2 replicates each), then compare groups statistically.
   - `ATAC_analysis.R` — coverage extraction per peak, box plots (log2 signal),
     Wilcoxon test with BH correction per peak set
3. **`03_signal_visualization/`** — deeptools heatmap and karyoploteR trackshot of
   ATAC signal, replicate BigWigs merged per condition first.
   - `atac_deeptools_heatmap_WTKO.sh` — merges 8 replicate BigWigs into 4 condition
     averages (WT_CTRL, KO_CTRL, WT_LPS, KO_LPS), heatmap over the 3 KDM3A peak sets
   - `ATAC_signal_ifih1.R` — trackshot at the IFIH1 locus, WT/KO x CTRL/LPS

## Figure map

| Figure | Script | Notes |
|---|---|---|
| Fig 2e (box plots + stats) | `02_signal_analysis/ATAC_analysis.R` | ATAC signal at KDM3A peaks, WT vs KO |
| Extended data Fig 2c (heatmap) | `03_signal_visualization/atac_deeptools_heatmap_WTKO.sh` | 4 conditions, 3 peak sets |
| Fig 2g (trackshot) | `03_signal_visualization/ATAC_signal_ifih1.R` | ifih1 locus |

## Notes

- Peak sets used here are generated in `KDM3A_HA_chipseq/`, not re-derived here.
- Raw BAMs and BigWigs are not included; see top-level README for GEO accession.
