# RNA Pol II ChIP-seq

RNA Pol II ChIP-seq used to assess transcriptional activity at KDM3A target loci.

## Pipeline order

1. **`01_hpc_mapping/`** — Trimming, alignment, and BAM filtering via the omnomnomics
   pipeline (see top-level README for version/commit). Exact command used is in this
   folder.
2. **`02_signal_visualization/`** — Merge replicate BigWigs per condition, then
   visualize signal.
   - `rnapolii_merge_bigwigs.sh` — merges replicate BigWigs into 4 condition
     averages (WT, KO, WT+LPS, KO+LPS) using deeptools bigwigAverage
   - *(add once available)* trackshot or heatmap script using the merged BigWigs
     above

## Figure map

| Figure | Script | Notes |
|---|---|---|
| Fig Xa | *(add once the visualization script for this experiment is finalized)* | |

## Notes

- Merged BigWigs from `02_signal_visualization/rnapolii_merge_bigwigs.sh` are the
  input for whatever trackshot/heatmap script is added here.
- Raw BAMs and BigWigs are not included; see top-level README for GEO accession.
