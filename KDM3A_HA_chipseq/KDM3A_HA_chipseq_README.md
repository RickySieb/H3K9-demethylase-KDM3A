# KDM3A-HA ChIP-seq

KDM3A-HA ChIP-seq in CTRL and LPS-stimulated conditions. Peaks are compared between
conditions to identify CTRL-specific, LPS-specific, and constitutive (shared) binding sites.

## Pipeline order

1. **`01_hpc_mapping/`** — Trimming, alignment, and BAM filtering via the omnomnomics
   pipeline. See top-level README for pipeline version/commit. Exact command used is in
   this folder.
2. **`02_hpc_peakcalling_IDR/`** — MACS3 peak calling on replicate BAMs, followed by IDR
   to identify reproducible peaks per condition (CTRL, LPS).
   - `chipseq_macs3_idr.sh`
   - `chipseq_bedtools_peaks.sh`
3. **`03_HOMER_peak_annotation/`** — HOMER genomic annotation and GO enrichment for each
   of the three peak sets above. Genome: mm10.
   - `homer_annotate_peaks.sh`
4. **`04_signal_visualization/`** — deeptools heatmaps and karyoploteR trackshots of
   KDM3A signal over the three peak sets / at specific loci.
   - `chipseq_deeptools_heatmap.sh` — merges CTRL and LPS replicate BigWigs, generates
     heatmap over the three peak sets
   - `KDM3AHA_signal_ifih1.R` — trackshot at the ifih1 locus

## Figure map

| Figure | Script | Notes |
|---|---|---|
| Fig 2d (heatmap) | `04_signal_visualization/chipseq_deeptools_heatmap.sh` | CTRL vs LPS, 3 peak sets |
| Fig 2g (trackshot) | `04_signal_visualization/KDM3AHA_signal_ifih1.R` | ifih1 locus |
| Extended data Fig 2b (annotation bar plot) | `03_HOMER_peak_annotation/KDM3AHA_peakannotationR.R` | genomic feature distribution |

## Notes

- Peak calling used MACS3 (p-value cutoff: see script) and IDR (threshold 0.05).
- Raw BAMs and BigWigs are not included here; see top-level README for accession.
