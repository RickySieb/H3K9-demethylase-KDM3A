# scRNA-seq: KDM3A expression in human plaques

KDM3A expression across leukocyte subtypes in a human atherosclerotic plaque
scRNA-seq dataset, using a pre-integrated Seurat object.

## Pipeline order

1. **`01_KDM3A_expression/`** — Map fine-grained cell identities to main leukocyte
   subtypes, then visualize and summarize KDM3A expression across those subtypes.
   - `scAnalysis.R` — UMAP colored by subtype, KDM3A feature plot, KDM3A dot plot
     across subtypes, summary statistics table

## Figure map

| Figure | Script | Notes |
|---|---|---|
| Fig Xa (UMAP + feature plot + dot plot) | `01_KDM3A_expression/scAnalysis.R` | combined 3-panel figure |

## Notes

- This script expects a pre-processed, pre-integrated Seurat object as input
  (`SEURAT_OBJECT_PATH` at the top of the script) — the integration/clustering
  pipeline itself is not part of this folder.
- Cell identity-to-subtype mapping is manually defined in the script
  (`subtype_map`); update this if the source object's cluster labels change.
