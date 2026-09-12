# H3K9-demethylase-KDM3A

Code accompanying [paper title], [journal], [year]. DOI: [add once available]

## Repository structure

- `KDM3A_HA_chipseq/` — KDM3A-HA ChIP-seq: mapping, peak calling, peak annotation, signal visualization
- `KDM3A_WTandKO_ATACseq/` — ATAC-seq in WT and KO, ± LPS: mapping, differential accessibility, signal visualization
- `KDM3A_WTandKO_H3K27ac/` — H3K27ac ChIP-seq in WT and KO, ± LPS: mapping, signal analysis, signal visualization
- `KDM3A_RNApolII_chipseq/` — RNA Pol II ChIP-seq: mapping, signal visualization
- `KDM3A_scRNAseq_humanplaques/` — KDM3A expression across leukocyte subtypes in human plaque scRNA-seq data

Each folder contains its own README describing the pipeline order and inputs/outputs for that analysis.

## Environment and dependencies

- R package versions: see `environment/session_info.txt`
- Conda environments: see `environment/*.yml`

## Raw data

Raw sequencing data are deposited at GEO under accession [add accession].

## External pipelines used

Mapping and initial processing used the omnomnomics pipeline (https://github.com/prangelab/omnomnomics_BASH, v0.4. Exact commands and module versions used are recorded in each experiment's `01_hpc_mapping/` folder.
