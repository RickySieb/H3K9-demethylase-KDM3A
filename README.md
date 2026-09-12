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

All public datasets utilized in the current study are available in the Gene Expression Omnibus:  
GSE70619. RNA-, ChIP-, and ATAC-seq data are available through Zenodo (DOIs:  10.5281/zenodo.19566015, 10.5281/zenodo.19567228). RNA-seq data and associated clinical metadata  
from part of the Athero-Express (AE) cohort are available via DataverseNL (DOIs:  https://doi.org/10.34894/D1MDKL, https://doi.org/10.34894/TYHGEF,  
https://doi.org/10.34894/4IKE3T ). Due to data governance and privacy regulations, access to private  patient data is controlled and can be requested through DataverseNL.

## External pipelines used

Mapping and initial processing used the omnomnomics pipeline (https://github.com/prangelab/omnomnomics_BASH, v0.4. Exact commands and module versions used are recorded in each experiment's `01_hpc_mapping/` folder.
