# SEA-AD Microglia scRNA‑seq dataset Immunometabolism Analysis Workflow
A reproducible workflow for pseudobulk analysis of the publicly available Seattle Alzheimer's Disease Brian Cell Atlas (SEA-AD) single-nucleus RNA-seq datasets. The repository includes data processing scripts and figure generation notebooks.
This analysis explores the SEA-AD HIP, DFC, and MTG datasets from the perspective of APLN and APLNR expression, aiming to characterise expression changes and their consequences across neurons, astrocytes, and microglia.
## Repository structure:
### workflow
* 01_pseudobulk_overall.py: takes HIP, MTG and DFC truncated pseudobulk objects and processes them into pseudobulk objects.
* 02_HIP_A_N_DE_GSEA.py: takes HIP h5ad object and performs cell-level differential expression (wilcoxon) and over-representation analysis.
* 03_table_1_AD_confounders.py: Describes confounding donor-level variables for the pseudobulk AD analysis
* 04_table_2_cell_confounders.py: Describes confounding cell-level variables for the differential expression.
### notebooks
* figures.ipynb: plot figures from the data objects generated in the scripts contained in workflow/
* h5ad_slicing.ipynb: load h5ad objects in small chunks and keep a subset of genes to allow processing of objects that would be too large (>30GB) otherwise. 
### src
* data_utils.py - Contains data processing helper functions for various files in this repository. 
* plot_utils.py - Contains plotting helper functions for various files in this repository.

## References and Data Availability
The SEA-AD consortium aims to build a map of brain cell types in aging and AD. It achieves this by generating publicly available multi-omic resources such as that used in this repository. 
* **Original Paper:** Travaglini et al. (2026) -[*Multiregional single-cell profiling reveals shared and specialized cellular vulnerability in Alzheimer’s disease*](https://www.biorxiv.org/content/10.64898/2026.07.01.734821v1)
* **Data available:** The full scRNAseq and scATACseq dataset can be found here: https://brain-map.org/consortia/sea-ad/our-data. The specific .h5ad object used for this workflow ("SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC-pre-release.2025-07-24.h5ad"), is now available at: https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#Multiregion_2026/previous_objects/

## Installation
* Running in google Colab (recommended). Open the notebook:notebooks/run_workflow.ipynb. Follow google drive comment instructions.
* Running locally (pip install -r requirements.txt)
