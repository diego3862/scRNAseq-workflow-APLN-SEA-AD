import gc
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from src.data_utils import *
from src.plot_utils import *

#Definitions
subtype_custom = "APJ_Cats"
chosen_subtypes = ["Astrocyte", "Neuronal"]  
genes_of_interest = ["APLN", "APLNR"]
lfc_threshold = 0.25
FIG_DIR = "figures"
os.makedirs("data", exist_ok=True)
set_paper_theme()

#Load Data & Assign Categories Once
file_choice = "SEAAD_HIP_RNAseq_final-nuclei.2026-06-22.h5ad"
print(f"Loading {file_choice}...")
adata = sc.read_h5ad(f"data/{file_choice}")
assign_APJ_Cats(adata)

gene_dict = {}

#Iterate over all requested cell subtypes
for subtype in chosen_subtypes:
  print(f"\n{'='*25} Processing: {subtype} {'='*25}")

  #Subset to current cell type
  sub_adata = adata[adata.obs[subtype_custom] == subtype].copy()
  print(f"Total cells: {sub_adata.n_obs}")

  #Assign binary positivity for each gene
  for gene in genes_of_interest:
    status_col = f"{gene}_status"
    gene_expr = sub_adata[:, gene].X

    if sp.issparse(gene_expr):
      gene_expr = gene_expr.toarray().flatten()
    elif hasattr(gene_expr, "get"):  #Handle CuPy/GPU arrays if using rapids
      gene_expr = gene_expr.get().flatten()
    else:
      gene_expr = np.asarray(gene_expr).flatten()

    sub_adata.obs[status_col] = np.where(gene_expr > 0, f"{gene}+", f"{gene}-")
    sub_adata.obs[status_col] = sub_adata.obs[status_col].astype("category")

    print(f"\n{gene} status breakdown in {subtype}:")
    print(sub_adata.obs[status_col].value_counts())

  #Run Wilcoxon DE, Save Unfiltered, and Filter DEGs
  for gene in genes_of_interest:
    status_col = f"{gene}_status"
    ref_level = f"{gene}-"

    print(f"\nRunning Wilcoxon for {subtype} - {gene}...")
    df_res = run_wilcoxon(
        sub_adata, design_factor=status_col, reference_level=ref_level
    )

    # Save unfiltered DEG table with dynamic name
    df_res.to_csv(f"data/{gene}_{subtype}_DEGs_unfiltered.csv")

    # Filter by significance and effect size
    res_clean = df_res.dropna(subset=["padj"])
    sig_degs = res_clean[res_clean["padj"] < 0.05]

    up_genes = sig_degs[
        sig_degs["log2FoldChange"] > lfc_threshold
    ].index.tolist()
    down_genes = sig_degs[
        sig_degs["log2FoldChange"] < -lfc_threshold
    ].index.tolist()

    print(
        f"{subtype} - {gene} DEGs (|LFC| > {lfc_threshold}):"
        f" {len(up_genes) + len(down_genes)}"
    )
    print(
        f"  -> Upregulated: {len(up_genes)} | Downregulated: {len(down_genes)}"
    )

    # Store with combined key for multi-panel GSEA (e.g. ('Astrocyte_APLN', 'Up'))
    dict_key = f"{subtype}_{gene}"
    gene_dict[(dict_key, "Up")] = up_genes
    gene_dict[(dict_key, "Down")] = down_genes

  # Clean up memory before next subtype
  del sub_adata
  gc.collect()

#Run GSEApy once across all collected conditions
print(f"\n{'='*25} Running GSEApy {'='*25}")
libraries = ["KEGG_2021_Human"]
gsea_df = get_signed_enrichment_long(gene_dict, libraries,top_n=50)

#Save outputs
if gsea_df is not None:
  # Save master combined dataframe
  gsea_df.to_csv("data/Combined_GSEA_results.csv", index=False)
  print("Saved master dataframe to 'data/Combined_GSEA_results.csv'")

  # Save individual files per subtype and gene
  for subtype in chosen_subtypes:
    for gene in genes_of_interest:
      key = f"{subtype}_{gene}"
      sub_gsea = gsea_df[gsea_df["Subtype"] == key].copy()
      if not sub_gsea.empty:
        sub_gsea.to_csv(f"data/{gene}_{subtype}_GSEA.csv", index=False)
        print(f"Saved: data/{gene}_{subtype}_GSEA.csv")
      elif sub_gsea.empty:
        print(f"{gene}, {subtype} gsea is empty.")
else:
  print("No significant pathways found across any condition.")