import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from src.data_utils import assign_APJ_Cats

#Load hippocampus data for APLN APLNR +/- comparisons.
file_choice = "SEAAD_HIP_RNAseq_final-nuclei.2026-06-22.h5ad"
print(f"Loading {file_choice}...")
adata_HIP = sc.read_h5ad(f"data/{file_choice}")
assign_APJ_Cats(adata_HIP)

# TABLE 2: Cell-Level Confounders (HIP dataset: Neurons & Astrocytes)
# Ensure binary positivity exists in adata_HIP
for gene in ["APLN", "APLNR"]:
  gene_expr = adata_HIP[:, gene].X
  gene_expr = (
      gene_expr.toarray().flatten()
      if sp.issparse(gene_expr)
      else np.asarray(gene_expr).flatten()
  )
  adata_HIP.obs[f"{gene}_status"] = np.where(
      gene_expr > 0, f"{gene}+", f"{gene}-"
  )

records_t2 = []
for cell_type in ["Neuronal", "Astrocyte"]:
  sub_ct = adata_HIP[adata_HIP.obs["APJ_Cats"] == cell_type]

  for gene in ["APLN", "APLNR"]:
    status_col = f"{gene}_status"

    for status in [f"{gene}+", f"{gene}-"]:
      sub = sub_ct[sub_ct.obs[status_col] == status]
      n_cells = len(sub)
      if n_cells == 0:
        continue

      n_donors = sub.obs["Donor ID"].nunique()
      umi_str = (
          f"{sub.obs['Number of UMIs'].mean():.1f} ±"
          f" {sub.obs['Number of UMIs'].std():.1f}"
          if "Number of UMIs" in sub.obs
          else "N/A"
      )
      genes_str = (
          f"{sub.obs['Genes detected'].mean():.1f} ±"
          f" {sub.obs['Genes detected'].std():.1f}"
          if "Genes detected" in sub.obs
          else "N/A"
      )

      if "Fraction mitochondrial UMIs" in sub.obs:
        mito_val = sub.obs["Fraction mitochondrial UMIs"].mean() * 100
        mito_sd = sub.obs["Fraction mitochondrial UMIs"].std() * 100
        mito_str = f"{mito_val:.2f}% ± {mito_sd:.2f}%"
      else:
        mito_str = "N/A"

      doublet_str = (
          f"{sub.obs['Doublet score'].mean():.3f} ±"
          f" {sub.obs['Doublet score'].std():.3f}"
          if "Doublet score" in sub.obs
          else "N/A"
      )

      records_t2.append({
          "Cell Lineage": cell_type,
          "Gene Target": gene,
          "Status": status,
          "Total Cells": n_cells,
          "Contributing Donors": n_donors,
          "Total UMIs / Cell": umi_str,
          "Genes Detected / Cell": genes_str,
          "Mitochondrial Reads (%)": mito_str,
          "Doublet Score": doublet_str,
      })

table2_df = pd.DataFrame(records_t2)
table2_df.to_csv("data/Supplementary_Table_2:_Cell_Quality_Confounders.csv", index=False)
print("\nTable 2: Cell-Level Confounders")
display(table2_df)