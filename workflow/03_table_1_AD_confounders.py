import numpy as np
import pandas as pd
import scanpy as sc

#Load multiregion data comparions.
print("loading data...")
adata_HIP=sc.read_h5ad(f"data/SEA-AD_HIP_Processed.h5ad")
adata_MTG=sc.read_h5ad(f"data/SEA-AD_MTG_Processed.h5ad")
adata_DFC=sc.read_h5ad(f"data/SEA-AD_DFC_Processed.h5ad")
brain_regions=['HIP','MTG','DFC']
adatas={
    "HIP":adata_HIP,
    "MTG":adata_MTG,
    "DFC":adata_DFC
}

# TABLE 1: Donor Cohort Demographics (Combined Donors across All Regions)
ad_col = "Overall AD neuropathological Change"
ad_tiers = ["Not AD", "Low", "Intermediate", "High"]

# Extract and combine unique donor metadata across HIP, MTG, and DFC
donor_df = (
    pd.concat(
        [
            adata_obj.obs.drop_duplicates(subset=["Donor ID"])
            for adata_obj in adatas.values()
        ],
        ignore_index=True,
    )
    .drop_duplicates(subset=["Donor ID"])
    .copy()
)

donor_df[ad_col] = pd.Categorical(
    donor_df[ad_col], categories=ad_tiers, ordered=True
)

records_t1 = {}
for tier in ad_tiers:
  sub = donor_df[donor_df[ad_col] == tier]
  n_donors = len(sub)
  if n_donors == 0:
    continue

  age_str = (
      f"{sub['Age at Death'].mean():.1f} ± {sub['Age at Death'].std():.1f}"
      if "Age at Death" in sub
      else "N/A"
  )
  n_female = (
      (sub["Sex"].str.lower() == "female").sum() if "Sex" in sub else 0
  )
  n_male = (sub["Sex"].str.lower() == "male").sum() if "Sex" in sub else 0
  sex_str = f"{n_female} ({(n_female/n_donors)*100:.1f}%) Female, {n_male} Male"

  pmi_str = (
      f"{sub['PMI'].mean():.1f} ± {sub['PMI'].std():.1f}"
      if "PMI" in sub
      else "N/A"
  )
  apoe_counts = (
      sub["APOE Genotype"].value_counts().to_dict()
      if "APOE Genotype" in sub
      else {}
  )
  apoe_str = (
      ", ".join([f"{k}: {v}" for k, v in apoe_counts.items()])
      if apoe_counts
      else "N/A"
  )

  records_t1[tier] = {
      "Donor Count": n_donors,
      "Age at Death": age_str,
      "Sex (% Female, Count)": sex_str,
      "PMI (hours)": pmi_str,
      "APOE Genotypes": apoe_str
  }

table1_df = pd.DataFrame(records_t1).T
table1_df.to_csv("data/Supplementary_Table_1:_Donor_Confounders.csv")
print("\nTable 1: Donor Cohort Demographics")
display(table1_df)