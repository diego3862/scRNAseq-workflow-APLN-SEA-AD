#pseudobulk.py
import pathlib
import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import gc
from src.data_utils import pseudobulk_cpm, prepare_data_multiple

#ESTABLISH DIRECTORY
DIRECTORY = "/content/drive/MyDrive/Colab_Notebooks/APJ_project_working_directory"

PROJECT = pathlib.Path(DIRECTORY)
os.chdir(PROJECT)
if PROJECT not in sys.path:
  sys.path.append(PROJECT)

#CHOOSE BRAIN REGION
options=[]
for file in os.listdir(DIRECTORY+"/data/"):
    if str(file[0:7]) =="SEA-AD_":
        options.append(str(file[7:10]))

print(f"pseudobulk processing the following sliced region files'{options}'")

#VARIABLES & SETUP
donor_col = 'Donor ID'
subtype_custom = 'APJ_Cats'
genes_of_interest = ['APLNR','APLN']

#READ ADATA
for region in options:
    print(f"loading processed adata... data/SEA-AD_{region}_Processed.h5ad")
    adata=sc.read_h5ad(f"data/SEA-AD_{region}_Processed.h5ad")

    #PSEUDOBULK GENERAL
    print("processing general pseudobulk")
    pdata_general = pseudobulk_cpm(adata, donor_col, None)
    df_general = prepare_data_multiple(pdata_general, genes_of_interest)
    print("saving APLN and APLNR plotting df to CSV")
    df_general.to_csv(f"data/pseudobulk_{region}_general", index=False)

    del pdata_general, df_general
    gc.collect()

    #PSEUDOBULK SUBTYPES
    print("Processing subtype pseudobulk")
    pdata_subtype=pseudobulk_cpm(adata,donor_col,subtype_custom)

    df_subtype = prepare_data_multiple(pdata_subtype, genes_of_interest)
    df_subtype.to_csv(f"data/pseudobulk_{region}_subtype", index=False)
    del pdata_subtype, df_subtype
    gc.collect()