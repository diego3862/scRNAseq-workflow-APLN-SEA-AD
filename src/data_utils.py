import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse
from scipy.stats import norm
import time
import decoupler as dc
import rapids_singlecell as rsc
import gseapy as gp

#used in most scripts

def assign_APJ_Cats(adata):

    #ASSIGN CUSTOM CELL TYPE CATEGORIES BASED ON AUTHOR SUBTYPES.
    class_col = 'Class' if 'Class' in adata.obs.columns else 'class'
    subclass_col = 'Subclass' if 'Subclass' in adata.obs.columns else 'subclass'

    # Define Boolean conditions
    cond_neuronal = adata.obs[class_col].isin(["Neuronal: Glutamatergic", "Neuronal: GABAergic"])
    cond_microglia = adata.obs[subclass_col] == "Immune"
    cond_astrocyte = adata.obs[subclass_col] == "Astrocyte"

    conditions = [
        cond_neuronal,
        cond_microglia,
        cond_astrocyte
    ]

    #Define corresponding category labels
    choices = [
        "Neuronal",
        "Microglia",
        "Astrocyte"
    ]

    #Assign to APJ_Cats ('Other' for any remaining cell types like Oligodendrocytes/Endothelial)
    adata.obs['APJ_Cats'] = np.select(conditions, choices, default="Other")

    #Convert to category for Scanpy plotting compatibility
    adata.obs['APJ_Cats'] = adata.obs['APJ_Cats'].astype('category')

#Used in 01_pseudobulk_overall.py

def pseudobulk_cpm(data, donor_col, group,return_cpm=True):
    # Create pseudobulk CPM from cell-level UMIs and log-transform
    pdata = dc.pp.pseudobulk(
        data,
        sample_col=donor_col,
        groups_col=group,
        layer='UMIs',
        mode='sum'
    )
    
    if return_cpm:
        sc.pp.normalize_total(pdata, target_sum=1e6)
        sc.pp.log1p(pdata)
        
    return pdata

def prepare_data_multiple(pdata, genes):
    # Extract expression for multiple genes into a dataframe; warns if gene missing
    plot_df = pdata.obs.copy()
    for gene in genes:
        if gene in pdata.var_names:
            gene_idx = pdata.var_names.get_loc(gene)
            if scipy.sparse.issparse(pdata.X):
                plot_df[gene] = pdata.X[:, gene_idx].toarray().flatten()
            else:
                plot_df[gene] = pdata.X[:, gene_idx].flatten()
        else:
            print(f"Warning: {gene} not found in object.")
    return plot_df


#used in 02_HIP_A_N_DE_GSEA.py

def run_wilcoxon(adata, design_factor, reference_level, test_level=None):
    if test_level is None:
        test_level = [
            x
            for x in adata.obs[design_factor].dropna().unique()
            if x != reference_level
        ][0]
    rsc.get.anndata_to_GPU(adata)
    # GPU-accelerated Wilcoxon rank-sum test
    rsc.tl.rank_genes_groups(
        adata, groupby=design_factor, reference=reference_level, method="wilcoxon"
    )
    rsc.get.anndata_to_CPU(adata)

    # Extract, rename, and format to expected output
    df = (
        sc.get.rank_genes_groups_df(adata, group=test_level)
        .set_index("names")
        .rename(
            columns={
                "logfoldchanges": "log2FoldChange",
                "scores": "stat",
                "pvals": "pvalue",
                "pvals_adj": "padj",
            }
        )
    )

    df["baseMean"] = np.nan
    return df[["baseMean", "log2FoldChange", "stat", "pvalue", "padj"]]

def get_signed_enrichment_long(gene_dict, libraries, p_cutoff=0.05, top_n=15):
    #Generates signed (based on +/- fold change) enrichment list for up/downregulation
    """
    Expects gene_dict formatted as: { ('Subtype/Group', 'Up'/'Down'): ['gene1', 'gene2'] }
    """
    results = []
    for (subtype, direction), genes in gene_dict.items():
        if len(genes) < 5: continue
        
        enr = gp.enrichr(gene_list=genes, gene_sets=libraries, organism='homo sapiens').results
        time.sleep(2)  
        enr = enr[enr['Adjusted P-value'] < p_cutoff].copy()
        if enr.empty: continue
        
        enr['Subtype'] = subtype
        sign = 1 if direction == 'Up' else -1
        enr['Signed_P'] = -np.log10(enr['Adjusted P-value'].replace(0, 1e-300)) * sign
        enr['Gene_Ratio'] = enr['Overlap'].apply(lambda x: float(x.split('/')[0])/float(x.split('/')[1]))
        results.append(enr)
            
    if not results: return None
    full_df = pd.concat(results)
    full_df['abs_P'] = full_df['Signed_P'].abs()
    full_df = (full_df.sort_values('abs_P', ascending=False)
               .drop_duplicates(subset=['Term', 'Subtype'])
               .drop(columns=['abs_P']))
    
    top_terms = full_df.sort_values('Adjusted P-value').groupby('Subtype').head(top_n)['Term']
    return full_df[full_df['Term'].isin(top_terms)]
