import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#used in figures.ipynb

def set_paper_theme(font_size=10):
    # Configure global plotting theme and return a custom palette
    sc.set_figure_params(dpi=200, frameon=False, facecolor='white', fontsize=font_size)
    F_TITLE, F_LABEL, F_TICK, F_ANNOTATION = 20, 16, 14, 12

    plt.rcParams.update({
        'font.size': F_LABEL,
        'axes.titlesize': F_TITLE,
        'axes.titleweight': 'bold',
        'axes.labelsize': F_LABEL,
        'axes.linewidth': 2,
        'xtick.labelsize': F_TICK,
        'ytick.labelsize': F_TICK,
        'legend.title_fontsize': F_LABEL,
        'legend.fontsize': F_TICK,
        'figure.titlesize': F_TITLE,
        'figure.titleweight': 'bold',
        'svg.fonttype': 'none'
    })
    sns.set_context("talk", font_scale=1.0)
    sns.set_style("ticks")

def finalize_panel(filename, fig=None, save=True, output_dir="../results/figures/"):
    # Save or show the current figure; exports SVG when saving
    if fig is None:
        fig = plt.gcf()

    if save:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, f"{filename}.svg")
        plt.savefig(path, format='svg', bbox_inches='tight', transparent=True)
        plt.close()
        print(f"Exported: {path}")
    else:
        plt.show()

def plot_pseudobulk_expression_category(
    df,
    genes_of_interest,
    strat_category,
    title,
    region,
    FIG_DIR,
    cell_type=None,
    subtype_col=None,
    order=None,
    save_status=False,
):
    #plot boxplots where 0 expression does not influence boxplot, and % expressing is displayed above. Allows clean comparison of distribution of expresing donors only.
    #for category comparison across 2 genes.
    plot_df = (
        df[df[subtype_col] == cell_type].copy()
        if cell_type and subtype_col in df
        else df.copy()
    )
    g1, g2 = genes_of_interest[:2]
    subtypes = order or list(plot_df[strat_category].dropna().unique())

    x = np.arange(len(subtypes))
    w, off = 0.30, 0.18
    c1, c2 = "#4C72B0", "#DD8452"
    rng = np.random.default_rng(42)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for i, cat in enumerate(subtypes):
        sub = plot_df[plot_df[strat_category] == cat]
        for g, col, pos in [
            (g1, c1, x[i] - off),
            (g2, c2, x[i] + off),
        ]:
            vals = sub[g].dropna()
            pct = (vals > 0).mean() * 100 if len(vals) else 0
            pos_vals = vals[vals > 0]

            # Boxplot based on only >0 expression
            if len(pos_vals):
                ax.boxplot(
                    pos_vals,
                    positions=[pos],
                    widths=w,
                    patch_artist=True,
                    showfliers=False,
                    boxprops=dict(facecolor=col, alpha=0.75),
                    medianprops=dict(color="black", lw=1.5),
                )

            # Dots on ALL values (including zeros)
            if len(vals):
                jitter = rng.uniform(-0.06, 0.06, size=len(vals))
                ax.scatter(pos + jitter, vals, color="black", alpha=0.6, s=20, zorder=3)

            # Text overlay if < 100%
            if pct < 100:
                ax.text(
                    pos,
                    0.92,
                    f"{pct:.0f}%",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    fontsize=9,
                    fontweight="bold",
                    color=col,
                )

    # Shared y-axis limits calculated across both genes
    top = max(plot_df[g1].max(), plot_df[g2].max())
    top = top if top > 0 else 10
    ax.set_ylim(0, top * 1.15)
    ax.set_ylabel("Expression (log1p CPM)", fontweight="bold")

    ax.set(
        xticks=x,
        xticklabels=subtypes,
        xlim=(-0.5, len(subtypes) - 0.5),
        xlabel="Subtype" if strat_category == "APJ_Cats" else strat_category,
    )
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)  # Keeps gridlines behind the boxes and dots
    ax.legend(
        handles=[
            mpatches.Patch(color=c1, label=g1),
            mpatches.Patch(color=c2, label=g2),
        ],
        loc="upper right",
        frameon=True,
    )

    plt.title(title)
    plt.tight_layout()
    finalize_panel(
        f"{region}_{strat_category}_{g1}_{g2}_single_axis_panel",
        output_dir=FIG_DIR,
        save=save_status,
    )
 
def plot_gene_across_regions_by_ad(
    dfs_by_region,
    gene,
    cell_type=None,
    subtype_col=None,
    ad_col="Overall AD neuropathological Change",
    ad_order=None,
    palette="Set2",
    figsize=(8, 5.5),
):
    #plot boxplots where 0 expression does not influence boxplot, and % expressing is displayed above. Allows clean comparison of distribution of expresing donors only.
    #for continuous comparison for one gene across multiple categories 
    ad_order = ad_order or ["Not AD", "Low", "Intermediate", "High"]
    regions = list(dfs_by_region.keys())

    #Combine regions (filter by cell_type only if provided)
    plot_df = pd.concat(
        [
            (
                df[df[subtype_col] == cell_type]
                if cell_type and subtype_col and subtype_col in df
                else df
            ).assign(Region=reg)
            for reg, df in dfs_by_region.items()
        ],
        ignore_index=True,
    )

    pcts = (
        plot_df.groupby([ad_col, "Region"], observed=True)[gene]
        .apply(lambda s: (s > 0).mean() * 100)
        .to_dict()
    )

    fig, ax = plt.subplots(figsize=figsize)
    pos_df = plot_df[plot_df[gene] > 0]

    #Boxplot on non-zeros only
    sns.boxplot(
        data=pos_df,
        x=ad_col,
        y=gene,
        hue="Region",
        order=ad_order,
        hue_order=regions,
        palette=palette,
        showfliers=False,
        ax=ax,
        boxprops=dict(alpha=0.75),
    )

    #Dots for all points (including zeros)
    sns.stripplot(
        data=plot_df,
        x=ad_col,
        y=gene,
        hue="Region",
        order=ad_order,
        hue_order=regions,
        dodge=True,
        color="black",
        alpha=0.6,
        jitter=0.2,
        s=4,
        ax=ax,
        legend=False,
    )

    #Text overlay if < 100% expressing
    w = 0.8 / len(regions)
    for i, stage in enumerate(ad_order):
        for j, reg in enumerate(regions):
            pct = pcts.get((stage, reg), 0)
            if pct < 100:
                ax.text(
                    i - 0.4 + (w / 2) + (j * w),
                    0.92,
                    f"{pct:.0f}%",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    fontsize=8.5,
                    fontweight="bold",
                    color="#333",
                )

    #Limits, labels & legend
    top = plot_df[gene].max() * 1.15 if plot_df[gene].max() > 0 else 10
    title_prefix = f"{cell_type} - " if cell_type else "Overall - "
    ax.set(
        ylim=(0, top),
        xlabel="Overall AD Neuropathological Change",
        ylabel=f"{gene} Expression (log1p CPM)",
        title=f"{title_prefix}{gene} Across Regions",
    )
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)  # Keeps gridlines behind the boxes and dots
    ax.legend(title="Region", frameon=True, loc="lower right")
    plt.tight_layout()
    return fig, ax

def plot_GSEA_dotplot(GENE, SUBTYPES, df, FIG_DIR='figures', save=False):
    """
    GENE: String (e.g., 'APLN' or 'APLNR')
    SUBTYPES: List of categories to plot (e.g., ['APLN'] or ['APLNR'])
    df: Pandas DataFrame from get_signed_enrichment_long
    """
    # Filter to requested subtypes and clean
    df = df[df['Subtype'].isin(SUBTYPES)].dropna(subset=['Term', 'Overlap']).copy()
    
    if df.empty:
        print(f"No data to plot for {GENE}.")
        return

    df['Count'] = df['Overlap'].astype(str).str.split('/').str[0].astype(float)

    # Compute Signed_P only if not already present
    if 'Signed_P' not in df.columns:
        sign = np.where(df.get('direction', 'Up') == 'Down', -1, 1)
        df['Signed_P'] = -np.log10(df['Adjusted P-value'].replace(0, 1e-300)) * sign

    # Sort pathways vertically (Down at bottom, Up at top)
    order = df.groupby('Term')['Signed_P'].mean().sort_values().index.tolist()
    df['Term'] = pd.Categorical(df['Term'], categories=order)
    df['Subtype'] = pd.Categorical(df['Subtype'], categories=SUBTYPES)

    # Plot using numeric codes
    fig, ax = plt.subplots(figsize=(4, len(order) * 0.5))
    vmax = np.nanpercentile(df['Signed_P'].abs(), 98)

    sc = ax.scatter(
        df['Subtype'].cat.codes, 
        df['Term'].cat.codes,
        s=df['Count']*10,
        c=df['Signed_P'],
        cmap='coolwarm', vmin=-vmax, vmax=vmax,
        edgecolors='black', linewidth=0.5
    )

    # Axis ticks & labels
    ax.set_xticks(range(len(SUBTYPES)))
    ax.set_xticklabels(SUBTYPES, rotation=90, ha='right')
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlim(-0.5, len(SUBTYPES) - 0.5)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_title(f'{GENE} Enriched Pathways')

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label(r'Signed $-\log_{10}(\mathrm{Adj\ } P)$')

    # Size legend
    for c in [5, 10, 20]:
        ax.scatter([], [], s=c*10, c='black', label=str(c))
    ax.legend(title="Gene Count", bbox_to_anchor=(1, 0.05), loc='upper left', frameon=False)

    # Hand off to finalize_panel
    finalize_panel(f"{GENE}_GSEA_dotplot", output_dir=FIG_DIR, save=save)