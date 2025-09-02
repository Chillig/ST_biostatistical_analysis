"""Evaluation metrics
    File name: plot_evaluations.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: September/01/2025
    Python Version: 3.8
"""

# import scripts
import utils_plots
import corr_statistics

# Plotting packages
import matplotlib.pyplot as plt
import seaborn as sns

# System specific
import os

# Calculation packages
import numpy as np
import pandas as pd
import random


# Figure params
fig_size, _, axis_label_fontsize, _, _, _, _, _ = utils_plots.figure_params()


def build_radius_vs_correlation_df(dict_corr_genes: dict, conditional_genes: list,
                                   response_type: str, corr_method: str) -> pd.DataFrame:
    """
    Build DataFrame mapping radius → correlation & pvalue & -log10(pvalue) for each gene of interest.

    Parameters
    ----------
    dict_corr_genes : dict
        Format:
            {
                radius: {
                    cytokine: {
                        'pearson': [(gene, gene_responder, correlation_value, pvalue), ...],
                        'spearman': [(gene, gene_responder, correlation_value, pvalue), ...]
                    },
                    ...
                },
                ...
            }
    conditional_genes : list of str
        Gene of interest to include
    response_type : list of str
        Close vicinity genes of each conditional gene
    corr_method : str
        'pearson' or 'spearman'

    Returns
    -------
    pd.DataFrame
    """
    rows = []
    for radius, gene_dict in dict_corr_genes.items():
        row = {"radius": radius}
        for gene in conditional_genes:
            if gene in gene_dict and corr_method in gene_dict[gene]:
                _, _, corr, pval = gene_dict[gene][corr_method][0]
                row["gene"] = gene
                row["responder type"] = response_type
                row["correlation"] = corr
                row["pvalue"] = pval
                row["-log10(pvalue)"] = -np.log10(pval)
                row["method"] = corr_method
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("radius").reset_index(drop=True)
    return df


def compute_responder_counts_by_radius_df(
        counts_dict: dict, conditionalgenes_responders: dict, radii: list) -> pd.DataFrame:
    """
    Compute normalized responder counts across radii and gene of interests.

    Parameters
    ----------
    counts_dict : dict
        Dictionary of DataFrames, keyed by radius. Each DataFrame must contain:
            - 'Specimen'
            - 'Cluster_num_spots'
            - '{gene}_responder' columns
    conditionalgenes_responders : dict
        Mapping of gene of interest to responder gene lists.
    radii : list
        Radii to compute statistics for.

    Returns
    -------
    df_respcounts_radius : pandas.DataFrame
        Normalized responder counts per specimen across radii and cytokines.
        Columns: ['normed_responder', 'radius', 'gene', 'num_clusters']
    """
    df_respcounts_radius = pd.DataFrame(columns=['normed_responder', 'radius', 'gene', 'num_clusters'])

    for radius in radii:
        for gene in conditionalgenes_responders.keys():
            df_temp = pd.DataFrame(columns=['normed_responder', 'radius', 'gene', 'num_clusters'])
            for specimen in counts_dict[radius]['Specimen'].unique():
                df_specimen = counts_dict[radius][counts_dict[radius]['Specimen'] == specimen].copy()

                num_clusters = df_specimen[~df_specimen[f'{gene}_responder'].isna()].shape[0]
                numspots_cluster = df_specimen.loc[~df_specimen[f'{gene}_responder'].isna(), 'Cluster_num_spots']

                if np.count_nonzero(numspots_cluster) > 0:
                    numspots_cluster = numspots_cluster.replace(0, np.nan)
                    normed_resp_counts = (df_specimen[f'{gene}_responder'].astype(float) / numspots_cluster).sum()
                else:
                    normed_resp_counts = np.nan

                df_temp = pd.concat([df_temp, pd.DataFrame([{
                    'normed_responder': normed_resp_counts,
                    'radius': radius,
                    'gene': gene,
                    'num_clusters': num_clusters
                }])])

            df_respcounts_radius = pd.concat([df_respcounts_radius, df_temp], axis=0)

    return df_respcounts_radius.reset_index(drop=True)


def plot_evaluate_distance(df_spearman: pd.DataFrame, cytokines: list, save_folder: str, corr_method: str,
                           min_radius: int = 0):
    """
    Elbow plot for best distance/ radius evaluation
    (expects pre-computed df_radius_vs_correlation)

    Parameters
    ----------
    df_spearman : pd.DataFrame
        Dataframe containing "radius", "gene", "responder type", "correlation", "pvalue", "-log10(pvalue)" columns
    cytokines : list of str
    save_folder : str
    corr_method : str
    min_radius : int, optional
        Used to trim the x-axis start

    Returns
    -------
    None
    """
    # load cytokine to color
    cyto_color = utils_plots.get_color_signaturegenes()
    for ind, cyto in enumerate(cytokines):
        if cyto not in cyto_color.keys():
            random.seed(ind)
            random_number = random.randint(0, 0xFFFFFF)
            cyto_color[cyto] = f'#{random_number:06x}'

    # --- x values (radius)
    x_vals = df_spearman["radius"].values
    if min_radius > 0:
        mask_r = x_vals >= min_radius
        x_vals = x_vals[mask_r]
    else:
        mask_r = np.ones_like(x_vals, dtype=bool)

    # --- set up plots
    fig_pval, ax_pval = plt.subplots(figsize=fig_size)
    fig_corr, ax_corr = plt.subplots(figsize=fig_size)

    ymin_s = -1 if df_spearman.filter(like="correlation").min().min() < 0 else 0
    ax_corr.set_ylim([ymin_s, 1])

    # --- loop over cytokines
    for cyto in cytokines:
        mask_gene = df_spearman['gene'] == cyto
        corr_vals = df_spearman.loc[mask_gene & mask_r, "correlation"].values
        pvals_log = df_spearman.loc[mask_gene & mask_r, "-log10(pvalue)"].values

        mask_sig = pvals_log > -np.log10(0.05)

        # --- p-value plot
        ax_pval.plot(x_vals, pvals_log, '-', c=cyto_color[cyto], label=cyto)
        ax_pval.scatter(x_vals[mask_sig], pvals_log[mask_sig], marker='^', c=cyto_color[cyto])
        ax_pval.scatter(x_vals[~mask_sig], pvals_log[~mask_sig], marker='o', c=cyto_color[cyto], facecolors='none')
        ax_pval.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_pval.set_ylabel(r'-log$_{10}$(p-values)', fontsize=axis_label_fontsize)
        ax_pval.set_xticks(x_vals)
        sns.despine(ax=ax_pval)

        # --- correlation plot
        ax_corr.plot(x_vals, corr_vals, '-', c=cyto_color[cyto], label=cyto)
        ax_corr.scatter(x_vals[mask_sig], corr_vals[mask_sig], marker='^', c=cyto_color[cyto])
        ax_corr.scatter(x_vals[~mask_sig], corr_vals[~mask_sig], marker='o', c=cyto_color[cyto], facecolor='white')
        ax_corr.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_corr.set_ylabel(r'Correlation value', fontsize=axis_label_fontsize)
        ax_corr.set_xticks(x_vals)
        sns.despine(ax=ax_corr)

    # Save plots
    leg = ax_pval.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    fig_pval.savefig(os.path.join(save_folder, f'{corr_method}_Pval_vs_Radius_Evaluation.pdf'),
                     bbox_inches='tight', bbox_extra_artists=(leg,))
    plt.close(fig_pval)

    leg = ax_corr.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    fig_corr.savefig(os.path.join(save_folder, f'{corr_method}_Corr_vs_Radius_Evaluation.pdf'),
                     bbox_inches='tight', bbox_extra_artists=(leg,))
    plt.close(fig_corr)


def plot_responder_vs_radius(df_respcounts_radius: pd.DataFrame, save_folder: str):
    """
    Plot normalized responder counts across radii for each gene of interest.

    Parameters
    ----------
    df_respcounts_radius : pd.DataFrame
        DataFrame returned by `compute_responder_counts_by_radius`.
        Columns: ['normed_responder', 'radius', 'gene', 'num_clusters']
    save_folder : str
        Path to save generated plots.

    Returns
    -------
    None
    """
    for gene in df_respcounts_radius['gene'].unique():
        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.pointplot(
            x="radius", y="normed_responder", errorbar="sd", capsize=0.1,
            data=df_respcounts_radius[df_respcounts_radius['gene'] == gene],
            dodge=True, join=False, ax=ax
        )
        sns.despine(ax=ax)
        ax.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax.set_ylabel(f'Normed {gene.split("_")[0]} responder counts', fontsize=axis_label_fontsize)
        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, f'{gene}__Radius_vs_normed_bynspots_Respcounts.pdf'), bbox_inches='tight')
        plt.close(fig)


def plot_responder_in_sdc_vs_outside_sdc_in_lesional_tissue(adata, conditionalgenes_responders, save_folder):
    """
    Generate boxplots comparing responder gene expression
    inside SDC (spatial density clusters) vs. outside SDC within lesional tissue.

    For each cytokine in `conditionalgenes_responders`, this function:
      - Extracts responder genes from `adata`
      - Aggregates responder counts inside and outside SDC regions
      - Performs statistical testing (Wilcoxon rank-sum test)
      - Calculates effect sizes (log2 fold-change, fold-change, group means)
      - Creates boxplots with statistical annotation
      - Saves plots to disk
      - Collects results into a summary DataFrame

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing expression values (in `.X` or `.layers['counts']`)
        and cell-level metadata in `.obs`. Must include columns of the form:
        ``{cytokine}_in_sdcc_r{replicate}`` and ``biopsy_type``.

    conditionalgenes_responders : dict
        Dictionary mapping cytokine names to lists of associated responder genes.
        Example:
            {
                'IFNG': ['geneA', 'geneB', ...],
                'IL13': ['geneX', 'geneY', ...],
                ...
            }

    save_folder : str
        Path to the folder where generated boxplots will be saved.

    Returns
    -------
    df_statistics : pandas.DataFrame
        Summary statistics for each cytokine, indexed by cytokine name.
        Columns:
            - 'p-value' : Wilcoxon test p-value
            - 'Log2FC'  : log2 fold-change (SDC vs outside SDC in lesional tissue)
            - 'FoldChange' : fold-change (linear scale)
            - 'Mean_cytokine_responders_in_SDC' : mean aggregated responder counts inside SDC
            - 'Mean_cytokine_responders_outside_SDC_lesion' : mean aggregated responder counts outside SDC

    Notes
    -----
    - Boxplots are saved as PDF files, one per cytokine.
    - The y-axis is log-scaled (base 2) to visualize wide expression ranges.
    - Significant differences are annotated with p-values and log2 fold-change,
      while non-significant comparisons are labeled "ns".
    """

    #  Transcript of responder genes: SCD vs L, non clustered boxplots
    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange',
                 'Mean_cytokine_responders_in_SDC', 'Mean_cytokine_responders_outside_SDC_lesion'],
        index=list(conditionalgenes_responders.keys()))

    # statistical annotation
    x1, x2 = 0, 1  # columns 'SDC' and 'NL'
    for cyto, r in zip(conditionalgenes_responders.keys(), '3'):
        mask_2L = (adata.obs['{}_in_sdcc_r{}'.format(cyto, r)] != 0)  # & (adata.obs['biopsy_type'] == 'LESIONAL')
        # Compare SCD vs L, non clustered
        mask_0L = (adata.obs['{}_in_sdcc_r{}'.format(cyto, r)] == 0) & (adata.obs['biopsy_type'] == 'LESIONAL')

        # Get responder genes of signature cytokine
        mask_genes = adata.var_names.isin(conditionalgenes_responders[cyto])
        adata_respondergenes = adata[:, mask_genes].copy()

        df_boxplot = pd.DataFrame(columns=['L'])
        if 'counts' not in adata_respondergenes.layers:
            df_boxplot['L'] = adata_respondergenes.X[mask_0L].sum(axis=1)
            df_scd = pd.DataFrame({'SDC of {}'.format(cyto): adata_respondergenes.X[mask_2L].sum(axis=1)})
            df_boxplot = pd.concat([df_boxplot, df_scd], axis=1)
        else:
            df_boxplot['L'] = adata_respondergenes.layers['counts'][mask_0L].sum(axis=1)
            df_scd = pd.DataFrame({
                'SDC of {}'.format(cyto): adata_respondergenes.layers['counts'][mask_2L].sum(axis=1)})
            df_boxplot = pd.concat([df_boxplot, df_scd], axis=1)

        df_melted = pd.melt(df_boxplot)
        df_melted['value'] = df_melted['value'].astype('float')
        y, h, col = df_melted['value'].max() + 20000, 20000, 'k'

        # Calculate statistical test
        p = corr_statistics.apply_wilcoxontest(
            df_highcounts=df_boxplot[
                ~df_boxplot['SDC of {}'.format(cyto)].isna()]['SDC of {}'.format(cyto)], df_lowcounts=df_boxplot['L'])
        # Calculate log2FC
        # IFNG: SCD NL L vs NL: 2.2317657 (mean) & 2.6739449766382893e-273 (pvalue)
        # IL13 SCD NL L vs NL:  0.77387166 (mean) & 2.7529357754941514e-38 (pvalue)
        # IL17A SCD NL L vs NL:  2.279413 (mean) & 9.077873202626242e-41 (pvalue)
        foldchange = np.nanmean(df_boxplot['SDC of {}'.format(cyto)]) / np.nanmean(df_boxplot['L'])
        log2fc = np.log2(foldchange)
        print("Log2FC {}: ".format(cyto), log2fc)
        # Log2FC IFNG:  2.2317657
        # Log2FC IL13:  0.77387166
        # Log2FC IL17A:  2.2804325
        # Calculate FoldChange
        print("FoldChange {}: ".format(cyto), foldchange)
        # FoldChange IFNG: 4.70
        # FoldChange IL13:  1.7098522
        # FoldChange IL17A:  4.854804

        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax)
        # sns.stripplot(x='variable', y='value', data=df_melted)
        ax.set_yscale('log', base=2)
        ax.set_ylabel('Responder gene counts', fontsize=axis_label_fontsize)
        ax.set_xlabel('')
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)
        # Add significance
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        if p < 0.05:
            plt.text((x1 + x2) * .5, y + h, "p-value = {:.2E}, Log2FC = {:.2f}".format(p, log2fc),
                     ha='center', va='bottom', color=col)
        else:
            plt.text((x1 + x2) * .5, y + h, "ns", ha='center', va='bottom', color=col)
        fig.savefig(os.path.join(save_folder, 'Boxplot_{}_r{}_SDC_LNL__vs__L.pdf'.format(cyto, r)),
                    bbox_inches='tight')
        plt.close(fig=fig)

        mean_responders_in_sdc = np.nanmean(df_boxplot['SDC of {}'.format(cyto)])
        print("Mean of {} responder counts inside all SDC: ".format(cyto), mean_responders_in_sdc)
        # Mean of IFNG responder counts inside all SDC:  37.811535
        # Mean of IL13 responder counts inside all SDC:  18.553255
        # Mean of IL17A responder counts inside all SDC:  2239.6333
        mean_responders_outside_sdc_lesion_skin = np.nanmean(df_boxplot['L'])
        print("Mean of {} responder counts outside SDC: ".format(cyto), mean_responders_outside_sdc_lesion_skin)
        # Mean of IFNG responder counts outside SDC:  8.05
        # Mean of IL13 responder counts outside SDC:  10.850795
        # Mean of IL17A responder counts outside SDC:  461.32312

        df_statistics.loc[cyto, :] = [p, log2fc, foldchange, mean_responders_in_sdc,
                                      mean_responders_outside_sdc_lesion_skin]

    return df_statistics


def plot_in_sdc_cytokine_vs_responder(adata, conditionalgenes_responders, radius, save_folder):
    """
    Compare expression of gene of interest to expression level of its
    hypothesised proximity genes inside SDC (spatial density clusters) and visualize results with boxplots.

    For each gene of interest in `conditionalgenes_responders`, this function:
      - Extracts counts for the gene of interest itself and its associated responder genes
      - Aggregates expression values within SDC cells
      - Performs Wilcoxon rank-sum test to compare gene of interest vs responder gene counts
      - Calculates effect sizes (log2 fold-change, fold-change, group means)
      - Creates boxplots with significance annotation
      - Saves plots as PDF
      - Returns summary statistics as a DataFrame

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing expression values (in `.X` or `.layers['counts']`)
        and metadata in `.obs`. Must include columns of the form:
        ``{gene}_in_sdcc_r{radius}``.

    conditionalgenes_responders : dict
        Dictionary mapping gene of interest names to lists of associated responder genes.
        Example:
            {
                'IFNG': ['geneA', 'geneB', ...],
                'IL13': ['geneX', 'geneY', ...],
                ...
            }

    radius : int
        Radius resulting in the correlation between gene of interest and its hypothesised proximity genes

    save_folder : str
        Path to the folder where generated boxplots will be saved.

    Returns
    -------
    df_statistics : pandas.DataFrame
        Summary statistics for each gene of interest, indexed by gene of interest name.
        Columns:
            - 'p-value' : Wilcoxon test p-value
            - 'Log2FC'  : log2 fold-change (Responder genes vs cytokine)
            - 'FoldChange' : fold-change (linear scale)
            - 'Mean_gene_in_SDC' : mean gene of interest counts across SDC
            - 'Mean_gene_responders_in_SDC' : mean aggregated responder counts across SDC

    Notes
    -----
    - Boxplots are saved as PDF files, one per gene of interest.
    - The y-axis is log-scaled (base 10) to capture wide dynamic ranges.
    - Significant comparisons are annotated with p-values and log2 fold-change,
      while non-significant ones are labeled "ns".
    """

    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange', 'Mean_gene_in_SDC', 'Mean_gene_responders_in_SDC'],
        index=list(conditionalgenes_responders.keys()))

    x1, x2 = 0, 1  # columns 'SDC' and 'NL'
    for cyto, r in zip(conditionalgenes_responders.keys(), str(radius)):
        # Read out spots outside the cluster
        mask_1L = (adata.obs['{}_in_sdcc_r{}'.format(cyto, r)] == 1)
        # Compare within SCD cytokine counts vs responder counts
        mask_2L = (adata.obs['{}_in_sdcc_r{}'.format(cyto, r)] != 0)

        # Get responder genes of signature cytokine
        mask_genes = adata.var_names.isin(conditionalgenes_responders[cyto])
        # Read out counts of responder genes and cytokine
        adata_respondergenes = adata[:, mask_genes].copy()
        adata_cytokine = adata[:, cyto].copy()

        df_boxplot = pd.DataFrame(columns=[cyto])
        if 'counts' not in adata_respondergenes.layers:
            df_boxplot[cyto] = adata_cytokine.X[mask_1L].sum(axis=1)
            df_scd = pd.DataFrame({'Responder genes of {}'.format(cyto): adata_respondergenes.X[mask_2L].sum(axis=1)})
            df_boxplot = pd.concat([df_boxplot, df_scd], axis=1)
        else:
            df_boxplot[cyto] = adata_cytokine.layers['counts'][mask_1L].sum(axis=1)
            df_scd = pd.DataFrame({
                'Responder genes of {}'.format(cyto): adata_respondergenes.layers['counts'][mask_2L].sum(axis=1)})
            df_boxplot = pd.concat([df_boxplot, df_scd], axis=1)

        df_melted = pd.melt(df_boxplot)
        df_melted['value'] = df_melted['value'].astype('float')
        y, h, col = df_melted['value'].max() + 20000, 20000, 'k'

        # Calculate statistical test
        p = corr_statistics.apply_wilcoxontest(
            df_highcounts=df_boxplot[
                ~df_boxplot['Responder genes of {}'.format(cyto)].isna()]['Responder genes of {}'.format(cyto)],
            df_lowcounts=df_boxplot[~df_boxplot[cyto].isna()][cyto])
        # IFNG: 2.6830858106296486e-72
        # IL13: 5.5251136752510094e-18
        # IL17A: 1.9717200682257908e-21

        # Calculate Log2FC and FoldChange
        foldchange = np.nanmean(df_boxplot['Responder genes of {}'.format(cyto)]) / np.nanmean(df_boxplot[cyto])
        log2fc = np.log2(foldchange)
        print("Log2FC {}: ".format(cyto), log2fc)
        # Log2FC IFNG:  4.9874682
        # Log2FC IL13:  4.0549026
        # Log2FC IL17A:  10.467101
        print("LogFC {}: ".format(cyto), foldchange)
        # LogFC IFNG:  31.723236
        # LogFC IL13:  16.620623
        # LogFC IL17A:  1415.5051

        mean_cytokine_sdc = np.nanmean(df_boxplot['Responder genes of {}'.format(cyto)])
        print("Mean of {} responder counts over all SDC: ".format(cyto), mean_cytokine_sdc)
        # Mean of IFNG responder counts over all SDC:  37.811535
        # Mean of IL13 responder counts over all SDC:  18.553255
        # Mean of IL17A responder counts over all SDC:  2241.2166
        mean_responder_sdc = np.nanmean(df_boxplot[cyto])
        print("Mean of {} counts over all SDC: ".format(cyto), mean_responder_sdc)
        # Mean of IFNG counts over all SDC:  1.1919192
        # Mean of IL13 counts over all SDC:  1.1162791
        # Mean of IL17A counts over all SDC:  1.5833334

        fig, ax = plt.subplots(figsize=fig_size)
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax)
        ax.set_yscale('log', base=10)
        ax.set_ylabel('Gene counts', fontsize=axis_label_fontsize)
        ax.set_xlabel('')
        # remove upper and right edge lines in plot
        sns.despine(ax=ax)
        # Add significance
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        if p < 0.05:
            plt.text((x1 + x2) * .5, y + h, "p-value = {:.2E}, Log2FC = {:.2f}".format(p, log2fc),
                     ha='center', va='bottom', color=col)
        else:
            plt.text((x1 + x2) * .5, y + h, "ns", ha='center', va='bottom', color=col)
        fig.savefig(os.path.join(save_folder, 'Boxplot_{}_r{}_SDC_Gene__vs__Responders.pdf'.format(cyto, r)),
                    bbox_inches='tight')
        plt.close(fig=fig)

        df_statistics.loc[cyto, :] = [p, log2fc, foldchange, mean_cytokine_sdc, mean_responder_sdc]

    return df_statistics
