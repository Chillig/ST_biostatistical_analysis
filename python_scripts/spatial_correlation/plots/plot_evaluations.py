"""Evaluation metrics
    File name: plot_evaluations.py
    Author: Christina Hillig
    Date created: May/01/2021
    Date last modified: May/01/2021
    Python Version: 3.7
"""

# import scripts
from python_scripts.spatial_correlation import helper_functions
from python_scripts.spatial_correlation import corr_statistics

# Plotting packages
import matplotlib.pyplot as plt
import seaborn as sns

# System specific
import os

# Calculation packages
import scanpy as sc
import numpy as np
import pandas as pd
import random


# Figure params
sc.set_figure_params(color_map='viridis')
fig_size, title_fontsize, axis_label_fontsize, legend_fontsize, fileformat, img_key, xy_ticks, text_fontsize = \
    helper_functions.figure_params()


def plot_evaluate_distance(significance: list, cytokines: list, min_radius: int, save_folder: str, corr_method: str):
    """Elbow plot for best distance/ radius evaluation

    Parameters
    ----------
    significance : list
    cytokines : list
    min_radius: int
    save_folder : str
    corr_method: str

    Returns
    -------

    """
    # Evaluate distance via elbow plot
    significance = np.array(significance).T
    sig_spearman = []
    for val in significance:
        sig_spearman.append(val[corr_method])

    # Correlation p-values
    sig_spearman = np.array(sig_spearman)

    # load cytokine to color
    cyto_color = helper_functions.get_color_signaturegenes()

    for cyto in cytokines:
        if cyto not in cyto_color.keys():
            # Create random hex color for missing keys
            random_number = random.randint(0, 16777215)
            hex_number = str(hex(random_number))
            hex_number = '#' + hex_number[2:]
            cyto_color[cyto] = hex_number

    if min_radius > 0:
        x_vals = np.arange(min_radius, sig_spearman.shape[0] + 1)
    else:
        x_vals = np.arange(min_radius, sig_spearman.shape[0])

    if sig_spearman.T[1:].astype('float').min() < 0:
        ymin_s = -1
    else:
        ymin_s = 0

    # Store p-values and radius in dataframe
    df_colnames = ['radius']
    cyto_correlation_colnames = ['correlation ' + x for x in cytokines]
    cyto_pvalue_colnames = ['-log10(pvalue) ' + x for x in cytokines]
    df_colnames.extend(cyto_correlation_colnames + cyto_pvalue_colnames)
    df_spearman = pd.DataFrame(columns=df_colnames)
    df_spearman['radius'] = x_vals

    fig_pval, ax_pval = plt.subplots(figsize=fig_size)
    # ax_pval.grid(False)
    fig_corr, ax_corr = plt.subplots(figsize=fig_size)
    ax_corr.set_ylim([ymin_s, 1])
    # ax_corr.grid(False)
    for ind, cyto in enumerate(cytokines):
        mask = sig_spearman.T[:, ind, :].T[:, 2].astype('float') < 0.05
        ind_notsigpval = np.where(mask == False)[0]
        ind_sigpval = np.where(mask == True)[0]

        # Save Correlation values to dataframe
        df_spearman["correlation {}".format(cyto)] = sig_spearman.T[:, ind, :].T[:, 1].astype('float')
        df_spearman["-log10(pvalue) {}".format(cyto)] = -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float'))

        ax_pval.plot(x_vals, -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float')),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_pval.scatter(x_vals[mask], -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float'))[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_pval.scatter(x_vals[~mask], -np.log10(sig_spearman.T[:, ind, :].T[:, 2].astype('float'))[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolors='none')
        ax_pval.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_pval.set_ylabel(r'-log$_{10}$(p-values)', fontsize=axis_label_fontsize)
        ax_pval.set_xticks(x_vals)
        sns.despine(ax=ax_pval)

        # Plot correlation vs radius
        ax_corr.plot(x_vals, sig_spearman.T[:, ind, :].T[:, 1].astype('float'),
                     linestyle='-', c=cyto_color[cyto], label=cyto)
        # Highlight significant markers with triangle and non significant ones with unfilled circle
        if len(ind_sigpval) > 0:
            ax_corr.scatter(x_vals[mask], sig_spearman.T[:, ind, :].T[:, 1].astype('float')[mask],
                            linestyle='-', marker='^', c=cyto_color[cyto])
        if len(ind_notsigpval) > 0:
            ax_corr.scatter(x_vals[~mask], sig_spearman.T[:, ind, :].T[:, 1].astype('float')[~mask],
                            linestyle='-', marker='o', c=cyto_color[cyto], facecolor='white')
        ax_corr.set_xlabel('Radius', fontsize=axis_label_fontsize)
        ax_corr.set_ylabel(r'Correlation value', fontsize=axis_label_fontsize)
        ax_corr.set_xticks(x_vals)
        sns.despine(ax=ax_corr)

    leg = ax_pval.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_pval.savefig(os.path.join(save_folder, '{}_Pval_vs_Radius_Evaluation{}'.format(corr_method, fileformat)),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_pval)
    leg = ax_corr.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=False)
    fig_corr.savefig(os.path.join(save_folder, '{}_Corr_vs_Radius_Evaluation{}'.format(corr_method, fileformat)),
                     bbox_inches='tight',  bbox_extra_artists=(leg,))
    plt.close(fig=fig_corr)

    return df_spearman


def plot_responder_vs_radius(counts_dict: dict, conditionalgenes_responders: dict, radii: list, save_folder: str):
    # Normalise Responder counts by by number of clusters
    df_respcounts_radius = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine'])
    for radius in radii:
        for cyto in conditionalgenes_responders.keys():
            df_temp = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine'])
            for ind, specimen in enumerate(counts_dict[radius]['Specimen'].unique()):
                df_temp_specimen = pd.DataFrame(columns=['normed_responder', 'radius', 'cytokine', 'num_clusters'])
                df_specimen = counts_dict[radius][counts_dict[radius].loc[:, 'Specimen'] == specimen].copy()
                # Number of clusters for a specific cytokine on a specimen
                num_clusters = df_specimen[~df_specimen['{}_responder'.format(cyto)].isna()].shape[0]
                numspots_cluster = df_specimen[~df_specimen['{}_responder'.format(cyto)].isna()][
                    'Cluster_num_spots']

                # Responder counts normed by clusters
                if np.count_nonzero(numspots_cluster) > 0:
                    numspots_cluster = numspots_cluster.replace(0, np.nan)
                    normed_resp_counts = df_specimen['{}_responder'.format(cyto)].astype(float) / numspots_cluster
                    normed_resp_counts = normed_resp_counts.sum()
                else:
                    normed_resp_counts = np.nan

                df_temp_specimen.loc[ind, 'normed_responder'] = normed_resp_counts
                df_temp_specimen.loc[ind, 'radius'] = radius
                df_temp_specimen.loc[ind, 'cytokine'] = cyto
                df_temp_specimen.loc[ind, 'num_clusters'] = num_clusters

                df_temp = pd.concat([df_temp, df_temp_specimen])

            df_respcounts_radius = pd.concat([df_respcounts_radius, df_temp], axis=0)

    # Draw radius vs Responder counts normed by clusters
    for cyto in conditionalgenes_responders.keys():
        fig, ax = plt.subplots()
        ax.grid(False)
        sns.pointplot(x="radius", y="normed_responder", ci="sd", capsize=0.1,
                      data=df_respcounts_radius[df_respcounts_radius['cytokine'] == cyto],
                      dodge=True, join=False, ax=ax)
        sns.despine(ax=ax)
        ax.set_xlabel('Radius')
        ax.set_ylabel('Normed {} responder counts'.format(cyto.split('_')[0]))
        plt.tight_layout()
        fig.savefig(os.path.join(save_folder, '{}__Radius_vs_normed_bynspots_Respcounts.pdf'.format(cyto)))
        plt.close(fig=fig)


def plot_responder_in_sdc_outside_l_tissue(adata, conditionalgenes_responders, save_folder):
    #  Transcript of responder genes: SCD vs L, non clustered boxplots
    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange',
                 'Mean_cytokine_responders_in_SDC', 'Mean_cytokine_responders_outside_SDC_lesion'],
        index=list(conditionalgenes_responders.keys()))

    # statistical annotation
    x1, x2 = 0, 1  # columns 'SDC' and 'NL'
    for cyto, r in zip(conditionalgenes_responders.keys(), ('4', '3', '0')):
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

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax)
        # sns.stripplot(x='variable', y='value', data=df_melted)
        ax.set_yscale('log', basey=2)
        ax.set_ylabel('Responder gene counts')
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
        fig.savefig(os.path.join(save_folder, 'Boxplot_{}_r{}_SDC_LNL__vs__L.pdf'.format(cyto, r)))
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


def plot_in_sdc_cytokine_vs_responder(adata, conditionalgenes_responders, save_folder):
    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange', 'Mean_cytokine_in_SDC', 'Mean_cytokine_responders_in_SDC'],
        index=list(conditionalgenes_responders.keys()))

    x1, x2 = 0, 1  # columns 'SDC' and 'NL'
    for cyto, r in zip(conditionalgenes_responders.keys(), ('4', '3', '0')):
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

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(False)
        sns.boxplot(x="variable", y="value", data=df_melted, ax=ax)
        ax.set_yscale('log', basey=10)
        ax.set_ylabel('Gene counts')
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
        fig.savefig(os.path.join(save_folder, 'Boxplot_{}_r{}_SDC_Cytokine__vs__Responders.pdf'.format(cyto, r)))
        plt.close(fig=fig)

        df_statistics.loc[cyto, :] = [p, log2fc, foldchange, mean_cytokine_sdc, mean_responder_sdc]

    return df_statistics
