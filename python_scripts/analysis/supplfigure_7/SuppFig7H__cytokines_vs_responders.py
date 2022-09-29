"""
Compare expression of cytokines vs their responder genes in lesion skin
1. Read out raw counts from cytokines and responders from lesion skin
2. Calculate mean expression of genes per specimen
3. Create dataframe with two groups cytokines and responders
4. Apply statistical test(s)
5. Calculate Log2FC
6. Plot boxplot
"""


from python_scripts.utils import gene_lists

import os
from datetime import date
import scanpy as sc
import anndata
import numpy as np
import pandas as pd

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt


def main(save_folder: str, pp_st_adata: anndata.AnnData):
    # Read out lesion skin
    pp_st_adata = pp_st_adata[pp_st_adata.obs['biopsy_type'] == 'LESIONAL'].copy()

    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange', 'Mean of cytokines per sample per specimen',
                 'Mean of responders per sample per specimen'], index=[1])

    # 1. Read out cytokines and responder genes only
    cytokines, allinone, cytoresps_dict = gene_lists.get_publication_cyto_resps()
    allinone = np.unique(allinone)
    pp_st_adata = pp_st_adata[:, pp_st_adata.var_names.isin(allinone)].copy()

    # 2. Calculate the mean expression over each specimen individually
    # by doing so, we take into account the number of spots per specimen
    array_mean_geneexpression_perspecimen = np.empty(
        shape=(len(pp_st_adata.obs['specimen'].cat.categories), pp_st_adata.shape[1]))
    for ind, specimen in enumerate(pp_st_adata.obs['specimen'].cat.categories):
        mask_specimen = pp_st_adata.obs['specimen'] == specimen
        if 'counts' in pp_st_adata.layers.keys():
            data_counts = pp_st_adata.layers['counts'][mask_specimen]
        else:
            data_counts = pp_st_adata.X[mask_specimen]
        array_mean_geneexpression_perspecimen[ind] = np.nanmean(data_counts, axis=0)

    # Get mask for cytokines
    mask_cytokines = pp_st_adata.var_names.isin(['IL17A', 'IL13', 'IFNG'])

    # Store Mean expression per specimen of cytokines and all other genes separately
    array_mean_geneexpression_cytokines = array_mean_geneexpression_perspecimen[:, mask_cytokines]
    array_mean_geneexpression_responders = array_mean_geneexpression_perspecimen[:, ~mask_cytokines]

    # Calculate Mean per sample per specimen of cytokines and all other genes separately
    responders_mean_samples_mean_specimen = np.nanmean(array_mean_geneexpression_responders, axis=0)
    cytokines_mean_samples_mean_specimen = np.nanmean(array_mean_geneexpression_cytokines, axis=0)

    # Create data for Boxplot -> for manuscript and reviewer
    df_cytokines = pd.DataFrame({'Cytokines': cytokines_mean_samples_mean_specimen})
    df_others = pd.DataFrame({'Others': responders_mean_samples_mean_specimen})
    print("\nMean in Dataframes: \nCytokines: {} \n Others:{}".format(df_cytokines.mean()[0], df_others.mean()[0]))
    # df_cytokines = pd.DataFrame({'Cytokines': data_cytokines})
    # df_others = pd.DataFrame({'Others': data_responders})
    df_data = pd.concat([df_others, df_cytokines], axis=1)
    # Reformat dataframe
    df_data = df_data.melt()
    df_data = df_data[~df_data['value'].isna()]

    # Test if two groups differ from another
    # 1. Check if data is normal distributed
    print("\nTest for normality")
    shapiro_test_cytokines = stats.shapiro(cytokines_mean_samples_mean_specimen)
    print('shapiro_test_cytokines: ', shapiro_test_cytokines)
    shapiro_test_others = stats.shapiro(responders_mean_samples_mean_specimen)
    print('shapiro_test_cytokines: ', shapiro_test_others)
    # Data is not normal distributed as both p-values < 0.05
    # 2. Check if data has equal variance
    print("\nTest for equal variance")
    stat, p = stats.levene(cytokines_mean_samples_mean_specimen, responders_mean_samples_mean_specimen)
    print("\np-value:", p)
    # The p-value of 0.003 suggests that the groups do not have equal variances
    # 3. Apply non-parametric version of t-test for two samples -> Mann-Whitney test
    res = stats.mannwhitneyu(cytokines_mean_samples_mean_specimen, responders_mean_samples_mean_specimen)
    print("Null Hypothesis that two related paired groups come from the same distribution is rejected: ", res)

    # 4. Calculate Log2FC
    # Calculate Mean over Mean expression (Mean over all specimen)
    mean_cyto_mean_samples_mean_specimen = np.nanmean(cytokines_mean_samples_mean_specimen, axis=0)
    mean_responders_mean_samples_mean_specimen = np.nanmean(responders_mean_samples_mean_specimen, axis=0)
    foldchange = np.divide(mean_responders_mean_samples_mean_specimen, mean_cyto_mean_samples_mean_specimen)
    # Get values for manuscript
    log2fc = np.log2(foldchange)
    print('Log2FC between average expression of cytokines and responder genes:', log2fc)
    print('Average expression of cytokines: {:.2E}'.format(mean_cyto_mean_samples_mean_specimen))
    print('Average expression of responder genes: {:.2E}'.format(mean_responders_mean_samples_mean_specimen))
    print('Fold change between average expression of cytokines and responder genes:', foldchange)

    # 5. Plot Boxplot
    # Draw Boxplots comparing the mean expression of cytokines against responder genes
    x1, x2 = 0, 1  # columns
    y, h, col = df_data['value'].max() + 1000, 1000, 'k'
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.grid(zorder=0)
    sns.boxplot(data=df_data, y='value',  x='variable', ax=ax, zorder=3)
    # Add significance asterix
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], linewidth=1, color='k')
    ax.text((x1 + x2) * .5, y + h + 200, "p-value = {:.2E}, Log2FC = {:.2f}".format(res.pvalue, log2fc),
            ha='center', va='bottom', color=col)
    sns.despine(ax=ax, fig=fig)
    ax.set_yscale('log')
    ax.set_ylabel('Mean expression of genes\nper sample and per specimen')
    ax.set_xlabel('Groups')
    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'Boxplots__Lesion_Mean_expression_cytokines_responders_over_specimen.pdf'))
    fig.savefig(os.path.join(save_folder, 'Boxplots__Lesion_Mean_expression_cytokines_responders_over_specimen.png'),
                dpi=300)
    plt.close(fig=fig)

    df_statistics.loc[1, :] = [res.pvalue, log2fc, foldchange,
                               mean_cyto_mean_samples_mean_specimen, mean_responders_mean_samples_mean_specimen]

    return df_statistics


if __name__ == '__main__':
    today = date.today()
    project_folder = os.path.join("..", "..", "..")
    output_path = os.path.join(
        project_folder, "output", "reviewers", 'SuppFig__cytokines_vs_responders', str(today))
    os.makedirs(output_path, exist_ok=True)

    # Read out preprocessed data
    adata_path = os.path.join(project_folder, "adata_storage")
    date_st_pp = '2022-04-08'  # "2020-12-04" -> 2020-12-04_Visium_Data_QC_BC_clustered.h5
    adata_pp = sc.read(os.path.join(adata_path, date_st_pp, 'st_QC_normed_BC_project_PsoADLP.h5'))
    adata_pp.obs.loc[(adata_pp.obs['basal EPIDERMIS'] == 1) & (adata_pp.obs['DERdepth1'] == 1),
                     'basal EPIDERMIS'] = [0, 0, 1]
    adata_pp.obs.loc[(adata_pp.obs['basal EPIDERMIS'] == 1) & (adata_pp.obs['DERdepth1'] == 1), 'DERdepth1'] = 0

    main(save_folder=output_path, pp_st_adata=adata_pp)
