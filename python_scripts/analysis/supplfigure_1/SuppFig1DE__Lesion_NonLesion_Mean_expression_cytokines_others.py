"""
This script produces Figure 2XX.
It compares the mean expression of all genes against the mean expression of cytokines over all specimen.
By doing so we take into account the different number of spots per specimen
"""
import os
from datetime import date
import scanpy as sc
import anndata
import numpy as np
import pandas as pd

from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt


def calculate_mean_expression_over_specimen(adata):
    array_mean_geneexpression_perspecimen = np.empty(shape=(len(adata.obs['specimen'].cat.categories), adata.shape[1]))

    for ind, specimen in enumerate(adata.obs['specimen'].cat.categories):
        mask_specimen = adata.obs['specimen'] == specimen
        if 'counts' in adata.layers.keys():
            data_counts = adata.layers['counts'][mask_specimen]
        else:
            data_counts = adata.X[mask_specimen]
        array_mean_geneexpression_perspecimen[ind] = np.nanmean(data_counts, axis=0)

    return array_mean_geneexpression_perspecimen


def split_mean_expression_group(adata, mean_geneexpression_perspecimen):
    # Get mask of cytokines
    mask_cytokines = adata.var_names.isin(['IL17A', 'IL13', 'IFNG'])
    # Subset data
    mean_geneexpression_cytokines = mean_geneexpression_perspecimen[:, mask_cytokines]
    mean_geneexpression_others = mean_geneexpression_perspecimen[:, ~mask_cytokines]

    return mean_geneexpression_cytokines, mean_geneexpression_others


def main(save_folder: str, pp_st_adata: anndata.AnnData):
    biopsy_types = ['LESIONAL', 'NON LESIONAL']

    # Calculate the mean expression over each specimen individually, split by biopsy_type
    # by doing so, we take into account the number of spots per specimen
    df_statistics = pd.DataFrame(
        columns=['p-value', 'Log2FC', 'FoldChange', 'Mean of cytokines per sample per specimen',
                 'Mean of other genes per sample per specimen'],
        index=biopsy_types)

    for biopsy_type in biopsy_types:
        print("\n================== Biopsy Type: {} ==================".format(biopsy_type))
        mask = pp_st_adata.obs['biopsy_type'] == biopsy_type
        pp_st_adata_temp = pp_st_adata[mask].copy()

        array_mean_geneexpression_perspecimen = calculate_mean_expression_over_specimen(adata=pp_st_adata_temp)

        # Store Mean expression of cytokines and all other genes separately
        array_mean_geneexpression_cytokines, array_mean_geneexpression_others = split_mean_expression_group(
            adata=pp_st_adata_temp, mean_geneexpression_perspecimen=array_mean_geneexpression_perspecimen)

        # Calculate Mean over Mean expression
        others_mean_samples_mean_specimen = np.nanmean(array_mean_geneexpression_others, axis=0)
        cytokines_mean_samples_mean_specimen = np.nanmean(array_mean_geneexpression_cytokines, axis=0)
        num_genes_above_maxvalue_cytokines = len(
            others_mean_samples_mean_specimen[others_mean_samples_mean_specimen >
                                              np.nanmax(cytokines_mean_samples_mean_specimen)])

        # Get values for manuscript abstract
        num_genes = len(others_mean_samples_mean_specimen)
        mean_cyto_mean_samples_mean_specimen = cytokines_mean_samples_mean_specimen.mean()
        mean_others_mean_samples_mean_specimen = others_mean_samples_mean_specimen.mean()
        foldchange = np.divide(mean_others_mean_samples_mean_specimen, mean_cyto_mean_samples_mean_specimen)
        log2fc = np.log2(foldchange)
        print('Percentage of data with larger values than cytokines:',
              num_genes_above_maxvalue_cytokines / num_genes * 100)
        print('Log2FC between average expression of cytokines and other genes:', log2fc)
        print('Average expression of cytokines: {:.2E}'.format(mean_cyto_mean_samples_mean_specimen))
        print('Average expression of other genes: {:.2E}'.format(mean_others_mean_samples_mean_specimen))
        print('Fold change between average expression of cytokines and other genes:', foldchange)

        # Create data for Boxplot -> for manuscript and reviewer
        df_cytokines = pd.DataFrame({'Cytokines': cytokines_mean_samples_mean_specimen})
        df_others = pd.DataFrame({'Others': others_mean_samples_mean_specimen})
        print("\nMean in Dataframes: \nCytokines: {} \n Others:{}".format(df_cytokines.mean()[0], df_others.mean()[0]))
        # df_cytokines = pd.DataFrame({'Cytokines': data_cytokines})
        # df_others = pd.DataFrame({'Others': data_others})
        df_data = pd.concat([df_others, df_cytokines], axis=1)
        # Reformat dataframe
        df_data = df_data.melt()
        df_data = df_data[~df_data['value'].isna()]

        # Test if two groups differ from another
        # 1. Check if data is normal distributed
        print("\nTest for normality")
        shapiro_test_cytokines = stats.shapiro(cytokines_mean_samples_mean_specimen)
        print('shapiro_test_cytokines: ', shapiro_test_cytokines)
        shapiro_test_others = stats.shapiro(others_mean_samples_mean_specimen)
        print('shapiro_test_others: ', shapiro_test_others)
        # Data is not normal distributed as both p-values < 0.05
        # 2. Check if data has equal variance
        print("\nTest for equal variance")
        stat, p = stats.levene(cytokines_mean_samples_mean_specimen, others_mean_samples_mean_specimen)
        print("p-value:", p)
        # The p-value of above 0.05 suggests that the groups have equal variances
        # 3. Apply non-parametric version of t-test for two samples -> Mann-Whitney test
        res = stats.mannwhitneyu(cytokines_mean_samples_mean_specimen, others_mean_samples_mean_specimen)
        print("\nNull Hypothesis that two related paired groups come from the same distribution is rejected: ", res)

        # Draw Boxplots comparing the mean expression of cytokines against all ohther genes
        x1, x2 = 0, 1  # columns
        y, h, col = df_data['value'].max() + 1000, 1000, 'k'
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.grid(zorder=0)
        sns.boxplot(data=df_data, y='value',  x='variable', ax=ax, zorder=3)
        # Add significance asterix
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], linewidth=1, color='k')
        ax.text((x1 + x2) * .5, y + h, "p-value = {:.2E}, Log2FC = {:.2f}".format(res.pvalue, log2fc),
                ha='center', va='bottom', color=col)
        sns.despine(ax=ax, fig=fig)
        ax.set_yscale('log')
        ax.set_ylabel('Mean expression of genes\nper sample and per specimen')
        ax.set_xlabel('Groups')
        plt.tight_layout()
        fig.savefig(
            os.path.join(
                save_folder, 'Boxplots__{}_Mean_expression_cytokines_others_over_specimen.pdf'.format(biopsy_type)))
        fig.savefig(
            os.path.join(
                save_folder, 'Boxplots__{}_Mean_expression_cytokines_others_over_specimen.png'.format(biopsy_type)),
            dpi=300)
        plt.close(fig=fig)

        df_statistics.loc[biopsy_type, :] = [
            res.pvalue, log2fc, foldchange, mean_cyto_mean_samples_mean_specimen,
            mean_others_mean_samples_mean_specimen]

    return df_statistics


if __name__ == '__main__':
    today = date.today()
    project_folder = os.path.join("..", "..", "..")
    output_path = os.path.join(
        project_folder, "output", "reviewers",
        'SuppFig1A__Lesion_NonLesion_Mean_expression_cytokines_others', str(today))
    os.makedirs(output_path, exist_ok=True)

    # Read out preprocessed data
    adata_path = os.path.join(project_folder, "adata_storage")
    date_st_pp = '2022-04-08'  # "2020-12-04" -> 2020-12-04_Visium_Data_QC_BC_clustered.h5
    adata_pp = sc.read(os.path.join(adata_path, date_st_pp, 'st_QC_normed_BC_project_PsoADLP.h5'))
    adata_pp.obs.loc[(adata_pp.obs['basal EPIDERMIS'] == 1) & (adata_pp.obs['DERdepth1'] == 1),
                     'basal EPIDERMIS'] = [0, 0, 1]
    adata_pp.obs.loc[(adata_pp.obs['basal EPIDERMIS'] == 1) & (adata_pp.obs['DERdepth1'] == 1), 'DERdepth1'] = 0

    main(save_folder=output_path, pp_st_adata=adata_pp)
