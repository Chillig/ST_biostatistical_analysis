import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import os
from datetime import date


def get_counts(data):
    # number of counts per spot
    data.obs['n_counts'] = data.X.sum(1)
    # number of genes per spot
    data.obs['n_genes'] = (data.X > 0).sum(1)

    data.obs['mean_counts'] = data.obs['n_counts'] / data.obs['n_genes']

    return data


def get_mean_std(data):
    np.std(data.X, axis=1)


def add_disease(data, disease_labels):
    data = data.copy()
    data.obs['disease'] = 'Unknown'
    data.obs['disease'] = data.obs['disease'].astype('<U16')
    for spot_label in disease_labels:
        m_spots = data.obs[spot_label] == 1
        data.obs['disease'][m_spots] = spot_label
    data.obs['disease'] = data.obs['disease'].astype('string').astype('category')

    return data


def plot(data, save_folder, key):
    df = data.obs

    # Plots
    # get mean counts of spots
    _, ax = plt.subplots(figsize=(12, 6))
    ax.set_yscale("log")
    ax = sns.boxplot(x="sample", y="mean_counts", data=df, ax=ax, fliersize=1, linewidth=1, hue="disease")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Mean UMI-counts per spot")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "{}_{}".format(key, "Mean_UMI_counts_per_spot.png")), dpi=300)
    plt.close()

    # mean counts per sample
    _, ax = plt.subplots(figsize=(12, 6))
    ax.set_yscale("log")
    ax = sns.boxplot(x="sample", y="n_counts", data=df, ax=ax, fliersize=1, linewidth=1, hue="disease")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Total UMI-counts per spot")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Put the legend out of the figure
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, "{}_{}".format(key, "Total_UMI_counts_per_spot.png")), dpi=300)
    plt.close()


def main(save_folder):
    adata_wo_4_7 = sc.read(
        '/Users/christina.hillig/Documents/Projects/annData_objects/spatial/2021-07-29/Spatial Transcriptomics_unpp.h5')

    print("Dimension of samples without 4 and 7", adata_wo_4_7.shape)

    # rename LICHEN to LP
    adata_wo_4_7.obs = adata_wo_4_7.obs.rename(columns={"LICHEN": "LP"})
    adata_wo_4_7 = get_counts(data=adata_wo_4_7)
    adata_wo_4_7 = add_disease(data=adata_wo_4_7, disease_labels=['PSO', 'AE', 'LP', 'PRP'])

    # per spot
    df_spot = pd.DataFrame(columns=['sample', 'mean', 'var', 'std'])
    df_spot.loc[:, 'sample'] = adata_wo_4_7.obs['sample'].values
    nan_data = adata_wo_4_7.X.copy()
    # 1. Set zero counts to nan
    nan_data[nan_data == 0] = np.nan
    # 2. Calculate mean, std, and var
    df_spot.loc[:, 'mean'] = np.nanmean(nan_data, axis=1)
    df_spot.loc[:, 'var'] = np.nanvar(nan_data, axis=1)
    df_spot.loc[:, 'std'] = np.nanstd(nan_data, axis=1)

    del nan_data

    # save in .csv file
    df_spot.to_csv(os.path.join(save_folder, "Statistics_per_spot_wo_4_7.csv"))
    df_spot.to_excel(os.path.join(save_folder, "Statistics_per_spot_wo_4_7.xlsx"))

    # TODO per sample
    df_sample = pd.DataFrame(
        columns=['sample', 'min_count_spot', 'max_count_spot', 'mean_count_spot', 'var_count_spot', 'std_count_spot',
                 'median_counts_spot', 'total_genes_detected', 'median_genes_spot'])
    # 1. Sum up counts per sample
    for ind, sample_id in enumerate(adata_wo_4_7.obs['sample'].cat.categories):
        df_sample.loc[ind, 'sample'] = sample_id

        mask_sample = adata_wo_4_7.obs['sample'] == sample_id

        matrix_sample = adata_wo_4_7.X[mask_sample, :].copy()

        sum_over_spots = np.sum(matrix_sample, axis=0)
        sum_over_spots[sum_over_spots > 0] = 1
        df_sample.loc[ind, 'total_genes_detected'] = sum_over_spots.sum()

        #  Calculate mean, std, and var
        # 1. Get counts per sample
        counts_per_sample = np.sum(matrix_sample, axis=0)

        # 2. Set zero counts to nan
        counts_per_sample[counts_per_sample == 0] = np.nan

        # Min and Max count pro sample
        df_sample.loc[ind, 'min_count_spot'] = np.nanmin(counts_per_sample)
        df_sample.loc[ind, 'max_count_spot'] = np.nanmax(counts_per_sample)

        # 3. Calculate mean, std, and var
        mean_count_gene_sample = np.nanmean(counts_per_sample)
        df_sample.loc[ind, 'mean_count_spot'] = mean_count_gene_sample

        var_count_gene_sample = np.nanvar(counts_per_sample)
        df_sample.loc[ind, 'var_count_spot'] = var_count_gene_sample

        std_count_gene_sample = np.nanstd(counts_per_sample)
        df_sample.loc[ind, 'std_count_spot'] = std_count_gene_sample

        # Median UMI Counts per Spot
        sum_over_genes = np.sum(matrix_sample, axis=1)
        median_counts_sample = np.median(sum_over_genes)
        df_sample.loc[ind, 'median_counts_spot'] = median_counts_sample

        # Median Genes per Spot
        matrix_sample[matrix_sample > 0] = 1
        num_genes_spot_sample = np.sum(matrix_sample, axis=1, dtype='int64')
        median_genes_sample = np.median(num_genes_spot_sample)
        df_sample.loc[ind, 'median_genes_spot'] = median_genes_sample

        np.sum(matrix_sample, axis=1, dtype='int64')

    df_sample.to_csv(os.path.join(save_folder, "Statistics_per_sample_wo_4_7.csv"))
    df_sample.to_excel(os.path.join(save_folder, "Statistics_per_sample_wo_4_7.xlsx"))

    plot(data=adata_wo_4_7, save_folder=save_folder, key='wo_4_7')

    del adata_wo_4_7

    adata = sc.read(
        '/Users/christina.hillig/Documents/Projects/annData_objects/spatial/2021-03-19/adata_P15509_P16357_unpp.h5')

    print("Dimension of all samples", adata.shape)

    adata.obs = adata.obs.rename(columns={"LICHEN": "LP"})
    adata = add_disease(data=adata, disease_labels=['PSO', 'AE', 'LP', 'PRP'])
    adata = get_counts(data=adata)
    plot(data=adata, save_folder=save_folder, key='all')


if __name__ == '__main__':
    today = date.today()
    path = os.path.join("..", "..")
    # create saving folder in current project path
    savepath = os.path.join(path, "output", "UMI-counts_per_spots", str(today))
    os.makedirs(savepath, exist_ok=True)

    main(save_folder=savepath)
