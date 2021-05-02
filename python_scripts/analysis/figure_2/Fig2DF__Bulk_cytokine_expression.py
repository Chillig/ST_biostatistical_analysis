#!/usr/bin/env python
"""Read out of cytokine counts in RNA-seq data
    File name: Fig2DF__Bulk_cytokine_expression.py
    Author: Christina Hillig
    Date created: 3/11/2020
    Date last modified: 3/21/2021
    Python Version: 3.7
"""

from datetime import date
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

fig_size = (8, 8)
xy_fontsize = 16
xy_ticks = 12
title_fontsize = 18
legend_fontsize = 14
text_fontsize = 14
fileformat = '.pdf'


def get_counts(data, metadata, cytokine, diagnosis):
    """

    Parameters
    ----------
    data : pandas.Dataframe
    metadata : pandas.Dataframe
    cytokine: str
    diagnosis : list or str

    Returns
    -------

    """
    if isinstance(diagnosis, list):
        m_diagnosis = np.where(np.array(metadata['diagnosis'])[:, np.newaxis] == np.array(diagnosis)[np.newaxis, :])[0]
    else:
        m_diagnosis = metadata['diagnosis'] == diagnosis
    biopsies_names = data.columns[m_diagnosis]
    diagnosis_data = data[biopsies_names]

    df = pd.DataFrame(columns=['counts', 'skin'])
    df['counts'] = diagnosis_data.loc[cytokine].values
    df['skin'] = metadata['skin'].values[m_diagnosis]

    return df


def plot_count_distribution(counts, cytokine, save_folder, key, title):
    """Plot count distribution of RNA-seq data

    Parameters
    ----------
    counts : pandas.Dataframe
    cytokine : str
    save_folder : str
    key : str
    title : str

    Returns
    -------

    """
    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed')
    # Lesional
    m_lesional = counts['skin'] == 'lesional'
    num_biopsies_lesional = len(counts['counts'].values[m_lesional])
    ax.scatter(x=np.arange(0, 2 * num_biopsies_lesional, 2), y=counts['counts'].values[m_lesional], s=8,
               c='k', label='lesional')
    # Non Lesional
    m_nonlesional = counts['skin'] == 'non-lesional'
    num_biopsies_nonlesional = len(counts['counts'].values[m_nonlesional])
    ax.scatter(x=np.arange(1, 2 * num_biopsies_nonlesional, 2), y=counts['counts'].values[m_nonlesional], s=8,
               c='grey', label='Non-lesional')
    ax.set_xlabel('Biopsy', fontsize=xy_fontsize)
    if cytokine == "IFNG":
        ax.set_ylabel(r'IFN$\gamma$ counts', fontsize=xy_fontsize)
    else:
        ax.set_ylabel(" ".join([cytokine, 'counts']), fontsize=xy_fontsize)
    ax.tick_params(labelsize=xy_ticks)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.annotate("Total No. biopsies: {}".format(int(len(counts.values))), xy=(0.68, 1.05),
                xycoords='axes fraction')
    ax.legend(loc='best', title_fontsize=title_fontsize, fontsize=legend_fontsize)

    plt.savefig(os.path.join(save_folder, "_".join([key, title, fileformat])))
    plt.close()


def main(save_folder, bulk_rnaseq, metadata):
    """Read out counts of cytokines in RNA-seq data and create PCA and UMAP plots

    Parameters
    ----------
    save_folder : str
    bulk_rnaseq : pandas.Dataframe
    metadata : pandas.Dataframe

    Returns
    -------

    """
    # Read out IL17A per pso biopsy
    il17a_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IL17A', diagnosis='psoriasis')

    # Read out IFNG per LICHEN biopsy
    ifng_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IFNG', diagnosis='lichen ruber')

    # Read out IL13 per AE biopsy
    il13_counts = get_counts(data=bulk_rnaseq, metadata=metadata, cytokine='IL13', diagnosis='Atopic eczema')

    # Read out all cytokines per PRP biospy
    il17a_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                  cytokine='IL17A', diagnosis='pytiriasis rubra pilaris')
    ifng_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                 cytokine='IFNG', diagnosis='pytiriasis rubra pilaris')
    il13_prp_counts = get_counts(data=bulk_rnaseq, metadata=metadata,
                                 cytokine='IL13', diagnosis='pytiriasis rubra pilaris')

    # Save as .csv file
    il17a_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PSO__IL17A.csv"))
    ifng_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__LICHEN__IFNG.csv"))
    il13_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__AD__IL13.csv"))
    # PRP
    il17a_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IL17A.csv"))
    ifng_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IFNG.csv"))
    il13_prp_counts.to_csv(os.path.join(save_folder, "Bulk_RNAseq__PRP__IL13.csv"))

    print("Plot")
    plot_count_distribution(counts=il17a_counts, cytokine='IL17A',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IL17A', 'PSO_distribution']))
    plot_count_distribution(counts=ifng_counts, cytokine='IFNG',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IFNG', 'LICHEN_distribution']))
    plot_count_distribution(counts=il13_counts, cytokine='IL13',
                            save_folder=save_folder, key="Bulk_RNAseq", title="_".join(['IL13', 'AD_distribution']))
    plot_count_distribution(counts=il17a_prp_counts, cytokine='IL17A', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IL17A', 'PRP_distribution']))
    plot_count_distribution(counts=ifng_prp_counts, cytokine='IFNG', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IFNG', 'PRP_distribution']))
    plot_count_distribution(counts=il13_prp_counts, cytokine='IL13', save_folder=save_folder, key="Bulk_RNAseq",
                            title="_".join(['IL13', 'PRP_distribution']))

    # plot log2(expr + 1) histogram
    fig_genes = plt.figure(facecolor='w', edgecolor='k', figsize=fig_size)
    fig_genes.subplots_adjust(bottom=0.025, left=0.025, top=0.975, right=0.975)
    ax = fig_genes.add_subplot(1, 1, 1)
    sns.distplot(np.log2(bulk_rnaseq.values + 1), kde=False, bins=50)
    ax.set_ylabel(r'log_{2}(counts + 1)', fontsize=xy_fontsize)
    ax.set_xlabel('Counts', fontsize=xy_fontsize)
    ax.tick_params(labelsize=xy_ticks)
    plt.savefig(os.path.join(save_folder, "Counts_Distribution.pdf"))
    plt.close()


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "Figure_2DF", str(today))
    os.makedirs(output_path, exist_ok=True)

    # input path
    input_path = os.path.join("..", "..", "..", "input", "bulk_RNAseq")

    # Read bulk-RNAseq count matrix
    bulk_data = pd.read_csv(os.path.join(input_path, "bulkRNA_countMat.txt"), sep='\t')
    # Read in metaData
    meta_data = pd.read_excel(os.path.join(input_path, "bulkRNA_metaData.xlsx"))

    main(save_folder=output_path, bulk_rnaseq=bulk_data, metadata=meta_data)
