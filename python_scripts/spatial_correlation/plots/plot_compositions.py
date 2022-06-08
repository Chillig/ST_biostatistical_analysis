import scanpy as sc
import numpy as np
import pandas as pd
import os
from operator import itemgetter
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.colors
import webcolors
from matplotlib.colors import ListedColormap
import seaborn as sns


def plot_stackedbars_perdiag(df_percentage_perdiag, responder_genes_cyto, cyto, save_folder):
    gene_color = dict(zip(responder_genes_cyto, sc.pl.palettes.default_102))
    for diag in df_percentage_perdiag.columns:
        fig, ax = plt.subplots(figsize=(8, 8))
        sorted_genes_temp = df_percentage_perdiag.sort_values(diag, ascending=False)
        df_percentage_perdiag.sort_values(diag, ascending=False).transpose().plot.bar(
            stacked=True, ax=ax, rot=0, color=itemgetter(*list(sorted_genes_temp.index))(gene_color))
        # remove grid
        ax.grid(False)
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Add y-axis label
        ax.set_ylabel('Gene frequency [%]')
        # put legend outside of axis & change the number of columns here & remove frame
        ax.legend(bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}sorted_{}_respondergene_composition.png'.format(diag, cyto)))
        plt.close(fig=fig)


def plot_bars_pergene(df_counts_responder_diag, cyto, save_folder):
    fig, ax = plt.subplots(figsize=(16, 8))
    sns.barplot(x="disease", y="value", data=df_counts_responder_diag, hue='variable',
                palette=sc.pl.palettes.default_102)
    # remove grid
    ax.grid(False)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Add y-axis label
    ax.set_ylabel('UMI-counts')
    # put legend outside of axis & change the number of columns here & remove frame
    ax.legend(bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Barplot_{}_respondergene_counts.png'.format(cyto)))
    plt.close(fig=fig)


def get_composition_percentage(df_counts, responder_genes_cyto, cyto):
    # diagnosis
    diagnosis = ['Pso', 'LP', 'AD']  # list(df_counts['disease'].unique())

    df_counts_cytorespspots = df_counts[~df_counts['{}_responder'.format(cyto)].isna()]
    # total number of counts per responder gene
    total_counts_pergene = df_counts_cytorespspots[responder_genes_cyto].sum(axis=0)

    # Get percentage and counts
    df_percentage_pergene = pd.DataFrame(columns=diagnosis, index=responder_genes_cyto)
    df_percentage_mean = pd.DataFrame(columns=diagnosis, index=responder_genes_cyto)
    df_percentage_perdiag = pd.DataFrame(columns=diagnosis, index=responder_genes_cyto)
    for diag in diagnosis:
        df_cyto_temp = df_counts_cytorespspots.groupby('disease')[responder_genes_cyto].get_group(diag)
        total_counts_perdiag = df_cyto_temp[responder_genes_cyto].sum(axis=0).sum()
        for resp_gene in responder_genes_cyto:
            # Calculating gene composition a disease
            df_percentage_perdiag.loc[resp_gene, diag] = np.divide(df_cyto_temp[resp_gene].sum(),
                                                                   total_counts_perdiag) * 100

            # Calculating Percentage of a gene over all diseases
            df_percentage_pergene.loc[resp_gene, diag] = np.divide(df_cyto_temp[resp_gene].sum(axis=0),
                                                                   total_counts_pergene[resp_gene]) * 100
            # Calculating Mean of a gene over all diseases
            df_percentage_mean.loc[resp_gene, diag] = df_cyto_temp[resp_gene].mean(axis=0)
            # # errors should be positive, and defined in the order of lower, upper
            # pso_errors = [[df_percentage_mean.loc[resp_gene, diag] - df_ifng_temp[resp_gene].min(),
            #                df_ifng_temp[resp_gene].max() - df_percentage_mean.loc[resp_gene, diag]] for c in df3.columns]

    return df_percentage_perdiag, df_percentage_pergene, df_percentage_mean


def get_gene_composition_percentage(df_counts, responder_genes_cyto, cyto):
    # diagnosis
    diagnosis = ['Pso', 'LP', 'AD']  # list(df_counts['disease'].unique())

    df_counts_cytorespspots = df_counts[~df_counts['{}_responder'.format(cyto)].isna()]
    # total number of counts per responder gene
    total_counts_pergene = df_counts_cytorespspots[responder_genes_cyto].sum(axis=0)


def main(df_counts: pd.DataFrame, cytokine_responders: dict, save_folder: str):
    # add sub-folder
    save_folder = os.path.join(save_folder, 'plot_compositions')
    os.makedirs(save_folder, exist_ok=True)
    # Cut-off in %, assuming that each gene is equally expressed
    cutoff = dict()
    for cyto in cytokine_responders.keys():
        cutoff[cyto] = (100 / len(cytokine_responders[cyto]))

    # 1.1 map each string to an integer value
    color_seq = pd.factorize(df_counts['disease'])[0]
    num_unique_colors = len(np.unique(color_seq))
    # 2.2 assign a color to each combination
    dict_disease_color = OrderedDict(zip(df_counts['disease'], color_seq))
    diseasecomb_colors = plt.cm.get_cmap("tab20", num_unique_colors)

    hfont = {'fontname': 'Calibri'}  # main font

    for cyto in cytokine_responders.keys():
        responder_genes_cyto = cytokine_responders[cyto]  # 36
        responder_genes_cyto = np.unique(responder_genes_cyto)  # 35

        df_respdiag_counts_temp = df_counts.melt(
            id_vars="disease", value_vars=responder_genes_cyto)
        df_respdiag_counts_temp['disease'] = df_respdiag_counts_temp['disease'].astype('category')
        df_respdiag_counts_temp['disease'].cat.reorder_categories(['Pso', 'LP', 'AD'], inplace=True)

        plot_bars_pergene(df_counts_responder_diag=df_respdiag_counts_temp, cyto=cyto, save_folder=save_folder)

        df_percentage_perdiag, df_percentage_pergene, df_percentage_mean = get_composition_percentage(
            df_counts=df_counts, responder_genes_cyto=responder_genes_cyto, cyto=cyto)
        plot_stackedbars_perdiag(df_percentage_perdiag=df_percentage_perdiag, responder_genes_cyto=responder_genes_cyto,
                                 save_folder=save_folder, cyto=cyto)

        # hue disease, x=gene, y=counts
        # fig, ax = plt.subplots(figsize=(16, 8))
        # sns.barplot(x="variable", y="value", data=df_respdiag_counts_temp, hue='disease',
        #             palette=sc.pl.palettes.default_102, ax=ax)
        # total count per gene per disease to create stacked barplot
        # hue disease, x=gene, y=percentage
        # test = df_respdiag_counts_temp.groupby(['variable', 'disease']).sum()['value'].value_counts(
        #     normalize=True).mul(100).rename('percent').reset_index()
        df_respdiag_counts_groupped_temp = df_respdiag_counts_temp.groupby(['variable', 'disease']).sum()
        test_percentage = df_respdiag_counts_groupped_temp.copy()
        for respgene_temp in responder_genes_cyto:
            total_count_temp = df_respdiag_counts_groupped_temp.loc[(respgene_temp)]['value'].sum()
            test_percentage.loc[(respgene_temp)]['value'] = np.divide(
                df_respdiag_counts_groupped_temp.loc[(respgene_temp)], total_count_temp)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax = test_percentage.unstack()['value'].plot(
            kind='bar', stacked=True, ax=ax, color=diseasecomb_colors.colors)
        # add text
        # .patches is everything inside of the chart
        for rect in ax.patches:
            # Find where everything is located
            height = rect.get_height()
            width = rect.get_width()
            x = rect.get_x()
            y = rect.get_y()

            # The height of the bar is the data value and can be used as the label
            height_temp = height * 100
            if height_temp < 0.01:
                label_text = "{:.2e}".format(height_temp)
            else:
                label_text = "{:.2f}".format(height_temp)  # decimal values

            # ax.text(x, y, text)
            label_x = x + width / 2
            label_y = y + height / 2

            # plot only when height is greater than specified value
            if height_temp > 0:
                ax.text(label_x, label_y, label_text, ha='center', va='center', fontsize=8)

        # remove grid
        ax.grid(False)
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Add y/x-axis label
        ax.set_ylabel('UMI-counts')
        ax.set_xlabel('')
        # put legend outside of axis & change the number of columns here & remove frame
        ax.legend(bbox_to_anchor=(1.0, 1.0), ncol=1, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}_disease_composition.png'.format(cyto)))
        plt.close(fig=fig)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax = df_respdiag_counts_groupped_temp.unstack()['value'].plot(
            kind='bar', stacked=True, ax=ax, color=diseasecomb_colors.colors)
        # .patches is everything inside of the chart
        for rect in ax.patches:
            # Find where everything is located
            height = rect.get_height()
            width = rect.get_width()
            x = rect.get_x()
            y = rect.get_y()

            # The height of the bar is the data value and can be used as the label
            height_temp = height * 100
            label_text = "{}".format(int(height_temp))  # decimal values

            # ax.text(x, y, text)
            label_x = x + width / 2
            label_y = y + height / 2

            # plot only when height is greater than specified value
            if height_temp > 0:
                ax.text(label_x, label_y, label_text, ha='center', va='center', fontsize=8)
        # for p in ax.patches:
        #     width, height = p.get_width(), p.get_height()
        #     x, y = p.get_xy()
        #     ax.text(x + width / 2,
        #             y + height / 2,
        #             '{:.0f}'.format(width),
        #             horizontalalignment='center',
        #             verticalalignment='center',
        #             color='white',
        #             fontsize=14,
        #             **hfont)
        # Scale y-axis
        ax.set_yscale("log")
        # remove grid
        ax.grid(False)
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Add y/x-axis label
        ax.set_ylabel('Gene frequency [%]')
        ax.set_xlabel('')
        # put legend outside of axis & change the number of columns here & remove frame
        ax.legend(bbox_to_anchor=(1.0, 1.0), ncol=1, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}_disease_composition_counts.png'.format(cyto)))
        plt.close(fig=fig)
