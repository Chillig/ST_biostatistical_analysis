import seaborn as sns
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import date


def plot_cytokine_specificresponder(counts_dict, radius, cytokine_responders):
    today = date.today()
    os.makedirs(
        os.path.join('/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today)),
        exist_ok=True)
    for radi in radius:
        df_counts = counts_dict[radi]
        # for cyto in cytokine_responders.keys():
        # read out counts of IL17A, IFNg, and IL13
        df_counts = df_counts[np.any(~df_counts[['IFNG', 'IL13', 'IL17A']].isna(), axis=1)]
        # read out counts specific for IL17A
        for cyto in cytokine_responders.keys():
            # df[resps] = df_counts_cytos[cytokine_responders[resps]].sum(axis=1)
            df_counts_cytos = df_counts[cytokine_responders[cyto]]
            df = pd.DataFrame.from_dict({
                'Responder counts': df_counts_cytos.sum(axis=1).values,
                'IL17A': df_counts['IL17A'], 'IL13': df_counts['IL13'], 'IFNG': df_counts['IFNG'],
                'Disease': df_counts['disease']})

            df = pd.melt(df, value_vars=['IL17A', 'IL13', 'IFNG'], id_vars=['Responder counts', 'Disease'])
            # remove clusters where cytokine was not measured
            df = df[~df['value'].isna()]
            df['value'] = df['value'].astype(int)
            df['variable'] = df['variable'].astype('category')
            df['Disease'] = df['Disease'].astype('category')
            df['Disease'].cat.reorder_categories(['LP', 'AD', 'Pso'], inplace=True)

            df.to_csv(os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today),
                'r{}__Cytokine_vs_{}_Responders.csv'.format(radi, cyto)))

            # fig, ax = plt.subplots(figsize=(6, 6))
            # sns.swarmplot(x="value", y="variable", data=df, hue='Disease',
            #               palette=['#ff7f00', '#e41a1c', '#377eb8'], ax=ax)
            # ax.set_ylabel('Cytokines')
            # ax.set_xlabel('Responder counts'.format(cyto))
            # # Put a legend to the right of the current axis
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            # sns.despine(ax=ax)
            # plt.tight_layout()
            # plt.savefig(os.path.join(
            #     '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts/r{}__Cytokines_vs_{}_Responders__Disease.pdf'.format(
            #         radi, cyto)))
            # plt.close(fig=fig)

            fig, ax = plt.subplots(figsize=(6, 6))
            sns.swarmplot(x="Responder counts", y="variable", data=df, palette=['#ff7f00', '#e41a1c', '#377eb8'], ax=ax)
            ax.set_ylabel('Cytokine')
            ax.set_xlabel('{} Responder counts'.format(cyto))
            sns.despine(ax=ax)
            plt.tight_layout()
            plt.savefig(os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today),
                'r{}__Cytokines_vs_{}_Responders.pdf'.format(
                    radi, cyto)))
            plt.close(fig=fig)


def plot_cytokine_responders(counts_dict, radius, cytokine_responders):
    today = date.today()
    os.makedirs(
        os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today)),
        exist_ok=True)
    for radi in radius:
        df_counts = counts_dict[radi]
        for cyto in cytokine_responders.keys():
            # read out counts of cytokines, density clusters
            df_counts_cytos = df_counts[np.any(~df_counts[[cyto]].isna(), axis=1)]
            # read out counts specific for IL17A
            df = pd.DataFrame(columns=['IFNG', 'IL13', 'IL17A'])
            for resps in cytokine_responders.keys():
                df[resps] = df_counts_cytos[cytokine_responders[resps]].sum(axis=1)

            df = pd.melt(df, value_vars=['IL17A', 'IL13', 'IFNG'])

            df.to_csv(os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today),
                'r{}__{}_vs_Responders.csv'.format(radi, cyto)))

            fig, ax = plt.subplots(figsize=(6, 6))
            sns.swarmplot(x="value", y="variable", data=df, palette=['#ff7f00', '#e41a1c', '#377eb8'], ax=ax)
            ax.set_ylabel('Cytokine')
            ax.set_xlabel('Responder counts'.format(cyto))
            sns.despine(ax=ax)
            plt.tight_layout()
            plt.savefig(os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today),
                'r{}__{}_vs_Responders.pdf'.format(radi, cyto)))
            plt.close(fig=fig)


def plot_density_clusters_respondergenes(counts_dict, cytokine_responders):
    today = date.today()
    save_folder = os.path.join(
                '/Volumes/CH__data/ST_immune_publication/Revision/Fig5/Cytokines_vs_Responder_counts', str(today))
    os.makedirs(save_folder, exist_ok=True)

    # Read out infos for optimal radii
    df_counts_il17a = pd.DataFrame(counts_dict[0][np.any(~counts_dict[0][['IL17A']].isna(), axis=1)])
    # IFNG density cluster
    df_counts_ifng = pd.DataFrame(counts_dict[4][np.any(~counts_dict[4][['IFNG']].isna(), axis=1)])
    # IL13 density cluster
    df_counts_il13 = pd.DataFrame(counts_dict[3][np.any(~counts_dict[3][['IL13']].isna(), axis=1)])

    for cyto in cytokine_responders.keys():
        counts_il17a = list(df_counts_il17a[cytokine_responders[cyto]].sum(axis=1).values)
        counts_ifng = list(df_counts_ifng[cytokine_responders[cyto]].sum(axis=1).values)
        counts_il13 = list(df_counts_il13[cytokine_responders[cyto]].sum(axis=1).values)

        df_counts_il17a["{}_Responders".format(cyto)] = counts_il17a
        df_counts_ifng["{}_Responders".format(cyto)] = counts_ifng
        df_counts_il13["{}_Responders".format(cyto)] = counts_il13

        # normalise by number of spots in cluster
        df_counts_il17a["normed_{}_Responders".format(cyto)] = np.divide(
            counts_il17a, df_counts_il17a['Cluster_num_spots'])
        df_counts_ifng["normed_{}_Responders".format(cyto)] = np.divide(
            counts_ifng, df_counts_ifng['Cluster_num_spots'])
        df_counts_il13["normed_{}_Responders".format(cyto)] = np.divide(
            counts_il13, df_counts_il13['Cluster_num_spots'])

    # pip install statannotations
    # Plots
    # IL17A responder genes in INFG, IL13 and IL17A density clusters
    dict_il17a_responders = {'IFNG_cluster': df_counts_ifng['normed_IL17A_Responders'].values,
                             'IL13_cluster': df_counts_il13['normed_IL17A_Responders'].values,
                             'IL17A_cluster': df_counts_il17a['normed_IL17A_Responders'].values}
    df_il17a_responders = pd.DataFrame({key: pd.Series(value) for key, value in dict_il17a_responders.items()})
    df_il17a = pd.melt(frame=df_il17a_responders, value_vars=['IFNG_cluster', 'IL13_cluster', 'IL17A_cluster'])
    # drop NaNs
    df_il17a = df_il17a[~df_il17a['value'].isna()]
    #
    plot_boxplot(df=df_il17a, cyto='IL17A', save_folder=save_folder)

    # IFNG responder genes in INFG, IL13 and IL17A density clusters
    dict_ifng_responders = {'IFNG_cluster': df_counts_ifng['normed_IFNG_Responders'].values,
                            'IL13_cluster': df_counts_il13['normed_IFNG_Responders'].values,
                            'IL17A_cluster': df_counts_il17a['normed_IFNG_Responders'].values}
    df_ifng_responders = pd.DataFrame({key: pd.Series(value) for key, value in dict_ifng_responders.items()})
    df_ifng = pd.melt(frame=df_ifng_responders, value_vars=['IFNG_cluster', 'IL13_cluster', 'IL17A_cluster'])
    plot_boxplot(df=df_ifng, cyto='IFNG', save_folder=save_folder)

    # IL13 responder genes in INFG, IL13 and IL17A density clusters
    dict_il13_responders = {'IFNG_cluster': df_counts_ifng['normed_IL13_Responders'].values,
                            'IL13_cluster': df_counts_il13['normed_IL13_Responders'].values,
                            'IL17A_cluster': df_counts_il17a['normed_IL13_Responders'].values}
    df_il13_responders = pd.DataFrame({key: pd.Series(value) for key, value in dict_il13_responders.items()})
    df_il13 = pd.melt(frame=df_il13_responders, value_vars=['IFNG_cluster', 'IL13_cluster', 'IL17A_cluster'])
    plot_boxplot(df=df_il13, cyto='IL13', save_folder=save_folder)


def plot_boxplot(df, cyto, save_folder):
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(x="variable", y="value", data=df, palette=['#ff7f00', '#e41a1c', '#377eb8'], ax=ax)
    ax.set_ylabel('{} Responder genes'.format(cyto))
    ax.set_xlabel('Density clusters')
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, 'Density_clusters_vs_{}_Responders.pdf'.format(cyto)))
    plt.close(fig=fig)
