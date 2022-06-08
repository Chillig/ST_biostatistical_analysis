from python_scripts.utils import add_observables

import scanpy as sc
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from datetime import date

# Load ST object with mapped cell types from reference data set
# Load mapping object
# Plot H&E image with Piecharts on cyto+ spots:
# https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_piecharts.html
# create barcharts
# Which cell types are the most abundant in cyto+ spots?


def merge_labels(df_map):
    col_names = list(np.arange(0, df_map.shape[1], 1))
    # df_map_celltypes = pd.DataFrame(
    #     index=['DC', 'Differentiated_KC', 'Fibroblasts', 'ILC1_3', 'ILC2', 'Inf_mac', 'LC', 'LE', 'Macro',
    #            'Mast_cell', 'Melanocyte', 'MigDC', 'Mono_mac', 'NK', 'Pericyte', 'Plasma', 'Proliferating_KC',
    #            'Schwann', 'Tc', 'Tc_IL13_IL22', 'Tc17_Th17', 'Th', 'Treg', 'Undifferentiated_KC', 'VE', 'moDC'],
    #     columns=col_names)
    # Show only immune cells and name other cell types others
    # df_map_celltypes = pd.DataFrame(
    #     index=['DC', 'Fibroblasts', 'ILC1_3', 'ILC2', 'Inf_mac', 'LC', 'LE', 'Macro', 'MigDC', 'Mono_mac', 'NK',
    #            'Pericyte', 'Schwann', 'Tc', 'Tc_IL13_IL22', 'Tc17_Th17', 'Th', 'Treg',  'moDC', 'Others'],
    #     columns=col_names)
    df_map_celltypes = pd.DataFrame(
        index=['APC', 'Innate cells', 'Epithelial cells', 'T cells', 'Endothelial cells', 'B cells', 'Schwann'],
        columns=col_names)
    df_map_celltypes['colors'] = sc.pl.palettes.zeileis_28[:len(df_map_celltypes.index)]
    for col in col_names:
        # df_map_celltypes.loc['DC', col] = df_map.loc[df_map.index.str.contains('DC[0-9]'), col].sum()
        # df_map_celltypes.loc['Differentiated_KC', col] = df_map.loc[
        #     df_map.index.str.contains('Differentiated_KC'), col].sum()
        # df_map_celltypes.loc['Fibroblasts', col] = df_map.loc[df_map.index.str.contains('F[0-9]'), col].sum()
        # df_map_celltypes.loc['ILC1_3', col] = df_map.loc[df_map.index.str.contains('ILC1_3'), col].sum()
        # df_map_celltypes.loc['ILC2', col] = df_map.loc[df_map.index.str.contains('ILC2'), col].sum()
        # df_map_celltypes.loc['Inf_mac', col] = df_map.loc[df_map.index.str.contains('Inf_mac'), col].sum()
        # df_map_celltypes.loc['LC', col] = df_map.loc[df_map.index.str.contains('LC_[0-9]'), col].sum()
        # df_map_celltypes.loc['LE', col] = df_map.loc[df_map.index.str.contains('LE[0-9]'), col].sum()
        # df_map_celltypes.loc['Macro', col] = df_map.loc[df_map.index.str.contains('Macro'), col].sum()
        # df_map_celltypes.loc['Mast_cell', col] = df_map.loc[df_map.index.str.contains('Mast_cell'), col].sum()
        # df_map_celltypes.loc['Melanocyte', col] = df_map.loc[df_map.index.str.contains('Melanocyte'), col].sum()
        # df_map_celltypes.loc['MigDC', col] = df_map.loc[df_map.index.str.contains('MigDC'), col].sum()
        # df_map_celltypes.loc['Mono_mac', col] = df_map.loc[df_map.index.str.contains('Mono_mac'), col].sum()
        # df_map_celltypes.loc['NK', col] = df_map.loc[df_map.index.str.contains('NK'), col].sum()
        # df_map_celltypes.loc['Pericyte', col] = df_map.loc[df_map.index.str.contains('Pericyte'), col].sum()
        # df_map_celltypes.loc['Plasma', col] = df_map.loc[df_map.index.str.contains('Plasma'), col].sum()
        # df_map_celltypes.loc['Proliferating_KC', col] = df_map.loc[
        #     df_map.index.str.contains('Proliferating_KC'), col].sum()
        # df_map_celltypes.loc['Schwann', col] = df_map.loc[df_map.index.str.contains('Schwann'), col].sum()
        # df_map_celltypes.loc['Tc', col] = df_map.loc[df_map.index.str.contains('Tc(?!_IL)'), col].sum()
        # df_map_celltypes.loc['Tc17_Th17', col] = df_map.loc[df_map.index.str.contains('Tc17_Th17'), col].sum()
        # df_map_celltypes.loc['Tc_IL13_IL22', col] = df_map.loc[
        #     df_map.index.str.contains('Tc_IL13_IL22'), col].sum()
        # df_map_celltypes.loc['Th', col] = df_map.loc[df_map.index.str.contains('Th(?!_)'), col].sum()
        # df_map_celltypes.loc['Treg', col] = df_map.loc[df_map.index.str.contains('Treg'), col].sum()
        # df_map_celltypes.loc['Undifferentiated_KC', col] = df_map.loc[
        #     df_map.index.str.contains('Undifferentiated_KC'), col].sum()
        # df_map_celltypes.loc['VE', col] = df_map.loc[df_map.index.str.contains('VE[0-9]'), col].sum()
        # df_map_celltypes.loc['moDC', col] = df_map.loc[df_map.index.str.contains('moDC'), col].sum()

        # df_map_celltypes.loc['Others', col] = df_map.loc[
        #     df_map.index.str.contains('KC|Mast_cell|Melanocyte|VE|Plasma'), col].sum()
        # df_map_celltypes.loc['DC', col] = df_map.loc[df_map.index.str.contains('DC[0-9]'), col].sum()
        # df_map_celltypes.loc['Fibroblasts', col] = df_map.loc[df_map.index.str.contains('F[0-9]'), col].sum()
        # df_map_celltypes.loc['ILC1_3', col] = df_map.loc[df_map.index.str.contains('ILC1_3'), col].sum()
        # df_map_celltypes.loc['ILC2', col] = df_map.loc[df_map.index.str.contains('ILC2'), col].sum()
        # df_map_celltypes.loc['Inf_mac', col] = df_map.loc[df_map.index.str.contains('Inf_mac'), col].sum()
        # df_map_celltypes.loc['LC', col] = df_map.loc[df_map.index.str.contains('LC_[0-9]'), col].sum()
        # df_map_celltypes.loc['LE', col] = df_map.loc[df_map.index.str.contains('LE[0-9]'), col].sum()
        # df_map_celltypes.loc['Macro', col] = df_map.loc[df_map.index.str.contains('Macro'), col].sum()
        # df_map_celltypes.loc['MigDC', col] = df_map.loc[df_map.index.str.contains('MigDC'), col].sum()
        # df_map_celltypes.loc['Mono_mac', col] = df_map.loc[df_map.index.str.contains('Mono_mac'), col].sum()
        # df_map_celltypes.loc['NK', col] = df_map.loc[df_map.index.str.contains('NK'), col].sum()
        # df_map_celltypes.loc['Pericyte', col] = df_map.loc[df_map.index.str.contains('Pericyte'), col].sum()
        # df_map_celltypes.loc['Schwann', col] = df_map.loc[df_map.index.str.contains('Schwann'), col].sum()
        # df_map_celltypes.loc['Tc', col] = df_map.iloc[list(df_map.index).index("Tc"), col].sum()
        # df_map_celltypes.loc['Tc17_Th17', col] = df_map.loc[df_map.index.str.contains('Tc17_Th17'), col].sum()
        # df_map_celltypes.loc['Tc_IL13_IL22', col] = df_map.loc[
        #     df_map.index.str.contains('Tc_IL13_IL22'), col].sum()
        # df_map_celltypes.loc['Th', col] = df_map.iloc[list(df_map.index).index("Th"), col].sum()
        # df_map_celltypes.loc['Treg', col] = df_map.loc[df_map.index.str.contains('Treg'), col].sum()
        # df_map_celltypes.loc['moDC', col] = df_map.loc[df_map.index.str.contains('moDC'), col].sum()

        df_map_celltypes.loc['APC', col] = df_map.loc[
            df_map.index.str.contains('DC[0-9]|Inf_mac|LC_[0-9]|Macro|Mast_cell|MigDC|Mono_mac|moDC'), col].sum()
        df_map_celltypes.loc['Innate cells', col] = df_map.loc[
            df_map.index.str.contains('ILC1_3|ILC2|NK'), col].sum()
        df_map_celltypes.loc['Epithelial cells', col] = df_map.loc[
            df_map.index.str.contains('F[0-9]|Melanocyte|KC'), col].sum()
        df_map_celltypes.loc['T cells', col] = df_map.loc[
            df_map.index.str.contains('Tc|TcIL13_IL22|Tc17_Th17|Th|Treg'), col].sum()
        df_map_celltypes.loc['Endothelial cells', col] = df_map.loc[
            df_map.index.str.contains('LE[0-9]|Pericyte|VE'), col].sum()
        df_map_celltypes.loc['B cells', col] = df_map.loc[
            df_map.index.str.contains('Plasma'), col].sum()
        df_map_celltypes.loc['Schwann', col] = df_map.loc[
            df_map.index.str.contains('Schwann'), col].sum()

    return df_map_celltypes


def get_cropped_sampleimg(img_key):
    samples = ['P15509_1004', 'P16357_1003',  'P16357_1020']
    if img_key == 'lowres':
        crops_img = [(49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 530, 120)]
    else:
        res_proportion = 2000 / 600
        crops_img = [(49, 349, 510, 100), (80, 380, 510, 100), (150, 450, 530, 120)]
        crops_img = (np.asarray(crops_img) * res_proportion).round()

    return samples, crops_img


def plot_barhgraph(df):
    df_bar_celltypes = df.copy()
    # df_bar = df_bar.reset_index()
    df_bar_celltypes = df_bar_celltypes.drop('colors', axis=1)
    df_bar_celltypes = df_bar_celltypes.transpose()

    x = df_bar_celltypes.index
    indexes = np.argsort(df_bar_celltypes.values).T
    widths = np.sort(df_bar_celltypes.values).T
    order = -1
    lefts = widths[::order].cumsum(axis=0)
    lefts = np.insert(lefts, 0, np.zeros(len(lefts[0])), axis=0)

    mpp_colors = dict(zip(df_bar_celltypes.columns, sc.pl.palettes.zeileis_28))

    fig, ax = plt.subplots(figsize=(12, 6))
    for k, (idxs, vals) in enumerate(list(zip(indexes, widths))[::order]):
        mps = np.take(np.array(df_bar_celltypes.columns), idxs)
        # ax.bar(x, height=vals, color=[mpp_colors[m] for m in mps])
        ax.barh(x, width=vals, left=lefts[k], color=[mpp_colors[m] for m in mps])

    ax.legend((np.take(np.array(df_bar_celltypes.columns), np.argsort(df_bar_celltypes.values)[0]))[::-1],
              bbox_to_anchor=(1.05, 1), loc='upper left', ncol=1)


def plot_unsorted_bargraph(df, df_map, patches, sample_name, save_folder):
    df_bar = df.copy()
    df_bar = df_bar.drop('colors', axis=1)
    df_bar = df_bar.transpose()
    df_bar['point'] = list(df_bar.index)

    df_bar = df_bar.sort_values(by='T cells', ascending=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    df_bar.plot(
        x='point', kind='bar', stacked=True, title=sample_name, mark_right=True, rot=0,
        color=sc.pl.palettes.zeileis_28[:df_map.shape[0]], ax=ax)
    ax.set_ylabel('Cell type composition [%]')
    ax.set_xlabel('Spot')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    plt.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    plt.tight_layout()

    plt.savefig(os.path.join(save_folder, '{}_pointcomposition_barplot.pdf'.format(sample_name)))
    plt.close(fig=fig)


def plot_sorted_byoccurence_bargraph(df_bar_celltypes, sample_name, df_colors, patches, save_folder):
    fig, ax = plt.subplots(figsize=(12, 6))
    df_bar_celltypes.plot(
        x='point', kind='bar', stacked=True, title=sample_name, mark_right=True, rot=0,
        color=df_colors.iloc[0, :].values, ax=ax)
    ax.set_ylabel('Cell type composition [%]')
    ax.set_xlabel('Spot')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    plt.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, '{}_sorted_by_mostoccuringcelltype_barplot.pdf'.format(sample_name)))
    plt.close(fig=fig)


def plot_he_celltype_piecharts(st_adata, st_adata_cyto, sample_name, ind_sample, df_map_celltypes, save_folder):
    scale_factor = st_adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']
    samples, crops_img = get_cropped_sampleimg(img_key='hires')
    fig, ax = plt.subplots(figsize=(8, 10))
    sc.pl.spatial(st_adata, color='tissue_layer', alpha=0.0, ax=ax, s=1,
                  show=False, crop_coord=crops_img[ind_sample])
    # -- for loop: piecharts
    for ind_drawpie in range(0, st_adata_cyto.shape[0], 1):
        ax.pie(df_map_celltypes[ind_drawpie].values.astype(float),
               colors=df_map_celltypes['colors'],
               center=(int(st_adata_cyto.obsm['spatial'][ind_drawpie][0] * scale_factor),
                       int(st_adata_cyto.obsm['spatial'][ind_drawpie][1] * scale_factor)),
               radius=10)
        # Add point number to spots / pichart
        ax.annotate(ind_drawpie, (int(st_adata_cyto.obsm['spatial'][ind_drawpie][0] * scale_factor) + 10,
                                  int(st_adata_cyto.obsm['spatial'][ind_drawpie][1] * scale_factor) + 10),
                    fontsize=5, rotation=180)
    ax.set_xlim(crops_img[ind_sample][0:2])
    ax.set_ylim(crops_img[ind_sample][2:4])

    patches = []
    for ind_label, val in enumerate(list(df_map_celltypes.index)):
        patches.append(mpatches.Patch(color=sc.pl.palettes.zeileis_28[ind_label], label=val))

    plt.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    plt.tight_layout()
    plt.savefig(os.path.join(save_folder, '{}_celltype_H&E.pdf'.format(sample_name)))
    plt.close(fig=fig)

    return patches


def main(save_folder):
    # Labels
    labels = ['DC1', 'DC2', 'Differentiated_KC', 'Differentiated_KC*', 'F1', 'F2', 'F3', 'ILC1_3', 'ILC2', 'Inf_mac',
              'LC_1', 'LC_2', 'LC_3', 'LC_4', 'LE1', 'LE2', 'Macro_1', 'Macro_2', 'Mast_cell', 'Melanocyte', 'MigDC',
              'Mono_mac', 'NK', 'Pericyte_1', 'Pericyte_2', 'Plasma', 'Proliferating_KC', 'Schwann_1', 'Schwann_2',
              'Tc', 'Tc_IL13_IL22', 'Th', 'Treg', 'Undifferentiated_KC', 'VE1', 'VE2', 'VE3',
              'moDC_1', 'moDC_2', 'moDC_3']
    samples = ['P15509_1004', 'P16357_1003', 'P16357_1020']
    samples_folder = {'P15509_1004': 'S4', 'P16357_1003': 'S7', 'P16357_1020': 'S24'}
    cytokine = {'P15509_1004': 'il17a', 'P16357_1003': 'il13', 'P16357_1020': 'ifng'}
    input_dir = '/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/adata_storage/2022-05-13'

    for ind_sample, sample_name in enumerate(samples):
        # Mapping object - cell type probabilities are stored in .X
        adata_mapping = sc.read(
            os.path.join(input_dir, '{}'.format(samples_folder[sample_name]), 'SCST_tangrammapping.h5'))
        st_adata = sc.read(os.path.join(input_dir, '{}'.format(samples_folder[sample_name]), 'ST_tangram.h5'))

        # Add annotations from obsm to obs
        for annotation in list(adata_mapping.obs['subclass_label'].cat.categories):
            st_adata.obs[annotation] = st_adata.obsm['tangram_ct_pred'][annotation].copy()

        # Add info of cyto+ spots
        st_adata, obs_name = add_observables.convert_variable_to_observable(
            adata=st_adata, gene_names=['il17a', 'ifng', 'il13'], task='cell_gene', label='celltype', condition=None)

        # Keep only cyto+ spots
        cyto_columns = list(st_adata.obs.columns[
                                st_adata.obs.columns.str.contains('_{}'.format(cytokine[sample_name]))])
        mask = st_adata.obs[cyto_columns].values == cytokine[sample_name]
        st_adata_cyto = st_adata[mask].copy()
        st_mapping_cyto = adata_mapping[:, mask].copy()

        # Calculate percentage for each cell type
        max_value = st_mapping_cyto.X.sum(axis=0)
        st_mapping_cytoper = st_mapping_cyto.X / max_value
        df_map = pd.DataFrame.from_dict(data=st_mapping_cytoper)
        df_map.index = list(st_mapping_cyto.obs['subclass_label'].cat.categories)

        labels = list(adata_mapping.obs['subclass_label'].cat.categories)
        df_sample = pd.melt(st_adata_cyto.obs, id_vars=['spot_type'], value_vars=labels)

        # Overall composition
        result = df_sample.groupby(["variable"])['value'].aggregate(np.mean).reset_index().sort_values('value')
        result = result.sort_values('value', ascending=False)

        # Plot cell type composition per cyto+ spot
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(data=df_sample, x="variable", y="value", order=result['variable'], log=True, hue='spot_type', ax=ax)
        ax.set_xticklabels(result['variable'], rotation=90, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}__cytospots_composition.png'.format(sample_name)))
        plt.close(fig=fig)

        # Merge sub-cell types
        df_map_celltypes = merge_labels(df_map=df_map)

        # Plot H&E image + Pie charts
        patches = plot_he_celltype_piecharts(
            st_adata=st_adata, st_adata_cyto=st_adata_cyto, sample_name=sample_name, ind_sample=ind_sample,
            df_map_celltypes=df_map_celltypes, save_folder=save_folder)

        print('Plot Barplot')
        # index are the points, columns are the cell types
        plot_unsorted_bargraph(
            df=df_map_celltypes, df_map=df_map, patches=patches, sample_name=sample_name, save_folder=save_folder)

        # Sort by number of cell types within a spot; start by least occuring and end by most occuring
        df_bar_celltypes = df_map_celltypes.copy()
        df_bar_celltypes = df_bar_celltypes.drop('colors', axis=1)
        df_bar_celltypes = df_bar_celltypes.transpose()
        # Save colors
        colors = sc.pl.palettes.zeileis_28[:df_map_celltypes.shape[0]]
        df_colors = pd.DataFrame(dict(zip(list(df_map_celltypes.index), colors)), index=[0])
        # start by least occuring and end by most occuring
        # 1. Get ranking
        ranked_celltypes = list(df_bar_celltypes.sum(axis=0).sort_values(ascending=True).index)
        for col in ranked_celltypes:
            # reorder spots
            df_bar_celltypes = df_bar_celltypes.sort_values(by=col, ascending=False)
            # reorder columns
        df_bar_celltypes = df_bar_celltypes[df_bar_celltypes.sum(axis=0).sort_values(ascending=False).index]
        # add point
        df_bar_celltypes['point'] = list(df_bar_celltypes.index)

        # Sort colors
        df_colors = df_colors[df_bar_celltypes.columns[:-1]]
        plot_sorted_byoccurence_bargraph(
            df_bar_celltypes=df_bar_celltypes, sample_name=sample_name, df_colors=df_colors,
            patches=patches, save_folder=save_folder)

        # Save cell type composition per spot for a specimen
        df_bar_celltypes.to_excel(
            os.path.join(save_folder, 'Tangram_{}_sample{}.xlsx'.format(cytokine[sample_name], sample_name)))

        # TODO check composition in nn spots
        # TODO check composition in others
        # TODO check composition of cyto+ spots in epidermis and dermis


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "SupplFig8A", str(today))
    os.makedirs(output_path, exist_ok=True)

    main(save_folder=output_path)
