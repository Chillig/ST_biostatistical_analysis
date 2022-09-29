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


class DeconvolutionAnalysis:
    def __init__(self, save_folder, adata_st, adata_mapping, sample_name):
        self.save_folder = save_folder
        self._adata_st = adata_st
        self._adata_mapping = adata_mapping

        self.n_spots = self._adata_st.shape[0]

        # Add info of cyto+ spots
        self._adata_st, obs_name = add_observables.convert_variable_to_observable(
            adata=self._adata_st, gene_names=['il17a', 'ifng', 'il13'],
            task='cell_gene', label='celltype', condition=None)

        self.mask_celltype = None
        self.df_sample_celltype_probabilities = None
        self.df_map = None
        self.df_map_celltypes = None

        self.mergedlables_celltypes = {
            'APCs': 'DC[0-9]|Inf_mac|LC_[0-9]|Macro|MigDC|Mono_mac|moDC',
            'Lymphoid & Mast': 'ILC1_3|NK|ILC2|Tc|TcIL13_IL22|Tc17_Th17|Th|Treg|Mast_cell|Plasma',
            'Dermal non-immune': 'Schwann|F[0-9]|Pericyte|VE|LE[0-9]',
            'Epidermal non-immune': 'Melanocyte|KC'}

    @property
    def adata_st(self):
        # Add annotations from obsm to obs
        self._adata_st.obs = pd.concat([self._adata_st.obs, self._adata_st.obsm['tangram_ct_pred']], axis=1)

        # Add column which states if spot belongs to epidermis or Dermis
        self._adata_st.obs['Epi-Dermis'] = 'Unknown'
        self._adata_st.obs.loc[self._adata_st.obs['EPIDERMIS'] == 1, 'Epi-Dermis'] = 'Epidermis'
        self._adata_st.obs.loc[self._adata_st.obs['DERMIS'] == 1, 'Epi-Dermis'] = 'Dermis'

        return self._adata_st.copy()

    @property
    def adata_mapping(self):
        return self._adata_mapping.copy()

    def determine_celltype_in_spot(self):
        if self.mask_celltype is None:
            self.mask_celltype = self.adata_st.obs[
                                     self.adata_mapping.obs['subclass_label'].cat.categories] > (1 / self.n_spots)
            self.mask_celltype = self.mask_celltype.astype(int)

        return self.mask_celltype

    def sort_df_celltype_ranking(self, df):
        # Get ranking
        ranked_celltypes = list(df.sum(axis=0).sort_values(ascending=True).index)
        for col in ranked_celltypes:
            # reorder spots
            df = df.sort_values(by=col, ascending=False)
            # reorder columns
        df = df[df.sum(axis=0).sort_values(ascending=False).index]

        return df

    def get_probability_per_cell_in_spot(self):
        # Calculate percentage for each cell type
        adata_temp = self.adata_mapping.copy()
        # True: st_adata_cyto.obsm['tangram_ct_pred']['Differentiated_KC'].values == st_mapping_cyto.X[0]
        max_value = adata_temp.X.sum(axis=0)  # get max prob per spot
        # norm by maximum value to retriev percentage for each cell type to be in a spot
        st_mapping_cytoper = np.empty(shape=adata_temp.X.shape)
        for cell in range(0, len(max_value)):
            st_mapping_cytoper[:, cell] = (adata_temp.X[:, cell] - np.min(adata_temp.X[:, cell], axis=0)) / (
                    max_value[cell] - np.min(adata_temp.X[:, cell], axis=0))
        # store as dataframe
        self.df_map = pd.DataFrame.from_dict(data=st_mapping_cytoper)
        self.df_map.columns = np.arange(0, len(max_value))
        # add epidermis and dermis labels to spots
        df_epi_dermis = pd.DataFrame(self.adata_st.obs['Epi-Dermis'].values).transpose()
        df_epi_dermis.columns = np.arange(0, len(max_value))
        df_epi_dermis.index.names = ['Epi-Dermis']
        self.df_map = self.df_map.append(df_epi_dermis, ignore_index=True)
        # rename index
        self.df_map.index = list(adata_temp.obs['subclass_label'].values) + ['Epi-Dermis']

        # get labels
        labels = list(adata_temp.obs['subclass_label'].values)
        self.df_sample_celltype_probabilities = pd.melt(self.adata_st.obs, id_vars=['Epi-Dermis'], value_vars=labels)

        return self.df_sample_celltype_probabilities

    def get_merged_celltypes(self):
        df_temp = self.df_map[:-1]
        col_names = list(np.arange(0, df_temp.shape[1], 1))
        df_map_cellgroups_elltypes = pd.DataFrame(
            index=['APCs', 'Lymphoid & Mast', 'Dermal non-immune', 'Epidermal non-immune'],
            columns=col_names)
        df_map_cellgroups_elltypes['colors'] = sc.pl.palettes.zeileis_28[:len(df_map_cellgroups_elltypes.index)]

        for col in col_names:
            # Use grouping like in HaniffaLab paper
            df_map_cellgroups_elltypes.loc['APCs', col] = df_temp.loc[
                df_temp.index.str.contains(self.mergedlables_celltypes['APCs']), col].sum()
            df_map_cellgroups_elltypes.loc['Lymphoid & Mast', col] = df_temp.loc[
                df_temp.index.str.contains(
                    self.mergedlables_celltypes['Lymphoid & Mast']), col].sum()
            df_map_cellgroups_elltypes.loc['Dermal non-immune', col] = df_temp.loc[
                df_temp.index.str.contains(self.mergedlables_celltypes['Dermal non-immune']), col].sum()
            df_map_cellgroups_elltypes.loc['Epidermal non-immune', col] = df_temp.loc[
                df_temp.index.str.contains(self.mergedlables_celltypes['Epidermal non-immune']), col].sum()

        # add infos about tissue layer
        self.df_map_celltypes = df_map_cellgroups_elltypes.append(self.df_map.tail(1), ignore_index=False)

        return self.df_map_celltypes

    def get_subcelltypes(self, subcelltypes):
        df = self.df_map.loc[
            self.df_map.index.str.contains(self.mergedlables_celltypes[subcelltypes])]
        df = df.reindex(sorted(df.columns), axis=1)
        # Get colors
        colors = sc.pl.palettes.default_20[:df.shape[0]]
        df_colors = pd.DataFrame(dict(zip(list(df.index), colors)), index=[0])
        # norm probabilities to 0 and 1 for each spot
        df = np.divide(df, df.sum())
        # Transpose df
        df = df.transpose()
        # Get ranking
        df = self.sort_df_celltype_ranking(df=df)
        # Add spot number
        df['spot number'] = list(df.index)

        return df, df_colors

    def sort_subcelltypes_by_tissuelayer(self, df, subcelltypes, sample_name):
        df['Tissue_layer'] = self.df_map_celltypes.loc['Epi-Dermis'][:-1][list(df.index)]
        df.to_excel(
            os.path.join(self.save_folder, 'Tangram_sample{}__{}.xlsx'.format(sample_name, subcelltypes)))

        # Sort by Epidermis and Dermis
        df = df.sort_values('Tissue_layer', ascending=[False])
        df_epidermis = df[df['Tissue_layer'].str.contains('Epi')]
        df_dermis = df[df['Tissue_layer'].str.contains('D')]
        # Sort by most frequent cell type
        df_epidermis = self.sort_df_celltype_ranking(df=df_epidermis.iloc[:, :-2])
        df_dermis = self.sort_df_celltype_ranking(df=df_dermis.iloc[:, :-2])
        df = pd.concat([df_epidermis, df_dermis])
        # Add spot number
        df['spot number'] = list(df.index)
        df['Tissue_layer'] = self.df_map_celltypes.loc['Epi-Dermis'][:-1][list(df.index)]

        return df

    def draw_tissuelayer_piecharts(self, df_map_celltypes, index, value, sample_name):
        df_epidermis = df_map_celltypes[df_map_celltypes.columns[df_map_celltypes.loc[index] == value]][:-1]
        df_epidermis = df_epidermis.transpose()
        df_piechart = df_epidermis.sum() / df_epidermis.sum().sum()

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.pie(df_piechart.values.astype(float),
               colors=df_map_celltypes['colors'], labels=list(df_piechart.index),
               center=(0, 0), autopct='%.1f%%',
               radius=1)
        ax.set_title('Celltype frequency in the {}'.format(str(value).lower()))
        plt.tight_layout()
        fig.savefig(os.path.join(self.save_folder, '{}_{}_piechart.pdf'.format(sample_name, str(value).lower())))
        plt.close(fig=fig)

    def plot_he_celltype_piecharts(self, st_adata, sample_name, ind_sample, df_map_celltypes):
        scale_factor = st_adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']
        samples, crops_img = get_cropped_sampleimg(img_key='hires')
        fig, ax = plt.subplots(figsize=(8, 10))
        sc.pl.spatial(st_adata, color='tissue_layer', alpha=0.0, ax=ax, s=1,
                      show=False, crop_coord=crops_img[ind_sample])
        # -- for loop: piecharts
        for ind_drawpie in range(0, self.adata_st.shape[0], 1):
            ax.pie(df_map_celltypes[ind_drawpie].values.astype(float),
                   colors=df_map_celltypes['colors'],
                   center=(int(self.adata_st.obsm['spatial'][ind_drawpie][0] * scale_factor),
                           int(self.adata_st.obsm['spatial'][ind_drawpie][1] * scale_factor)),
                   radius=10)
            # Add point number to spots / piechart
            ax.annotate(ind_drawpie, (int(self.adata_st.obsm['spatial'][ind_drawpie][0] * scale_factor) + 10,
                                      int(self.adata_st.obsm['spatial'][ind_drawpie][1] * scale_factor) + 10),
                        fontsize=5, rotation=180)
        ax.set_xlim(crops_img[ind_sample][0:2])
        ax.set_ylim(crops_img[ind_sample][2:4])

        patches = []
        for ind_label, val in enumerate(list(df_map_celltypes.index)):
            patches.append(mpatches.Patch(color=sc.pl.palettes.zeileis_28[ind_label], label=val))

        plt.legend(handles=patches, loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
        plt.tight_layout()
        plt.savefig(os.path.join(self.save_folder, '{}_celltype_H&E.pdf'.format(sample_name)))
        plt.close(fig=fig)

        return patches

    def plot_barcharts(self, df, subcelltype, sample_name, df_colors, order):
        df_temp = df.reset_index()
        vline_position = df_temp[df_temp['Tissue_layer'] == 'Dermis'].index[0] - 0.5

        fig, ax = plt.subplots(figsize=(12, 6))
        df.plot.bar(
            x='spot number', stacked=True, title=sample_name, mark_right=True,
            color=df_colors.iloc[0, :].values, rot=0, ax=ax)
        ax.set_ylim([0, 1])
        ax.axvline(vline_position, c='black', ls='--')
        ax.set_ylabel('Cell type composition [%]')
        ax.set_xlabel('Spot')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, title=subcelltype)
        # get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()
        # add legend to plot with specified order of items in legend
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
                   loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, title=subcelltype)
        plt.tight_layout()
        plt.savefig(os.path.join(
            self.save_folder, '{}_{}_sorted_tissuelayer_barplot.pdf'.format(subcelltype, sample_name)))
        plt.savefig(os.path.join(
            self.save_folder, '{}_{}_sorted_tissuelayer_barplot.png'.format(subcelltype, sample_name)))
        plt.close(fig=fig)


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

        # # TODO norm by number of sub cell types included in a main cell type?
        # df_map_celltypes.loc['APC', col] = df_map.loc[
        #     df_map.index.str.contains('DC[0-9]|Inf_mac|LC_[0-9]|Macro|Mast_cell|MigDC|Mono_mac|moDC'), col].sum()
        # df_map_celltypes.loc['Innate cells', col] = df_map.loc[
        #     df_map.index.str.contains('ILC1_3|ILC2|NK'), col].sum()
        # df_map_celltypes.loc['Epithelial cells', col] = df_map.loc[
        #     df_map.index.str.contains('F[0-9]|Melanocyte|KC'), col].sum()
        # df_map_celltypes.loc['T cells', col] = df_map.loc[
        #     df_map.index.str.contains('Tc|TcIL13_IL22|Tc17_Th17|Th|Treg'), col].sum()
        # df_map_celltypes.loc['Endothelial cells', col] = df_map.loc[
        #     df_map.index.str.contains('LE[0-9]|Pericyte|VE'), col].sum()
        # df_map_celltypes.loc['B cells', col] = df_map.loc[
        #     df_map.index.str.contains('Plasma'), col].sum()
        # df_map_celltypes.loc['Schwann', col] = df_map.loc[
        #     df_map.index.str.contains('Schwann'), col].sum()

        df_map_celltypes.loc['APC', col] = df_map.loc[
            df_map.index.str.contains('APC'), col].sum()
        df_map_celltypes.loc['Innate cells', col] = df_map.loc[
            df_map.index.str.contains('Innate cells'), col].sum()
        df_map_celltypes.loc['Epithelial cells', col] = df_map.loc[
            df_map.index.str.contains('Epithelial cells'), col].sum()
        df_map_celltypes.loc['T cells', col] = df_map.loc[
            df_map.index.str.contains('T cells'), col].sum()
        df_map_celltypes.loc['Endothelial cells', col] = df_map.loc[
            df_map.index.str.contains('Endothelial cells'), col].sum()
        df_map_celltypes.loc['B cells', col] = df_map.loc[
            df_map.index.str.contains('B cells'), col].sum()
        df_map_celltypes.loc['Schwann', col] = df_map.loc[
            df_map.index.str.contains('Schwann'), col].sum()

    return df_map_celltypes


def get_cropped_sampleimg(img_key):
    samples = ['P15509_1001']
    if img_key == 'lowres':
        crops_img = [(49, 349, 580, 170)]
    else:
        res_proportion = 2000 / 600
        crops_img = [(49, 349, 480, 70)]
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

    df_bar = df_bar.sort_values(by='Lymphoid & Mast', ascending=False)

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


def plot_he_celltype_piecharts(st_adata, sample_name, ind_sample, df_map_celltypes, save_folder):
    scale_factor = st_adata.uns['spatial'][sample_name]['scalefactors']['tissue_hires_scalef']
    samples, crops_img = get_cropped_sampleimg(img_key='hires')
    fig, ax = plt.subplots(figsize=(8, 10))
    sc.pl.spatial(st_adata, color='tissue_layer', alpha=0.0, ax=ax, s=1,
                  show=False, crop_coord=crops_img[ind_sample])
    # -- for loop: piecharts
    for ind_drawpie in range(0, st_adata.shape[0], 1):
        ax.pie(df_map_celltypes[ind_drawpie].values.astype(float),
               colors=df_map_celltypes['colors'],
               center=(int(st_adata.obsm['spatial'][ind_drawpie][0] * scale_factor),
                       int(st_adata.obsm['spatial'][ind_drawpie][1] * scale_factor)),
               radius=10)
        # Add point number to spots / pichart
        ax.annotate(ind_drawpie, (int(st_adata.obsm['spatial'][ind_drawpie][0] * scale_factor) + 10,
                                  int(st_adata.obsm['spatial'][ind_drawpie][1] * scale_factor) + 10),
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


def sort_df_celltype_ranking(df):
    # Get ranking
    ranked_celltypes = list(df.sum(axis=0).sort_values(ascending=True).index)
    for col in ranked_celltypes:
        # reorder spots
        df = df.sort_values(by=col, ascending=False)
        # reorder columns
    df = df[df.sum(axis=0).sort_values(ascending=False).index]

    return df


def main(save_folder):
    # Labels
    # labels = ['DC1', 'DC2', 'Differentiated_KC', 'Differentiated_KC*', 'F1', 'F2', 'F3', 'ILC1_3', 'ILC2', 'Inf_mac',
    #           'LC_1', 'LC_2', 'LC_3', 'LC_4', 'LE1', 'LE2', 'Macro_1', 'Macro_2', 'Mast_cell', 'Melanocyte', 'MigDC',
    #           'Mono_mac', 'NK', 'Pericyte_1', 'Pericyte_2', 'Plasma', 'Proliferating_KC', 'Schwann_1', 'Schwann_2',
    #           'Tc', 'Tc_IL13_IL22', 'Th', 'Treg', 'Undifferentiated_KC', 'VE1', 'VE2', 'VE3',
    #           'moDC_1', 'moDC_2', 'moDC_3']
    samples = ['P15509_1001']
    samples_folder = {'P15509_1001': 'S1'}
    input_dir = '/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/adata_storage/2022-08-24'

    for ind_sample, sample_name in enumerate(samples):
        # Mapping object - cell type probabilities are stored in .X
        adata_mapping = sc.read(
            os.path.join(input_dir, '{}'.format(samples_folder[sample_name]), 'SCST_tangrammapping.h5'))
        st_adata = sc.read(os.path.join(input_dir, '{}'.format(samples_folder[sample_name]), 'ST_tangram.h5'))

        # Add annotations from obsm to obs
        for annotation in list(adata_mapping.obs['subclass_label'].cat.categories):
            st_adata.obs[annotation] = st_adata.obsm['tangram_ct_pred'][annotation].copy()

        df_test = st_adata.obsm['tangram_ct_pred']
        pd.DataFrame({'Cell type spot': df_test.columns[
            np.argmax(np.asarray(df_test), axis=1)]})['Cell type spot'].value_counts().hist()
        df_occurency = pd.DataFrame(pd.DataFrame({'Cell type spot': df_test.columns[np.argmax(np.asarray(df_test), axis=1)]})[
            'Cell type spot'].value_counts())
        df_occurency = df_occurency.reset_index()
        df_occurency.columns = ['Cell type', 'Frequency']
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(data=df_occurency, x=df_occurency['Cell type'],  y=df_occurency['Frequency'], ax=ax)
        ax.set_xticklabels(df_occurency['Cell type'].values, rotation=90, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}__Most_occuring_celltype_in_spot.png'.format(sample_name)))
        plt.close(fig=fig)

        deconv_obj = DeconvolutionAnalysis(save_folder=save_folder, adata_st=st_adata, adata_mapping=adata_mapping,
                                           sample_name=sample_name)

        df_sample = deconv_obj.get_probability_per_cell_in_spot()

        # Overall composition
        result = df_sample.groupby(["variable"])['value'].aggregate(np.mean).reset_index().sort_values('value')
        result = result.sort_values('value', ascending=False)

        # Plot cell type composition per cyto+ spot
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.barplot(
            data=df_sample, x="variable", y="value", order=result['variable'], log=True, hue='Epi-Dermis', ax=ax)
        ax.set_xticklabels(result['variable'], rotation=90, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(save_folder, '{}__cytospots_composition.png'.format(sample_name)))
        plt.close(fig=fig)

        # Merge sub-cell types
        df_map_celltypes = deconv_obj.get_merged_celltypes()

        # Draw epidermis and Dermis piecharts
        deconv_obj.draw_tissuelayer_piecharts(
            df_map_celltypes=df_map_celltypes, index='Epi-Dermis', value='Epidermis', sample_name=sample_name)
        deconv_obj.draw_tissuelayer_piecharts(
            df_map_celltypes=df_map_celltypes, index='Epi-Dermis', value='Dermis', sample_name=sample_name)

        # Get lymphoid cell probabilities
        df_lymphoid_mastcells, _ = deconv_obj.get_subcelltypes(subcelltypes='Lymphoid & Mast')
        # Group lymphoid Mast cells into groups
        df_lymphoid_mastcells_grouped = pd.DataFrame(
            columns=['ILC', 'T-cells', 'NK', 'Mast cells', 'Plasma'],
            index=df_lymphoid_mastcells.index)
        for col in df_lymphoid_mastcells_grouped.index:
            # Use grouping like in HaniffaLab paper
            df_lymphoid_mastcells_grouped.loc[col, 'ILC'] = df_lymphoid_mastcells.loc[
                col, df_lymphoid_mastcells.columns.str.contains('ILC')].sum()
            df_lymphoid_mastcells_grouped.loc[col, 'T-cells'] = df_lymphoid_mastcells.loc[
                col, df_lymphoid_mastcells.columns.str.contains('Tc|Tc17|Tc_IL13|Th|Treg')].sum()
            df_lymphoid_mastcells_grouped.loc[col, 'NK'] = df_lymphoid_mastcells.loc[
                col, df_lymphoid_mastcells.columns == 'NK'].sum()
            df_lymphoid_mastcells_grouped.loc[col, 'Mast cells'] = df_lymphoid_mastcells.loc[
                col, df_lymphoid_mastcells.columns.str.contains('Mast')].sum()
            df_lymphoid_mastcells_grouped.loc[col, 'Plasma'] = df_lymphoid_mastcells.loc[
                col, df_lymphoid_mastcells.columns.str.contains('Plasma')].sum()

        df_lymphoid_mastcells_grouped['spot number'] = df_lymphoid_mastcells['spot number']
        df_colors = pd.DataFrame(columns=['ILC', 'T-cells', 'NK', 'Mast cells', 'Plasma'],
                                 data=[['#1f77b4', '#ff7f0e', '#d62728', '#b5bd61', '#ffbb78']])

        # Save Lymphoid & Mast cell composition per spot for a specimen
        df_lymphoid_mastcells_grouped = deconv_obj.sort_subcelltypes_by_tissuelayer(
            df=df_lymphoid_mastcells_grouped, subcelltypes='Lymphoid & Mast', sample_name=sample_name)
        # sort columns of df_color by columns of other dataframe
        df_colors = df_colors[df_lymphoid_mastcells_grouped.columns[:-2]]
        deconv_obj.plot_barcharts(
            df=df_lymphoid_mastcells_grouped, subcelltype='Lymphoid & Mast', sample_name=sample_name,
            df_colors=df_colors, order=[list(df_lymphoid_mastcells_grouped.columns).index('T-cells'),
                                        list(df_lymphoid_mastcells_grouped.columns).index('NK'),
                                        list(df_lymphoid_mastcells_grouped.columns).index('ILC'),
                                        list(df_lymphoid_mastcells_grouped.columns).index('Mast cells'),
                                        list(df_lymphoid_mastcells_grouped.columns).index('Plasma')])
        # Save Groups to excel file
        df_lymphoid_mastcells_grouped.to_excel(os.path.join(save_folder, 'Lymphoid_Mastcells.xlsx'))

        # Get APCs probabilities
        df_apcs, df_colors = deconv_obj.get_subcelltypes(subcelltypes='APCs')
        df_apcs_grouped = pd.DataFrame(
            columns=['DC', 'Mo_DC', 'Macro', 'Mono_mac', 'LC'],
            index=df_apcs.index)
        for col in df_apcs_grouped.index:
            df_apcs_grouped.loc[col, 'DC'] = df_apcs.loc[
                col, df_apcs.columns.str.contains('DC[0-9]|MigDC')].sum()
            df_apcs_grouped.loc[col, 'LC'] = df_apcs.loc[
                col, df_apcs.columns.str.contains('LC_[0-9]')].sum()
            df_apcs_grouped.loc[col, 'Macro'] = df_apcs.loc[
                col, df_apcs.columns.str.contains('Macro_[1-2]|Inf_mac')].sum()
            df_apcs_grouped.loc[col, 'Mono_mac'] = df_apcs.loc[
                col, df_apcs.columns.str.contains('Mono_mac')].sum()
            df_apcs_grouped.loc[col, 'Mo_DC'] = df_apcs.loc[
                col, df_apcs.columns.str.contains('moDC_[0-9]')].sum()

        df_apcs_grouped['spot number'] = df_apcs['spot number']
        df_colors = pd.DataFrame(columns=['DC', 'LC', 'Macro', 'Mono_mac', 'Mo_DC'],
                                 data=[['#458B00', '#CD2626', '#8EE5EE', '#FF8C00', '#FF1493']])

        # Save APCs composition per spot for a specimen
        df_apcs_grouped = deconv_obj.sort_subcelltypes_by_tissuelayer(
            df=df_apcs_grouped, subcelltypes='APCs', sample_name=sample_name)
        # sort columns of df_color by columns of other dataframe
        df_colors = df_colors[df_apcs_grouped.columns[:-2]]
        deconv_obj.plot_barcharts(
            df=df_apcs_grouped, subcelltype='APCs', sample_name=sample_name, df_colors=df_colors,
            order=[0, 3, 1, 4, 2])
        # Save Groups to excel file
        df_apcs_grouped.to_excel(os.path.join(save_folder, 'APCs.xlsx'))

        # Plot H&E image + Pie charts
        patches = plot_he_celltype_piecharts(
            st_adata=st_adata, sample_name=sample_name, ind_sample=ind_sample,
            df_map_celltypes=df_map_celltypes[:-1], save_folder=save_folder)

        print('Plot Barplot')
        # index are the points, columns are the cell types
        plot_unsorted_bargraph(
            df=df_map_celltypes[:-1], df_map=deconv_obj.df_map, patches=patches, sample_name=sample_name, save_folder=save_folder)

        # Sort by number of cell types within a spot; start by least occuring and end by most occuring
        df_bar_celltypes = df_map_celltypes[:-1].copy()
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
            os.path.join(save_folder, 'Tangram_NL_sample{}.xlsx'.format(sample_name)))

        # TODO check composition in nn spots
        # TODO check composition in others
        # TODO check composition of cyto+ spots in epidermis and dermis


if __name__ == '__main__':
    today = date.today()
    # create saving folder
    output_path = os.path.join("..", "..", "..", "output", "SupplFig3GH", str(today))
    os.makedirs(output_path, exist_ok=True)

    main(save_folder=output_path)
