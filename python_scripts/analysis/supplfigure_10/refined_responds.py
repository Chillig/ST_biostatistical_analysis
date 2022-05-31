import os
import pandas as pd
import numpy as np
import scanpy as sc
from datetime import date
from itertools import chain

import pylab as plt
from matplotlib_venn import venn3, venn3_circles, venn2

from python_scripts.utils import gene_lists
from python_scripts.spatial_correlation import density_clustering

import matplotlib.patches as patches

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)


def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1], fontsize=14, ha="center", va="center"):
    ax.text(
        x, y, text,
        horizontalalignment=ha,
        verticalalignment=va,
        fontsize=fontsize,
        color="black")


def get_labels(data, fill):
    """
    get a dict of labels for groups in data
    @type data: list[Iterable]
    @type fill: list
    @rtype: dict[str, str]
    input
      data: data to get label for
      fill: ["number"|"logic"|"percent"]
    return
      labels: a dict of labels for different sets
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
    Out[12]:
    {'001': '0',
     '010': '5',
     '011': '0',
     '100': '3',
     '101': '2',
     '110': '2',
     '111': '3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                     # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        data_size = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)

    return labels


def venn5(labels, names, cytokine, save_folder, **options):
    """
    plots a 5-set Venn diagram
    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)
    input
      labels: a label dict where keys are identified via binary codes ('00001', '00010', '00100', ...),
              hence a valid set could look like: {'00001': 'text 1', '00010': 'text 2', '00100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi, fontsize
    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(5)])
    figsize = options.get('figsize', (13, 13))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.428, 0.449, 0.87, 0.50, 155.0, colors[0])
    draw_ellipse(fig, ax, 0.469, 0.543, 0.87, 0.50, 82.0, colors[1])
    draw_ellipse(fig, ax, 0.558, 0.523, 0.87, 0.50, 10.0, colors[2])
    draw_ellipse(fig, ax, 0.578, 0.432, 0.87, 0.50, 118.0, colors[3])
    draw_ellipse(fig, ax, 0.489, 0.383, 0.87, 0.50, 46.0, colors[4])
    draw_text(fig, ax, 0.27, 0.11, labels.get('00001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.72, 0.11, labels.get('00010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.55, 0.13, labels.get('00011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.91, 0.58, labels.get('00100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.78, 0.64, labels.get('00101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.84, 0.41, labels.get('00110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.76, 0.55, labels.get('00111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.90, labels.get('01000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.39, 0.15, labels.get('01001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.42, 0.78, labels.get('01010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.50, 0.15, labels.get('01011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.67, 0.76, labels.get('01100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.70, 0.71, labels.get('01101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.74, labels.get('01110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.64, 0.67, labels.get('01111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.10, 0.61, labels.get('10000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.20, 0.31, labels.get('10001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.76, 0.25, labels.get('10010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.65, 0.23, labels.get('10011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.18, 0.50, labels.get('10100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.21, 0.37, labels.get('10101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.81, 0.37, labels.get('10110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.74, 0.40, labels.get('10111', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.27, 0.70, labels.get('11000', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.34, 0.25, labels.get('11001', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.33, 0.72, labels.get('11010', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.22, labels.get('11011', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.25, 0.58, labels.get('11100', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.28, 0.39, labels.get('11101', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.36, 0.66, labels.get('11110', ''), fontsize=fontsize)
    draw_text(fig, ax, 0.51, 0.47, labels.get('11111', ''), fontsize=fontsize)

    # legend
    draw_text(fig, ax, 0.02, 0.72, names[0], colors[0], fontsize=fontsize, ha="right")
    draw_text(fig, ax, 0.72, 0.94, names[1], colors[1], fontsize=fontsize, va="bottom")
    draw_text(fig, ax, 0.97, 0.74, names[2], colors[2], fontsize=fontsize, ha="left")
    draw_text(fig, ax, 0.88, 0.05, names[3], colors[3], fontsize=fontsize, ha="left")
    draw_text(fig, ax, 0.12, 0.05, names[4], colors[4], fontsize=fontsize, ha="right")
    leg = ax.legend(names, loc='center left', bbox_to_anchor=(1.0, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'R1-5__VennDiaggramm_{}.pdf'.format(cytokine)))
    plt.close(fig=fig)


def plot_venndiagram(r1_genes, r2_genes, r3_genes, cyto, save_folder):
    set1 = set(r1_genes)
    set2 = set(r2_genes)
    set3 = set(r3_genes)

    # Change one group only
    fig, ax = plt.subplots()
    venn3([set1, set2, set3], ('Radius 1', 'Radius 2', 'Radius 3'), ax=ax)
    c = venn3_circles([set1, set2, set3],
                      linestyle='dashed', linewidth=1, color="grey", alpha=0.5, ax=ax)
    c[2].set_lw(8.0)
    c[2].set_ls('dotted')
    c[2].set_color('skyblue')

    fig.savefig(os.path.join(save_folder, '{}_VennDiagram.png'.format(cyto)))
    plt.close(fig=fig)


def plot_venn2(set1, set2, save_folder, cyto):
    # Second way
    fig, ax = plt.subplots(figsize=(6, 6))
    venn2(subsets=[set1, set2], set_labels=['Experimental', 'Computational'],
          set_colors=('darkorange', 'r'), ax=ax)
    plt.tight_layout()
    fig.savefig(os.path.join(save_folder, 'Keras_vs_Datadriven__{}_VennDiagram.pdf'.format(cyto)))
    plt.close(fig=fig)


def main():
    # Load adata
    project_folder = os.path.join("..", "..", "..")
    adata_folder = os.path.join(project_folder, "adata_storage")
    date_st_unpp = '2022-04-08'  # "2020-10-06" -> "st_adata_P15509_P16357_wo_4_7_unpp.h5"
    unpp_st_adata = sc.read(os.path.join(adata_folder, date_st_unpp, "Spatial Transcriptomics_unpp_cleaned.h5"))
    unpp_st_adata = unpp_st_adata[unpp_st_adata.obs['DISEASE'] != 'PRP'].copy()

    input_dir = '/Volumes/CH__data/ST_immune_publication/Revision/data/DGE_cyto_vs_others_Spearman'

    # TODO adjust cut-offs -> make more stringent than log2fc > 1
    log2fc_cutoff = 1.
    padj_cutoff = 0.001

    # create empty dict
    dict_newresponders = dict()
    # Per cytokine
    for cytokine in ['IL17A', 'IFNG', 'IL13']:
        save_folder = os.path.join(
            os.path.join(
                '/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers',
                'refined_response', str(date.today()), "log2FC_{}".format(log2fc_cutoff)))
        os.makedirs(save_folder, exist_ok=True)

        # Load in for radius 1-3 DEGs derived from refine_responder_genes.py and density_clustering
        df_r1 = pd.read_csv(os.path.join(input_dir, '1', '{}_in_sdcc_wilcoxon__radius1.csv'.format(cytokine)))
        df_r2 = pd.read_csv(os.path.join(input_dir, '2', '{}_in_sdcc_wilcoxon__radius2.csv'.format(cytokine)))
        df_r3 = pd.read_csv(os.path.join(input_dir, '3', '{}_in_sdcc_wilcoxon__radius3.csv'.format(cytokine)))
        df_r4 = pd.read_csv(os.path.join(input_dir, '4', '{}_in_sdcc_wilcoxon__radius4.csv'.format(cytokine)))
        df_r5 = pd.read_csv(os.path.join(input_dir, '5', '{}_in_sdcc_wilcoxon__radius5.csv'.format(cytokine)))

        # Find up-regulated genes in cyto+
        r1_up_genes = df_r1.loc[(df_r1['1_logfoldchanges'] > log2fc_cutoff) & (df_r1['1_pvals_adj'] < padj_cutoff),
                                '1_names'].values
        r2_up_genes = df_r2.loc[(df_r1['1_logfoldchanges'] > log2fc_cutoff) & (df_r2['1_pvals_adj'] < padj_cutoff),
                                '1_names'].values
        r3_up_genes = df_r3.loc[(df_r3['1_logfoldchanges'] > log2fc_cutoff) & (df_r3['1_pvals_adj'] < padj_cutoff),
                                '1_names'].values
        r4_up_genes = df_r4.loc[(df_r4['1_logfoldchanges'] > log2fc_cutoff) & (df_r4['1_pvals_adj'] < padj_cutoff),
                                '1_names'].values
        r5_up_genes = df_r5.loc[(df_r5['1_logfoldchanges'] > log2fc_cutoff) & (df_r5['1_pvals_adj'] < padj_cutoff),
                                '1_names'].values

        # drop signature cytokine
        r1_up_genes.tolist().remove(cytokine)
        r2_up_genes.tolist().remove(cytokine)
        r3_up_genes.tolist().remove(cytokine)
        r4_up_genes.tolist().remove(cytokine)
        r5_up_genes.tolist().remove(cytokine)

        subset_12 = list(np.intersect1d(r1_up_genes, r2_up_genes))
        intersect = list(np.intersect1d(subset_12, r3_up_genes))
        intersect = list(np.intersect1d(intersect, r4_up_genes))
        intersect = list(np.intersect1d(intersect, r5_up_genes))

        # plot_venndiagram(r1_genes=r1_up_genes, r2_genes=r2_up_genes, r3_genes=r3_up_genes, r4_genes=r4_up_genes,
        #                  r5_genes=r5_up_genes, cyto=cytokine, save_folder=save_folder)
        # https://github.com/tctianchi/pyvenn
        labels = get_labels([r1_up_genes, r2_up_genes, r3_up_genes, r4_up_genes, r5_up_genes], fill=['number'])
        venn5(labels, names=['Radius 1', 'Radius 2', 'Radius 3', 'Radius 4', 'Radius 5'],
              cytokine=cytokine, save_folder=save_folder)

        # compare to responder genes
        _, _, cyto_resps = gene_lists.get_publication_cyto_resps()

        print('Intersection with old responder genes:')
        print(cytokine, np.intersect1d(cyto_resps[cytokine], intersect))
        print('=========================================\n')

        dict_newresponders[cytokine] = intersect

        # TODO add Venndiagarmm responder new vs keratonicyte experiment
        plot_venn2(set1=set(cyto_resps[cytokine]), set2=set(intersect), save_folder=save_folder, cyto=cytokine)

    # try density clustering with new responder genes..
    # parameter
    save_folder = os.path.join(
        os.path.join(
            '/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/output/reviewers',
            'refined_response', str(date.today()), "log2FC_{}".format(log2fc_cutoff)))
    os.makedirs(save_folder, exist_ok=True)
    tissue_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    epidermis_layers = ['upper EPIDERMIS', 'middle EPIDERMIS', 'basal EPIDERMIS']
    radius = list(np.arange(0, 10, 1))  # 1 or list [1, 2, 3, ..]
    density_clustering.main(
        adata=unpp_st_adata, save_folder=save_folder, tissue_types=tissue_layers, epidermis_layers=epidermis_layers,
        radii=radius, corr_method='spearman', get_plots=False, conditional_genes=list(dict_newresponders.keys()),
        conditionalgenes_responders=dict_newresponders)

    # TODO Run Fig 4A-C with new responder genes


if __name__ == '__main__':
    main()


