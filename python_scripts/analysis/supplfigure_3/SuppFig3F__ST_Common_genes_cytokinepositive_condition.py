#!/usr/bin/env python
"""Find common up-regulated genes in the cytokine positive groups across IL-17A, IL-13, and IFN-g
    File name: SuppFig3F__ST_Common_genes_cytokinepositive_condition.py
    Author: Christina Hillig
    Date created: March/xx/2021
    Date last modified: April/30/2021
    Python Version: 3.7
"""

import pandas as pd
import os
from datetime import date


def get_intersection_genes(a, b, c):
    intersection_genes = set(a) & set(b) & set(c)

    return intersection_genes


def load_files(file):
    # Read out only those driver and responder genes which are specific for a cytokine
    df = pd.read_csv(file, error_bad_lines=False)
    # Remove column Unnamed: 0
    # df = df.drop(['Unnamed: 0'], axis=1)

    return df


def save_genes(genes, output_path, title):
    df = pd.DataFrame()
    df['gene_symbol'] = list(genes)
    df.to_csv(os.path.join(output_path, "".join([title, '.csv'])))
    df.to_excel(os.path.join(output_path, "".join([title, '.xlsx'])))


def main(save_folder):
    today = date.today()
    # create output path
    output_path = os.path.join(save_folder, "Pan_ncISDgenes", str(today))
    os.makedirs(output_path, exist_ok=True)

    input_path = os.path.join(save_folder, "output", "Figure_3B", str(today),
                              "spatial", "Whole_T_cell_matrix__cdr_patient_annotation_cyto")
    df_il17a = load_files(os.path.join(input_path, "IL17A", "IL17A_glmGamPoi_DGE.csv"))
    df_il13 = load_files(os.path.join(input_path, "IL13", "IL13_glmGamPoi_DGE.csv"))
    df_ifng = load_files(os.path.join(input_path, "IFNG", "IFNG_glmGamPoi_DGE.csv"))

    # get significant in cyto+ group
    m_sig_il17a = (df_il17a['pval'].values <= 0.05) & (df_il17a['log2fc'].values <= -1)
    m_sig_il13 = (df_il13['pval'].values <= 0.05) & (df_il13['log2fc'].values <= -1)
    m_sig_ifng = (df_ifng['pval'].values <= 0.05) & (df_ifng['log2fc'].values <= -1)
    pan_genes = get_intersection_genes(df_il17a['gene_symbol'][m_sig_il17a], df_il13['gene_symbol'][m_sig_il13],
                                       df_ifng['gene_symbol'][m_sig_ifng])

    print(len(list(pan_genes)))

    # Save Pan-inflamamtory skin disease genes
    save_genes(genes=pan_genes, output_path=output_path,
               title='Pan-inflammatory_skin_disease_genes__cdr_patient_annotation_cyto')

    save_genes(genes=df_il17a['gene_symbol'][m_sig_il17a], output_path=output_path,
               title='IL17A_SignificantGenes__cdr_patient_annotation_cyto')
    save_genes(genes=df_il13['gene_symbol'][m_sig_il13], output_path=output_path,
               title='IL13_SignificantGene__cdr_patient_annotation_cyto')
    save_genes(genes=df_ifng['gene_symbol'][m_sig_ifng], output_path=output_path,
               title='IFNG_SignificantGene__cdr_patient_annotation_cyto')


if __name__ == '__main__':
    savepath = os.path.join("..", "..", "..", "output")
    main(save_folder=savepath)


# todo path and copy input

