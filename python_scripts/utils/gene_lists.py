"""Genes used in Publication
    File name: gene_lists.py
    Author: Christina Hillig
    Date created: 2/01/2021
    Date last modified: 3/21/2021
    Python Version: 3.7
"""
from collections import OrderedDict


def leukocyte_markers():
    """Get leukocytes positive spots by 'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD247' and 'PTPRC' surface markers

    Returns
    -------

    """
    marker = ['CD2', 'CD3D', 'CD3E', 'CD3G', 'CD247', 'PTPRC']

    return marker


def highlight_genes():
    """Driver and responder genes for the cytokines IL-17A, IFN-g, and IL-13

    Genes which shall be highlighted in Volcano plot and Boxplots for each cytokine

    Returns
    -------

    """

    highlighting_genes = OrderedDict()
    highlighting_genes['IL17A'] = OrderedDict()
    highlighting_genes['IFNG'] = OrderedDict()
    highlighting_genes['IL13'] = OrderedDict()

    # TODO check with Steffi if we should remove also PPARD
    highlighting_genes['IL17A']['Driver_genes'] = ["IL17F", 'IL26', 'IL22', 'CCL3', 'CD274', 'PPARD', 'RAB27A']
    highlighting_genes['IL17A']["Responder_genes"] = ['NOS2', 'IL19', 'CCL20', 'LCN2', 'SPPRR2F', 'SPPRR2D', 'SPPRR2B',
                                                      'DEFB4A', 'S100A7A', 'IL36G', 'VNN3', 'CXCL8', 'CXCL1']

    # removed KLRG1, VACM1
    highlighting_genes['IFNG']['Driver_genes'] = ['GZMB', "FASLG", 'CD70', 'CRTAM', "CXCR6", 'CXCR3']
    highlighting_genes['IFNG']["Responder_genes"] = ['CXCL13', 'CXCL10', 'CXCL11', 'CXCL9', 'FGF9']

    # removed CXCR4
    highlighting_genes['IL13']['Driver_genes'] = ['IL2', 'IL10', 'CD48']
    highlighting_genes['IL13']["Responder_genes"] = ['CCL17', 'CCL19', 'CCL26', 'OSM']

    return highlighting_genes


def get_publication_cyto_resps():
    """
    Get cytokines and their responder genes from the Kera Array data for figure 4A and 4B

    :return: list, list, dict
    """
    cytoresps_dict = OrderedDict()
    # Type 1 LICHEN
    cytoresps_dict["IFNG"] = ["UBD", "CXCL9", "CXCL10", "CXCL11", "IL32", "ICAM1",
                              "BATF2", "GBP5", "GBP4", "CCL8", "IRF1", "IRF1", "RSAD2", "MUC1"]
    # Type 2 AE
    cytoresps_dict["IL13"] = ["CCL26", "NTRK1", "HSD3B1", "BET3L", "SERPINB13", "CH25H",
                              "EML5", "P2RY1", "LOXL4", "FAM26D", "THBS1"]
    # Type 3 PSO
    cytoresps_dict["IL17A"] = ["SPRR2B", "LCN2", "IL19", "CXCL6", "SPRR2E",
                               "IGF2", "SPRR2A", "SPRR2F", "SPRR2D", "DEFB4A"]

    cytokines = list(cytoresps_dict.keys())
    allinone = ["IL17A", "SPRR2B", "LCN2", "IL19", "CXCL6", "SPRR2E", "IGF2",
                "SPRR2A", "SPRR2F", "SPRR2D", "DEFB4A",
                "IFNG", "UBD", "CXCL9", "CXCL10", "CXCL11", "IL32", "ICAM1",
                "BATF2", "GBP5", "GBP4", "CCL8", "IRF1", "IRF1", "RSAD2", "MUC1",
                "IL13", "CCL26", "NTRK1", "HSD3B1", "BET3L", "SERPINB13", "CH25H",
                "EML5", "P2RY1", "LOXL4", "FAM26D", "THBS1"]

    return cytokines, allinone, cytoresps_dict


def get_exclusive_cd4cd8():
    celltypes = OrderedDict()
    celltypes['CD4+'] = ["CD4"]
    celltypes['CD8+'] = ['CD8A', 'CD8B']  # can only occur as dimere: NO

    return celltypes


def get_color_signaturegenes():
    # Group cytokines into Type 1-3 (Signature T-cells)
    signatures = OrderedDict()
    # publication
    signatures["IFNG"] = "#ff7f00"  # orange LICHEN
    signatures["IL13"] = "#e41a1c"  # red AE
    signatures["IL17A"] = "#377eb8"  # blue PSO
    signatures["HkG"] = '#4daf4a'  # green GAPDH

    return signatures


def cyto_asdict():
    cytokines = OrderedDict()
    cytokines["IFNG"] = ["IFNG"]  # LICHEN
    cytokines["IL13"] = ["IL13"]  # AE
    cytokines["IL17A"] = ["IL17A"]  # PSO

    return cytokines


def get_permuted_respondergenes():
    """Get cytokines and their responder genes from the human Keratinocyte experiment

    Returns
    -------
    cytokines : list
    cytoresps_dict : dict

    """
    cytoresps_dict = OrderedDict()
    # Type 1 LICHEN
    cytoresps_dict["IFNG"] = ["UBD", "CXCL9", "CXCL10", "CXCL11", "IL32", "ICAM1",
                              "BATF2", "GBP5", "GBP4", "CCL8", "IRF1", "IRF1", "RSAD2", "MUC1",
                              "CCL5", "IL4I1", "CCL4", "MEI1", "CD70", "SAMD9L", "CSF2",
                              "EBI3", "APOBEC3G", "TNFSF13B", "NR4A3", "TNFRSF9", "BATF3", "KLHDC7B", "ADAM19",
                              "SP140", "ZBP1", "IDO1", "IFI30", "MUC16", "LAYN", "BCL2L14"]

    # Type 2 AE
    cytoresps_dict["IL13"] = ["CCL26", "NTRK1", "HSD3B1", "BET3L", "SERPINB13", "CH25H",
                              "EML5", "P2RY1", "LOXL4", "FAM26D", "THBS1", "CCL2", "C1QTNF1"]

    # Type 3 PSO
    cytoresps_dict["IL17A"] = ["SPRR2B", "LCN2", "IL19", "CXCL6", "SPRR2E",
                               "IGF2", "SPRR2A", "SPRR2F", "SPRR2D", "DEFB4A", "PI3", "S100A7A", "SERPINB3",
                               "IL36G", "SERPINB4", "RHCG", "SLC6A14", "ZC3H12A", "DEFB103B", "PRSS22", "PLAT"]

    # permuted cytokines and responders
    cytoresps_dict["IFNG_IL13_responder"] = cytoresps_dict["IL13"]
    cytoresps_dict["IFNG_IL17A_responder"] = cytoresps_dict["IL17A"]

    cytoresps_dict["IL13_IFNG_responder"] = cytoresps_dict["IFNG"]
    cytoresps_dict["IL13_IL17A_responder"] = cytoresps_dict["IL17A"]

    cytoresps_dict["IL17A_IFNG_responder"] = cytoresps_dict["IFNG"]
    cytoresps_dict["IL17A_IL13_responder"] = cytoresps_dict["IL13"]

    cytokines = list(cytoresps_dict.keys())

    return cytokines, cytoresps_dict
