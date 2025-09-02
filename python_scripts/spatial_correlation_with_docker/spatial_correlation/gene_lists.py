"""Genes used in Publication and  more
    File name: gene_lists.py
    Author: Christina Hillig
    Date created: 8/17/2023
    Python Version: 3.8
"""
from collections import OrderedDict
import numpy as np


def get_publication_spatial_transcriptomics_cyto_resps(respone_type: str):
    """
    Get cytokines and their responder genes from the Kera Array data for figure 4A and 4B

    Source: Spatial transcriptomics landscape of lesions from non-communicable inflammatory skin diseases
    doi: https://doi.org/10.1038/s41467-022-35319-w

    :return: list, dict
    """
    cytoresps_dict = OrderedDict()

    # Type 1 LICHEN
    if respone_type == 'IFNG':
        cytoresps_dict["IFNG"] = list(
            np.unique(["UBD", "CXCL9", "CXCL10", "CXCL11", "IL32", "ICAM1", "BATF2", "GBP5", "GBP4", "CCL8", "IRF1",
                       "RSAD2", "MUC1", "CCL5", "IL4I1", "CCL4", "MEI1", "CD70", "SAMD9L", "CSF2",
                       "EBI3", "APOBEC3G", "TNFSF13B", "NR4A3", "TNFRSF9", "BATF3", "KLHDC7B", "ADAM19",
                       "SP140", "ZBP1", "IDO1", "IFI30", "MUC16", "LAYN", "BCL2L14"]))
    elif respone_type == 'IL13':
        # Type 2 AE - BET3L and FAM26D not in gene list
        cytoresps_dict["IL13"] = list(
            np.unique(["CCL26", "NTRK1", "HSD3B1", "BET3L", "SERPINB13", "CH25H",
                       "EML5", "P2RY1", "LOXL4", "FAM26D", "THBS1", "CCL2", "C1QTNF1", "CISH"]))
    elif respone_type == 'IL17':
        # Type 3 PSO
        cytoresps_dict["IL17A"] = list(
            np.unique(["SPRR2B", "LCN2", "IL19", "CXCL6", "SPRR2E", "IGF2", "SPRR2A", "SPRR2F", "SPRR2D", "DEFB4A", "PI3",
                       "S100A7A", "SERPINB3", "IL36G", "SERPINB4", "RHCG", "SLC6A14", "ZC3H12A", "DEFB103B",
                       "PRSS22", "PLAT"]))
    else:
        NotImplementedError()

    return list(cytoresps_dict.keys()), cytoresps_dict
