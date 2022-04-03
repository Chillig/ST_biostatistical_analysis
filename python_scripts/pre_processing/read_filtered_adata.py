import scanpy as sc
import pandas as pd
from python_scripts.utils import helper_tools as ht, add_observables

if __name__ == '__main__':
    adata = sc.read(
        '/Users/christina.hillig/Documents/Projects/annData_objects/spatial/2021-03-19/adata_P15509_P16357_QC.h5')

    adata.to_df().T.to_csv('/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/input/P15509_P16357_counts_QC.csv')

    adata = add_observables.add_tissue_obs(adata=adata)
    adata.obs['tissue_type'].to_csv('/Users/christina.hillig/PycharmProjects/ST_Immune_publication/Publication_analysis/input/Tissue_layers.csv')