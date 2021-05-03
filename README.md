# Few specific cytokine transcripts drive inflammatory skin diseases by initiating amplification cascades in localized epidermal cluster
Schäbitz A*, Hillig C*, Farnoud A, Jargosch M, Scala E, Pilz C, Bhalla N, Mubarak M, Thomas J, Stahle M, Biedermann T, 
Schmidt-Weber CB, Theis F, Garzorz-Stark N, Eyerich K*, Menden M*, Eyerich S* <br>

You can find the preprint version of the manuscript here: <br>
[add URL](https://..??)


## Abstract
Chronic inflammatory diseases are characterized by tissue infiltration of T cells which are typically abundant and 
polyclonal. However, it remains challenging to distinguish disease-driving from bystander immune cells. 
To address this, we investigated more than 40,000 cutaneous transcriptomes of non-communicable inflammatory skin 
diseases (ncISD) at spatial resolution. Despite the expected T cell infiltration, only 1-10 transcripts per skin 
section were observed for T cell cytokines central for the pathogenic reaction, namely IL-17A for psoriasis, IL-13 for 
atopic dermatitis or IFNG for lichen planus. <br>
These findings were confirmed in independent patient cohorts leveraging *in situ* hybridisation, bulk and single cell 
sequencing data. Expression of cytokine transcripts was limited to lesional skin in a disease-specific pattern, 
indicating a functional relevance. In fact, using a newly developed algorithm that identifies responder signatures in 
the direct proximity of cytokines, single cytokine transcripts were shown to initiate amplification cascades of 
thousands of specific responder transcripts forming localized inflammatory micro-environments (e.g. IL-17A inducing
 NOS2, DEFB4A, CXCL8, IL-36G, or IL-19). <br>
Thus, within the abundant and polyclonal T cell infiltrates of ncISD, only very few T cells are disease-driving by 
initiating an inflammatory amplification cascade in their local microenvironment. Specific targeting of these T cells is 
a prerequisite for causative and individualised therapies of ncISD. <br>

*key words*: Spatial transcriptomics, chronic inflammatory skin diseases, cytokines, T-cells, 
leukocytes, density clustering, .. <br>

Short summary of key message <br>

add Workflow figure <br>
![Workflow_final.pdf](/Users/christina.hillig/Documents/Projects/IGSSE-TUM_Projects/K_Eyerich_skin_spatial_transcriptomics/Paper/Low_mRNA_cytokines/Figure_1/Workflow_final.pdf)


## Lets get started!
The source code used to generate the results and figures are in the py_scripts and r_scripts folder. 
Further, the developed *conditional density clustering algorithm* can be found in folder 
python_scripts/spatial_correlation. <br> 
All preprocessing steps and analysis are run in python. Solely the differential gene expression (DGE) analysis and 
pathway enrichment analysis are run in R. <br> 
The processed data used in this study is provided in [Add URL](https://..??) and the results generated by the code are 
saved in the output folder. <br>
You can clone the repository by typing in the terminal the command: 
```{bash}
git clone https://github.com/Chillig/ST_biostatistical_analysis.git
```

### Dependencies
In order to run the preprocessing and analysis in python you will have to create a conda environment from the 
python_conda_env.yml. Before doing so, you will have to manually set the prefix variable in the python_conda_env.yml 
file to your directory. Now, you can run the following commands in the terminal:
```{bash}
# The following command will create an env with the name py37_sc_rpy2_diffxpy
conda env create -f python_conda_env.yml
# Activate the conda env with
conda activate py37_sc_rpy2_diffxpy
```
For the DGE and pathway enrichment analysis the R version R 4.0.3 Patched (2020-10-23 r79366) is required. 
Additionally, the following cran and Bioconductor packages are needed: <br> 
```{r}
# cran packages
rlibs <- c("dplyr", "gtools", "hash", "kableExtra", "knitr", "stringr", "tibble", "xlsx")
invisible(lapply(rlibs, require, character.only = TRUE))
# Bioconductor packages (Bioconductor version 3.12 (BiocManager 1.30.12))
bioclibs <- c("glmGamPoi", "pathview",  "org.Hs.eg.db", "ReactomePA",  "enrichplot", "clusterProfiler")
invisible(lapply(bioclibs, require, character.only = TRUE))
```


## Tutorial
The results generated in this study have been created from un-preprocessed and pre-processed data. <br>
In a first step, the quality of all 64 samples has been checked and 56 samples passed the Quality Control (QC). 
These 56 samples were then used for the analysis.  
In a second step the data was filtered for low expressed genes and putative batch effects were identified.



### Preprocessing
Before the samples are analysed, they have been preprocessed by applying the standard procedures 
such as quality control (QC), filtering, normalisation and batch correction. 
You can find the scripts here: [Source code - Preprocessing](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/pre_processing) <br>
The preprocessing is initiated by running [main_preprocessing.py](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/pre_processing/main_preprocessing.py):
```{python}
/path/to/conda_dir/py37_sc_rpy2_diffxpy/bin/python /path/to/Publication_analysis/python_scripts/pre_processing/main_preprocessing.py
``` 
It starts with creating a config file [Config - init_variables.py](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/pre_processing/init_variables.py) 
and continues with loading the count matrix together with the images, manual annotations, spot positions, 
and further information. The input is saved in [adata_storage](https://github.com/Chillig/ST_biostatistical_analysis/blob/main/python_scripts/adata_storage). <br>
The config file contains information about, for instance, data set type, applied preprocessing steps, and path variables.
This gives you the possibility to check later which parameters and thresholds you have set. <br>
During the preprocessing phase, the script will ask you which thresholds for the QC and how many 
Principal Components (PCs) to use.<br>
#### Overview of used Threshold and algorithms
Applied Thresholds: <br>
QC: 
1. 
1. 
1. 
PCs: 8 n_pcs ST and 7 n_pcs for singel cell data
Batch correction algorithm: [scanorama]()


### Analysis
Figures generated with the source code can be recreated by running 
[main_analysis.py - Analysis](https://github.com/Chillig/ST_biostatistical_analysis/blob/main/python_scripts/analysis/main_analysis.py):
```{python}
/path/to/conda_dir/py37_sc_rpy2_diffxpy/bin/python /path/to/Publication_analysis/python_scripts/analysis/main_analysis.py
```
List Parameters: <br>
1. p-value cut-off: 0.05
1. log2FC cut: 1.0
1. adjusted p-value: 0.05

### DGE analysis
In order to determine characteristic genes associated with cytokine-expressing leukocytes, 
a DGE analysis was performed between spots / cells containing cytokine-positive leukocytes and 
cytokine-negative leukocytes. A Vignette can be found here 
[Vignette - DGE Analysis](https://github.com/Chillig/ST_biostatistical_analysis/blob/main/r_scripts/dge_analysis/Vignette__DGE_Analysis.Rmd).

### Pathway enrichment analysis
The output of the DGE analysis containing the gene names, p-values, and log2FC are used as input for 
this type of analysis. A Vignette can be found here 
[Vignette - Pathway Enrichment Analysis](https://github.com/Chillig/ST_biostatistical_analysis/blob/main/r_scripts/pathway_analysis/Vignette__Pathway_Enrichment_Analysis.Rmd).

### Correlation between cytokines and their signature responder genes

#### Conditional density clustering

#### Spatially weighted correlation 

## References

