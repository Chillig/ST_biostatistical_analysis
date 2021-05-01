# Few specific cytokine transcripts drive inflammatory skin diseases by initiating amplification cascades in localized epidermal cluster
Schäbitz A*, Hillig C*, Farnoud A, Jargosch M, Pilz C, Bhalla N, Mubarak M, Stahle M, Biedermann T, Schmidt-Weber C, Theis F, Garzorz-Stark N, Eyerich K*, Menden M*, Eyerich S* <br>

You can find the preprint version of the manuscript here: <br>
[add URL](https://..??)


## Abstract
Chronic inflammatory diseases are characterized by tissue infiltration of T cells which are typically abundant and polyclonal. However, it remains challenging to distinguish disease-driving from bystander immune cells. To address this, we investigated more than 40,000 cutaneous transcriptomes of non-communicable inflammatory skin diseases (ncISD) at spatial resolution. Despite the expected T cell infiltration, only 1-10 transcripts per skin section were observed for T cell cytokines central for the pathogenic reaction, namely IL-17A for psoriasis, IL-13 for atopic dermatitis or IFNG for lichen planus. <br>
These findings were confirmed in independent patient cohorts leveraging *in situ* hybridisation, bulk and single cell sequencing data. Expression of cytokine transcripts was limited to lesional skin in a disease-specific pattern, indicating a functional relevance. In fact, using a newly developed algorithm that identifies responder signatures in the direct proximity of cytokines, single cytokine transcripts were shown to initiate amplification cascades of thousands of specific responder transcripts forming localized inflammatory micro-environments (e.g. IL-17A inducing NOS2, DEFB4A, CXCL8, IL-36G, or IL-19). <br>
Thus, within the abundant and polyclonal T cell infiltrates of ncISD, only very few T cells are disease-driving by initiating an inflammatory amplification cascade in their local microenvironment. Specific targeting of these T cells is a prerequisite for causative and individualised therapies of ncISD. <br>

key words: Spatial transcriptomics, chronic inflammatory skin diseases, cytokines , .. <br>

Short summary of key message <br>

add Workflow figure


## Getting started
The source code used to generate the results and figures are in the py_scripts and r_scripts folder. Fruther, the developed *conditional density clustering algorithm* can be found in folder scripts/spatial_correlation. <br> 
All preprocessing and analysis are run in python. Solely the differential gene expression (DGE) analysis and pathway enrichment analysis are run in R. <br> 
The processed data used in this study is provided in [Add URL](https://..??) and the results generated by the code are saved in the output folder. 
You can clone the repository by typing in the terminal the command: 
```{bash}
git clone https://github.com/Chillig/ST_biostatistical_analysis.git
```

1. Dependencies
In order to run the preprocessing and analysis in python you will have to create a conda environment from the python_conda_env.yml. Before doing so, you will have to manually set the prefix variable in the python_conda_env.yml file to your directory. Now, you can run the follwoing commands in the terminal:
```{bash}
# The following command will create an env with the name py37_sc_rpy2_diffxpy
conda env create -f python_conda_env.yml
# Activate the conda env with
conda activate py37_sc_rpy2_diffxpy
```
For the DGE and pathway enrichment analysis R4.0.3 Patched (2020-10-23 r79366) is required.<br>


## Preprocessing

## Clustering

## DGE Analysis

## Pathway Enrichemnt Analysis

## Conditional Density Clustering
