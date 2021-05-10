# Low numbers of cytokine transcripts drive inflammatory skin diseases by initiating amplification cascades in localized epidermal clusters
Schäbitz A*, Hillig C*, Farnoud A, Jargosch M, Scala E, Pilz AC, Bhalla N, Mubarak M, Thomas J, Stahle M, 
Biedermann T, Schmidt-Weber CB, Theis F, Garzorz-Stark N, Eyerich K*, Menden MP*, Eyerich S*§ <br>

You can find the preprint version of the manuscript here: <br>
[add URL](https://..??)


## Abstract
Chronic inflammatory diseases are characterized by tissue infiltration of abundant and polyclonal T cells. 
It remains challenging to distinguish disease-driving from bystander immune cells. Here, we used spatial transcriptomics 
to investigate 52,000 human cutaneous transcriptomes of non-communicable inflammatory skin diseases (ncISD), including 
psoriasis, atopic dermatitis, and lichen planus. <br>
Despite the expected T cell infiltration, per spatial spot we observed only 1-10 transcripts for T cell cytokines 
central for the pathogenic reaction. We confirmed these findings in independent patient cohorts using in 
situ hybridization, bulk and single-cell sequencing, and in vitro T cell stimulation assays. 
Expression of cytokine transcripts was limited to lesional skin and presented in a disease-specific pattern, 
indicating functional relevance. In fact, we identified responder signatures in the direct proximity of cytokines, 
and showed that single cytokine transcripts initiate amplification cascades of thousands of specific responder 
transcripts forming localized inflammatory microenvironments. <br>
Thus, within the abundant and polyclonal T cell infiltrates of ncISD, a few T cells drive disease by initiating an 
inflammatory amplification cascade in their local microenvironment. Specifically targeting these T cells could form the 
basis of causative and individualized therapies of ncISD. <br>

*key words*: Spatial transcriptomics, chronic inflammatory skin diseases, immunological disorders, 
immunopathogenesis, adaptive immunity <br>


## Lets get started!
The source code used to generate the results and figures are in the py_scripts and r_scripts folder. 
Further, the developed *conditional density clustering algorithm* can be found in folder 
python_scripts/spatial_correlation. <br> 
All preprocessing steps and analysis are run in python. Solely the differential gene expression (DGE) analysis and 
pathway enrichment analysis are run in R. <br> 
A processed example data is provided in [Add URL](https://..??) and by running the pipeline the generated results will 
be saved in the output folder. <br>
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
The results generated in this study have been created from raw count matrices and pre-processed data. <br>
An example with a subset of our data on how to run the pre-processing pipeline and our analysis can be found as a 
[Add URL - google colab notebook]() or [Add URL - R markdown]().


### Preprocessing
Before the samples are analysed, they have been preprocessed by applying the standard procedures 
such as quality control (QC), filtering, normalisation and batch correction. 
You can find the scripts here: [Source code - Preprocessing](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/pre_processing) <br>
The preprocessing is initiated by running [main_preprocessing.py](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/pre_processing/main_preprocessing.py):
```{python}
/path/to/conda_dir/py37_sc_rpy2_diffxpy/bin/python /path/to/Publication_analysis/python_scripts/pre_processing/main_preprocessing.py
``` 


### Analysis
Figures generated with the source code can be recreated by running 
[main_analysis.py - Analysis](https://github.com/Chillig/ST_biostatistical_analysis/blob/main/python_scripts/analysis/main_analysis.py):
```{python}
/path/to/conda_dir/py37_sc_rpy2_diffxpy/bin/python /path/to/Publication_analysis/python_scripts/analysis/main_analysis.py
```


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
We investigated the functional relevance of the few cytokine transcripts in lesional ncISD skin by looking at the 
correlation between cytokine-positive spots and their responder signatures in the epidermis. 
To do so we applied a density-based clustering customised for ST data and calculated the spatially weighted Pearson 
correlation between the UMI-counts of cytokines and the responder counts. 
[Source code - spatial weighted correlation](https://github.com/Chillig/ST_biostatistical_analysis/tree/main/python_scripts/spatial_correlation) 
and a [Add URL - tutorial]() is provided.
This analysis can be run by:
```{python}
/path/to/conda_dir/py37_sc_rpy2_diffxpy/bin/python /path/to/Publication_analysis/python_scripts/spatial_correlation/main.py
```

## References

