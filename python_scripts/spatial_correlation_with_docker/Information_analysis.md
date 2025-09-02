# Spatial Transcriptomics Density Clustering and Correlation Analysis

This project analyses gene-of-interest (GOI) to responder relationships in spatial transcriptomics data using density-based clustering.

The pipeline:
1. Identifies GOI-positive clusters using a spatial density clustering algorithm.
2. Quantifies responder gene counts within and outside clusters.
3. Computes weighted Pearson and Spearman correlations with corrected p-values.
4. Evaluates cluster size and spatial resolution across different radii.
5. Generates plots highlighting tissue layers, biopsy types, diseases, and patients.


## Cluster Formation
Clusters are built step by step:
1. GOI positive spots serve as initial cluster seeds / biological anchors.
2. Nearest-neighbour GOI-positive spots are grouped.
3. For each anchor, a fixed radius search identifies responder-positive spots.
4. Responder spots within this radius are added to the cluster.
5. Clusters are merged if they share one or more overlapping spots.

This allows quantification of GOI–responder spatial interactions at different spatial scales.

---

## Table of Contents

- [Overview](#overview)
- [Outputs](#outputs)
- [Workflow Diagram](#workflow-diagram)

---

## Overview

The script _main.py_ allows you to:

1. Subset spatial transcriptomics data based on disease type and tissue layers.
2. Identify clusters of cytokine-positive spots using a **spatial density clustering** approach.
3. Compute counts of **responder genes** inside and outside clusters.
4. Calculate **weighted Pearson and Spearman correlations** between cytokine expression and responder genes.
5. Generate plots evaluating correlations, cluster sizes, and responder gene distributions.

The workflow is built around the core classes:
- **`SpatialDensityCluster`** 
- **`ClusterDataCollector`**

---

## Outputs
The script modul_density_clustering.py produces:
1. Cluster statistics for each radius and GOI.
2. Counts of GOI and responder genes.
3. Weighted correlations (Pearson/Spearman + p-values) between GOI and responder genes.
4. Plots:
   * Correlation vs. radius
   * Normalised responder counts across radii
   * GOI counts vs. responder counts within a cluster
5. Excel file summarizing all counts and correlations.
   * Cluster counts per radius 
   * Radius vs. correlation results 
   * Normalized responder counts


## Workflow Diagram
1. Raw Spatial Transcriptomics Data (.h5)
2. Subset by Disease & Tissue Layer
3. Identify GOI+ Spots (SpatialDensityCluster)
4. KDTree / cKDTree Clustering
5. Connected Components → Clusters of Spots
6. Count Responder Genes inside/outside clusters
7. Weighted Correlation (Pearson/Spearman + p-values)
8. Generate Plots & Export Excel
