# calcium-binding-site-structural-clustering

# Structural Clustering of Calcium-Binding Sites Using Minimal Functional Sites

This repository contains the scripts and processed results used in the structural analysis of calcium-binding sites based on the **Minimal Functional Site (MFS)** representation.

The pipeline implements a two-stage clustering strategy to identify structurally recurrent calcium-binding environments across proteins belonging to different structural superfamilies.

The methodology combines structural alignment using **MM-align**, clustering based on **TM-score–derived distances**, and representative structure selection using **experimental resolution obtained from the RCSB PDB**.

---

# Pipeline Overview

The workflow consists of two main clustering stages.

## 1. Intra-superfamily clustering

Structures are grouped within each **CATH** and **SCOPe** superfamily.

Structural similarity is evaluated using TM-score values obtained from structural alignments. Clustering is performed using a similarity threshold corresponding to:

TM-score ≥ 0.5

This step groups structurally similar calcium-binding sites within homologous structural superfamilies.

Representative structures are then selected based on **experimental resolution**, retrieved from the RCSB PDB.

---

## 2. Global clustering of representative sites

All representative MFSs are combined into a single dataset and subjected to a second clustering stage.

This stage applies a stricter structural similarity threshold:

distance < 0.3  
(TM-score ≥ 0.7)

Clusters containing fewer than **10 members** are excluded in order to focus the analysis on structurally recurrent calcium-binding architectures.

---

# Scripts

## similarity_matrix_and_clustering.py

Performs the **first-stage intra-superfamily clustering**.

### Function

- Reads structural similarity matrices derived from MM-align
- Converts similarity values into distances
- Performs clustering using a **TM-score threshold of 0.5**

### Output files

CATH_StructuralClusters_TMscore_0.5.txt  
SCOPe_StructuralClusters_TMscore_0.5.txt

These files contain the structural clusters identified within each superfamily.

---

## APIresolutionsCATH.py/APIresolutionsSCOPe.py

Retrieves experimental resolution values from the **RCSB PDB REST API** for structures belonging to **CATH and SCOPe clusters**.

### Function

- Queries experimental resolution values from the PDB
- Identifies the structures with the best experimental resolution
- Selects representative structures for each cluster

### Output files

experimental_resolutionsCATH.txt  
experimental_resolutionsSCOPe.txt  
CATH_Representatives_SelectedByBestResolution_TMscore_0.5.txt
SCOPe_Representatives_SelectedByBestResolution_TMscore_0.5.txt

---

## global_clustering_tm07.py

Performs the **second-stage global clustering** using the representative structures.

### Function

- Reads a pairwise distance matrix derived from average TM-scores
- Applies clustering using a strict similarity threshold

distance < 0.3  
(TM-score ≥ 0.7)

Clusters with fewer than **10 members** are removed.

### Output files

SecondStage_Clustering_MembershipList_0.3.txt  
SecondStage_Clustering_MembershipList_0.5.txt

---

# Data Files

## Representative structures

Global_RepresentativeStructures.txt

Contains the list of **2960 representative Minimal Functional Sites (MFSs)** selected from CATH and SCOPe clusters.

---

## Intra-superfamily clustering results

CATH_StructuralClusters_TMscore_0.5.txt  
SCOPe_StructuralClusters_TMscore_0.5.txt

These files contain the structural clusters generated during the first stage of the pipeline.

---

## Selected representatives

CATH_Representatives_SelectedByBestResolution_TMscore_0.5.txt  
SCOPe_Representatives_SelectedByBestResolution_TMscore_0.5.txt

These files contain the representative MFSs selected based on:

- best experimental resolution
- structural similarity within clusters

---

## Global clustering results

SecondStage_Clustering_MembershipList_0.3.txt  
SecondStage_Clustering_MembershipList_0.5.txt

These files contain the final cluster membership lists obtained in the second-stage clustering.

---

## Structural statistics

cluster_statistics_final_11_clusters.xlsx

This file contains summary statistics of the **11 final clusters**, including cluster sizes and structural features.

---

## Coordination geometry distributions

Geometry_PercentageDistribution_11Clusters_0.3.txt  
Geometry_PercentageDistribution_25Clusters_0.5.txt

These files contain the percentage distribution of coordination geometries observed within clusters.

---

# Data Availability

The calcium-binding site structures used in this study are **not distributed in this repository**.

All Minimal Functional Sites (MFSs) analyzed in this work were retrieved from the **MetalPDB database**, a curated resource of metal-binding sites extracted from structures deposited in the Protein Data Bank (PDB).

The dataset can be obtained directly from:

https://metalpdb.cerm.unifi.it

Users interested in reproducing the analysis should download calcium-binding MFS structures from MetalPDB and generate the structural similarity matrices used in the clustering pipeline.

This repository provides the **analysis scripts and processed clustering results**, but not the original structural data.

---

# Requirements

Python 3.x

Required libraries:

pandas  
requests

Install with:

pip install pandas requests
