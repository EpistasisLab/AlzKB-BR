# AlzKB-BR: Extending the Alzheimer's Knowledge Base with Brain Region Transcriptomics for Interpretable Drug Discovery

## Overview

AlzKB-BR extends the Alzheimer's Knowledge Base (AlzKB) by integrating brain region-specific transcriptomic information from the Allen Human Brain Atlas (ABA). The goal is to improve the biological interpretation of knowledge graph-based drug repurposing predictions for Alzheimer's disease (AD) by incorporating neuroanatomical context.

<p align="center">
  <img src="figures/study overview.png" width="1000">
</p>

---

## Study Workflow

### 1. Brain Region Data Integration

* Obtain transcriptomic data from the Allen Human Brain Atlas (ABA)
* Perform one-versus-all differential expression analysis using **limma**
* Generate gene–brain region relationships
* Extend AlzKB to create **AlzKB-BR**

### 2. Permuted AlzKB-BR Negative Control KG

* We imported a negative control triples into Memgraph through the process from https://github.com/EpistasisLab/AlzKB/blob/master/BUILD.org
  1) Set up Memgraph locally with Docker and import the AlzKB Cypher export
  2) Remove the original brain-region–gene relationships
  3) Copy the permuted triples CSV into the Memgraph container
  4) Load the permuted triples and recreate the overexpression and underexpression relationships


### 2. Knowledge Graph Evaluation

* Compare **AlzKB**,  **AlzKB-BR**, and **Permuted AlzKB-BR Negative Control KG**
* Train the **RotatE** knowledge graph embedding model
* Evaluate performance using **MRR**, **Hits@1**, **Hits@3**, and **Hits@10**

### 3. Drug Prediction and Prioritization

* Predict missing **Drug–Treats–Disease** relationships
* Aggregate predictions using **Robust Rank Aggregation (RRA)**
* Apply statistical and biological filtering

### 4. Biological Interpretation

* Interpret candidate drugs using AD-associated brain regions
* Analyze associated gene sets using **GeneAgent**
* Review supporting clinical trial and literature evidence

---

## Repository Structure

```text
.
├── data/
│   ├── Allen_Human_Brain_Atlas/
│   ├── AlzKB/
│   │   └── CYPHERL/
│   └── AlzKB-BR/
│       └── CYPHERL/
│
├── output/
│   ├── DEG_results/
│   ├── Validation_results/
│   └── Drug_selection_results/
│
├── scripts/
│   ├── data_acquisition_curation.R
│   ├── probe_to_ensembl_mapping.R
│   ├── differential_expression_analysis.R
│   ├── calculate_z_scores.R
│   ├── validate_rotate_model.py
│   ├── robust_rank_aggregation.R
│   ├── construct_brain_region_relationships.ipynb
│   └── biological_refinement.ipynb,permuted_AlzKB-BR_negative-control_KG.py

│
├── figures/
├── supplementary/
└── README.md
```

---

## Methods

### Brain Region Integration

* Allen Human Brain Atlas (ABA) transcriptomic data
* Probe-to-gene mapping
* Differential expression analysis using **limma**
* BrainRegion ontology construction
* Integration into **AlzKB-BR**

### Knowledge Graph Completion

* **PyKEEN**
* **RotatE** knowledge graph embedding
* **Optuna** hyperparameter optimization
* 10 independent random seeds

### Drug Prioritization

* Robust Rank Aggregation (RRA)
* Median + Median Absolute Deviation (MAD) filtering
* Biological refinement

  * FDA-approved drugs
  * AD-associated genes
  * Pharmacogene exclusion

### Biological Interpretation

* GeneAgent (GPT-4 based)
* AD-associated brain regions
* Drug indications
* ClinicalTrials.gov
* Literature review

---

## Key Results

* Added **330 BrainRegion** nodes
* Added **25,174** brain region-related relationships
* Preserved RotatE knowledge graph completion performance
* Identified **36** prioritized AD drug candidates
* Identified **22** candidate drugs unique to AlzKB-BR
* Improved biological interpretation through neuroanatomical context

---

## Software

* Python
* R
* Neo4j
* PyKEEN
* limma
* Optuna
* GeneAgent
