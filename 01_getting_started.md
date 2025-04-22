# Getting Started

## Prerequisites
- Install Conda (Anaconda or Miniconda).
- Basic familiarity with R (>=4.1) and Python (>=3.8).

## Environment Setup

We provide a Conda environment file to install required packages:
```bash
conda env create -f environment.yml
conda activate scrna_env
```

Alternatively, install packages manually:

### Python (pip)
```bash
pip install scanpy anndata matplotlib seaborn scvelo cellrank python-igraph louvain
```

### R (in R console)
```r
install.packages("Seurat")
BiocManager::install(c("scater","scran","monocle3","SingleCellExperiment","CellChat"))
```

## Example Dataset

We will use the 10x Genomics PBMC 3K dataset as an example. Download:
```bash
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_3k/pbmc_3k_filtered_feature_bc_matrix.h5 -P data/
```
