# Trajectory & Pseudotime Analysis

Trajectory inference and pseudotime analysis allow reconstruction of continuous biological processes from snapshot scRNA-seq data. Instead of discrete clusters, cells are placed along one or more trajectories that represent differentiation paths, developmental stages, or response to stimuli.

## Understanding Trajectory Inference

- **Goal:** Model dynamic processes (e.g., development, cell fate decisions) by connecting cells in reduced-dimensional space to reflect lineage relationships.
- **Key concepts:**
  - **Lineage topology:** The branching structure of cell state transitions (trees, graphs).
  - **Branch points:** Moments where a progenitor population diverges into multiple fates.
  - **Gene expression dynamics:** How gene expression changes along trajectories.
- **Use cases:**
  - Reconstructing hematopoietic differentiation paths.
  - Mapping immune cell activation states over time.

## Understanding Pseudotime

- **Definition:** A relative measure that orders cells along a trajectory, reflecting their progress through a biological process. It does not correspond to actual time but to a continuum of cell states.
- **Applications:**
  - Identifying genes with dynamic expression patterns.
  - Comparing progression across conditions or treatments.
  - Inferring regulatory events (e.g., transcription factor activation).

## Common Tools
- **Monocle3** (R): Learns trajectories and orders cells using reversed graph embedding.
- **Slingshot** (R): Fits minimum spanning trees on reduced-dimensional embeddings to infer lineages.
- **PAGA** (Python/Scanpy): Uses graph abstraction to capture cluster connectivity and visualize topology.
- **scVelo** (Python): Estimates RNA velocity and latent time to predict future cell states and infer dynamic processes.

## Example Workflows

### R (Monocle3)
```r
library(monocle3)
# Convert Seurat to cell_data_set
cds <- as.cell_data_set(pbmc)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = 'pseudotime')
```

### Python (Scanpy + scVelo)
```python
import scanpy as sc
import scvelo as scv
# Preprocess and compute moments
sc.pp.filter_and_normalize(adata)
sc.pp.pca(adata)
scv.pp.moments(adata)
# RNA velocity and latent time
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='viridis')
```

## References
- Trapnell C, et al. (2014). The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. _Nat Biotechnol_, 32(4):381–386.
- Street K, et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. _BMC Genomics_, 19:477.
- Wolf FA, et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. _Genome Biol_, 20:59.
- Bergen V, et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. _Nat Biotechnol_, 38:1408–1414.
