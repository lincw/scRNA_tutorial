# Cell-Cell Communication Analysis

Cell-cell communication analysis aims to infer how different cell types interact within a tissue by examining the expression of ligand and receptor pairs in single-cell RNA-seq data. This helps researchers understand signaling networks, tissue organization, immune responses, and disease mechanisms.

## Why is it important?
- Reveals how cells coordinate their functions in complex tissues
- Identifies key signaling pathways and cellular interactions in development, disease, and therapy
- Can uncover novel therapeutic targets or biomarkers

## Common Tools
- **CellPhoneDB**: Statistical framework for inferring cell-cell communication networks using known ligand-receptor interactions.
- **CellChat**: R package for quantitative inference and analysis of intercellular communication networks.
- **NicheNet**: Predicts ligand-target links between interacting cells, integrating gene expression with prior knowledge of signaling and gene regulatory networks.
- **SingleCellSignalR**: R package for analyzing intercellular communication from scRNA-seq data.

## Typical Workflow
1. Annotate cell types or clusters in your data
2. Use a communication analysis tool to infer ligand-receptor interactions between cell types
3. Visualize interaction networks and interpret biological significance
4. Validate key interactions with orthogonal methods (e.g., immunostaining, FACS, spatial transcriptomics)

## Integrating Custom Protein-Protein Interaction Networks

For researchers interested in protein interaction networks, advanced analyses can incorporate external PPI databases with scRNA-seq data:

### 1. Using Custom Interaction Databases with CellPhoneDB

```python
# CellPhoneDB v3.0+ supports custom interaction databases
import cellphonedb

# Export your Seurat or Scanpy object to compatible format
# For Scanpy:
adata_meta = adata.obs[['cell_type']].reset_index()
adata_meta.columns = ['cell_id', 'cell_type']
adata_meta.to_csv('metadata.txt', sep='\t', index=False)

# Convert counts to compatible format
counts = adata.raw.X.T  # Transpose to genes x cells
var_names = adata.raw.var_names.values
obs_names = adata.obs_names.values
pd.DataFrame(counts, index=var_names, columns=obs_names).to_csv('counts.txt', sep='\t')

# Run with custom database
!cellphonedb database generate --user-interactions custom_interactions.csv --user-gene-input custom_genes.csv
!cellphonedb method statistical_analysis metadata.txt counts.txt --counts-data gene_name --database custom_interactions
```

### 2. Using NicheNet to Connect Ligands with Target Genes

```r
library(nichenetr)

# Load your protein interaction networks
# You can use built-in networks or custom ones
ligand_target_matrix = readRDS("ligand_target_matrix.rds")

# Identify potential ligands from sender cells
potential_ligands = sender_cells %>% 
  get_expressed_genes(min_percentage = 0.10) %>% 
  intersect(ligand_target_matrix %>% rownames())

# Identify target genes that respond to ligands
target_genes = receiver_cells %>% 
  get_expressed_genes() %>% 
  .[. %in% colnames(ligand_target_matrix)]

# Run NicheNet analysis
ligand_activities = predict_ligand_activities(target_genes, potential_ligands, ligand_target_matrix)

# Visualize results
ligand_target_heatmap(ligand_activities, ligand_target_matrix, top_ligands, target_genes)
```

### 3. Building Custom Network Analysis Pipeline

```python
# Generic workflow to integrate custom PPI networks
import networkx as nx
import numpy as np
import pandas as pd

# 1. Load your PPI database (e.g., STRING, BioGRID, IntAct)
ppi_df = pd.read_csv('ppi_database.csv')

# 2. Filter interactions to include only genes expressed in your scRNA-seq data
expressed_genes = set(adata.var_names)
ppi_df_filtered = ppi_df[
    ppi_df['protein1'].isin(expressed_genes) & 
    ppi_df['protein2'].isin(expressed_genes)]

# 3. Create a graph of interactions
G = nx.from_pandas_edgelist(ppi_df_filtered, 'protein1', 'protein2', edge_attr='score')

# 4. Identify cell-type specific subnetworks
def get_cell_type_network(adata, cell_type, G, min_expr=0.1):
    """Extract subnetwork of genes expressed in a specific cell type"""
    cells = adata.obs[adata.obs['cell_type'] == cell_type].index
    cell_data = adata[cells]
    
    # Find genes expressed in at least min_expr fraction of cells
    expr_mask = (np.sum(cell_data.X > 0, axis=0) / cell_data.n_obs) > min_expr
    expr_genes = cell_data.var_names[expr_mask]
    
    # Extract subgraph
    return nx.subgraph(G, expr_genes)

# 5. Analyze intercellular communication
def find_interactions(G_sender, G_receiver):
    """Find potential interactions between two cell type networks"""
    sender_genes = set(G_sender.nodes())
    receiver_genes = set(G_receiver.nodes())
    
    interactions = []
    for gene1 in sender_genes:
        for gene2 in receiver_genes:
            if G.has_edge(gene1, gene2):
                interactions.append((gene1, gene2, G[gene1][gene2]['score']))
    
    return pd.DataFrame(interactions, columns=['ligand', 'receptor', 'score'])
```

## References
- Efremova M, et al. (2020). CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. Nat Protoc, 15(4), 1484–1506. [https://www.nature.com/articles/s41596-020-0292-x](https://www.nature.com/articles/s41596-020-0292-x)
- Efremova M, et al. (2022). CellPhoneDB v3.0: new features facilitate dissection of cell-cell communication by enabling the importation of new databases. bioRxiv. [https://doi.org/10.1101/2022.06.20.496774](https://doi.org/10.1101/2022.06.20.496774)
- Jin S, et al. (2021). Inference and analysis of cell-cell communication using CellChat. Nat Commun, 12, 1088. [https://www.nature.com/articles/s41467-021-21246-9](https://www.nature.com/articles/s41467-021-21246-9)
- Browaeys R, et al. (2020). NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods, 17, 159–162. [https://www.nature.com/articles/s41592-019-0667-5](https://www.nature.com/articles/s41592-019-0667-5)
