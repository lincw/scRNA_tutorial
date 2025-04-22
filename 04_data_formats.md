# Data Formats

Understanding file formats is essential for handling and converting single-cell RNA-seq data. Common formats include:

## FASTQ
- Raw sequencing reads with quality scores.
- Contains read sequences and per-base quality information.

## BAM/SAM
- Aligned reads mapped to a reference genome.
- Binary (BAM) or text (SAM) format including alignment coordinates and flags.

## Count Matrices
- **Matrix Market (MTX):** Sparse matrix format with separate files for matrix (`matrix.mtx`), features/genes (`features.tsv`/`genes.tsv`), and barcodes (`barcodes.tsv`).
- **HDF5-based (H5AD):** AnnData format for Python/Scanpy (`.h5ad`).
- **Loom:** HDF5 format for sparse single-cell data, used by loompy.
- **10x Genomics HDF5:** Filtered feature-barcode matrix (`filtered_feature_bc_matrix.h5`).

## Other Formats
- **Cell Ranger Outputs:** Directory structure including `raw_feature_bc_matrix/` and `filtered_feature_bc_matrix/` with matrices and analysis files.
- **SingleCellExperiment (R):** R data structure saved as `.rds`, storing counts, metadata, and analysis results.
- **Seurat Objects:** RDS files (`.rds` or `.rda`) containing Seurat objects with raw and processed data.

Use appropriate loading functions to import these formats into your analysis environment:
- Python/Scanpy: `sc.read_10x_h5()`, `sc.read_loom()`, `sc.read_h5ad()`
- R/Seurat: `Read10X()`, `Read10X_h5()`, `readRDS()`
