# Example Data

This directory contains example data files used in the tutorial. The main example dataset is the PBMC 3K dataset from 10x Genomics.

## Downloading the PBMC 3K Dataset

You can download the dataset using the following commands:

### Using wget
```bash
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

### Using curl
```bash
curl -O https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

### Using R
```r
download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", 
              destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
```

### Using Python
```python
import requests
import tarfile
import os

url = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
response = requests.get(url, stream=True)
with open("pbmc3k_filtered_gene_bc_matrices.tar.gz", "wb") as file:
    for chunk in response.iter_content(chunk_size=1024):
        if chunk:
            file.write(chunk)

with tarfile.open("pbmc3k_filtered_gene_bc_matrices.tar.gz") as tar:
    tar.extractall()
```

## Expected Directory Structure After Download

After downloading and extracting, you should have the following directory structure:

```
data/
├── filtered_gene_bc_matrices/
│   └── hg19/
│       ├── barcodes.tsv
│       ├── genes.tsv
│       └── matrix.mtx
└── pbmc3k_filtered_gene_bc_matrices.tar.gz
```

## Alternative Formats

10x Genomics also provides this dataset in other formats, such as HDF5. If you prefer working with HDF5, you can download it using:

```bash
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_feature_bc_matrix.h5
```