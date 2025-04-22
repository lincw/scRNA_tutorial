# Analysis Examples

This directory contains example analysis scripts for the scRNA-seq tutorial. These examples demonstrate complete analysis workflows for the PBMC 3K dataset using both R (Seurat) and Python (Scanpy).

## Files

- **pbmc_analysis.ipynb**: Jupyter notebook with Python/Scanpy analysis
- **pbmc_analysis.Rmd**: R Markdown document with R/Seurat analysis

## Data Requirements

These examples require the PBMC 3K dataset from 10x Genomics. See the `../data/README.md` file for instructions on downloading the dataset.

## Running the Examples

### Python Notebook

To run the Jupyter notebook, use:

```bash
jupyter notebook pbmc_analysis.ipynb
```

Requirements:
- Python 3.6+
- scanpy
- numpy
- pandas
- matplotlib
- seaborn

### R Markdown

To run the R Markdown document, open it in RStudio or use:

```bash
Rscript -e "rmarkdown::render('pbmc_analysis.Rmd')"
```

Requirements:
- R 4.0+
- Seurat (v4+)
- dplyr
- ggplot2
- knitr
- rmarkdown

## Output

Both examples produce visualizations and processed data files that will be saved to the `../data` directory.