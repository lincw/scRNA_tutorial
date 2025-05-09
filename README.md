# Beginner's Guide to Single-Cell Transcriptomics Data Analysis

Welcome to the world of single-cell transcriptomics! This guide is designed for beginners who want to understand and perform single-cell RNA sequencing (scRNA-seq) data analysis.

---

> **Note:** This tutorial was prepared in cooperation with GPT-4.1, o4-mini-high, and Claude 3.7 Sonnet. Last updated: 2025-04-22.

---

## Table of Contents
1. [Getting Started](./01_getting_started.md)
2. [Introduction](./02_introduction.md)
3. [Key Concepts](./03_key_concepts.md)
4. [Data Formats](./04_data_formats.md)
5. [scRNA-seq Library Preparation](./05_library_preparation.md)
6. [Typical Workflow & Analysis Pipeline](./06_typical_workflow.md)
7. [Essential Tools & Software](./07_essential_tools.md)
8. [Step-by-Step Analysis Pipeline](./08_analysis_pipeline.md)
9. [Trajectory & Pseudotime Analysis](./09_trajectory_pseudotime.md)
10. [Visualization & Reporting](./10_visualization_and_reporting.md)
11. [Cell-Cell Communication Analysis](./11_cell_communication.md)
12. [Tips & Pitfalls](./12_tips_pitfalls.md)
13. [Hands-On Examples](./13_hands_on_examples.md)
14. [Advanced Topics](./14_advanced_topics.md)
15. [Glossary & FAQ](./15_glossary_and_faq.md)
16. [Further Reading & References](./16_further_resources.md)

---

Each topic is in its own file for easy editing and expansion. Click on any section above to get started!

## Key Features and Best Practices

This tutorial emphasizes rigorous analytical approaches and best practices for scRNA-seq data analysis:

- **Comprehensive Quality Control**: Including doublet detection with tools like DoubletFinder and Scrublet
- **Modern Normalization Methods**: SCTransform and scran pooling-based normalization alongside traditional approaches
- **Batch Effect Correction**: Clear guidance on correct workflow ordering for batch integration
- **Multi-resolution Clustering**: Testing different resolutions to find optimal cell type granularity
- **Robust Cell Type Annotation**: Using multiple markers and orthogonal validation methods
- **Trajectory Analysis**: Detailed prerequisites and proper workflows for RNA velocity analysis
- **Protein Interaction Network Integration**: Advanced methods for incorporating PPI data into cell communication analysis

These advanced approaches will help ensure more reliable and reproducible results from your scRNA-seq analysis.
