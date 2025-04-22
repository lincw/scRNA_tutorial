# scRNA-seq Library Preparation

Library preparation is a critical step that determines data quality, sensitivity, and cost. Below we compare popular protocols, outline key steps, and highlight QC metrics.

## Protocol Comparison

| Protocol               | Platform      | Throughput      | Sensitivity   | Cost per cell | Coverage                |
|------------------------|---------------|-----------------|---------------|---------------|-------------------------|
| **Smart-seq2**         | Plate-based   | ~384 cells/day  | High          | ~$5           | Full-length transcripts |
| **10x Genomics**       | Droplet-based | ~10K cells/run  | Medium        | ~$1           | 3′ or 5′ end counting   |
| **Drop-seq**           | Droplet-based | ~10K cells/run  | Low–Medium    | Low           | 3′ end                   |

## Key Steps

1. **Cell Isolation**  
   - Plate-based (e.g., FACS into 96/384-well plates) or droplet-based microfluidics (10x, Drop-seq).
2. **Cell Lysis & RNA Capture**  
   - Cells are lysed; mRNA captured via oligo-dT primers attached to beads or wells.
3. **Reverse Transcription & Barcoding**  
   - cDNA synthesis with incorporation of cell barcodes and UMIs.
4. **Amplification**  
   - PCR to amplify cDNA, balancing yield and bias.
5. **Library Construction**  
   - Adapters added, libraries size-selected and quantified.

## QC Metrics at Library Prep

- **mRNA Capture Rate**: Proportion of transcripts captured; measured with ERCC spike-ins or known housekeeping genes.
- **Library Complexity**: Number of unique transcripts (UMIs) per cell.
- **Amplification Bias**: Variation introduced during PCR, assessed by UMI duplication rates.
- **Cell Viability**: Pre-lysis viability (>85% recommended) impacts cDNA yield.

## References

- Picelli S, et al. (2013). Smart-seq2 for sensitive full-length transcriptome profiling in single cells. Nat Methods, 10(11), 1096–1098.
- Macosko EZ, et al. (2015). Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell, 161(5), 1202–1214.
- Zheng GX, et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nat Commun, 8, 14049.
