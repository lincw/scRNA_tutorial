# Tips & Pitfalls

This section highlights common challenges and best practices to avoid misinterpretation or technical issues in scRNA-seq analysis.

## Quality Control Pitfalls
- **Over-filtering:** Removing cells with true biological variation when applying strict thresholds (e.g., low gene count, high mitochondrial content).
- **Under-filtering:** Retaining low-quality cells or doublets that introduce noise.

## Batch Effects
- **Ignoring batch effects:** Can lead to clusters driven by technical rather than biological differences.
- **Over-correction:** Using aggressive integration may remove genuine biological signals.

## Clustering
- **Over-clustering vs. Under-clustering:**
  - High resolution may split similar cells into artificial clusters.
  - Low resolution may merge distinct cell types.
- **Parameter tuning:** Test multiple resolutions and neighbor parameters; validate clusters with marker genes.

## Normalization & Scaling
- **Method choice:** Different methods (log-normalization, SCTransform) affect downstream results.
- **Scaling pitfalls:** Incorrect scaling may remove real biological signals or inflate noise.

## Annotation
- **Marker bias:** Relying solely on known markers may miss novel cell types.
- **Validation:** Cross-validate annotations with external references or orthogonal assays (e.g., FACS, immunostaining).

## Pseudotime & Trajectory
- **Root selection:** Pseudotime ordering depends on chosen root cell or cluster.
- **Relative measure:** Pseudotime reflects cell-state progression, not actual chronological time.

## Interpretation & Reporting
- **Overinterpretation:** Avoid conclusions from small or poorly-supported clusters.
- **Plot aesthetics:** Color scales, axis limits, and point size influence perception; use consistent themes.

## Best Practices
- **Document parameters:** Record all thresholds, normalization, clustering, and integration settings.
- **Version control:** Track code and analysis steps with Git for reproducibility.
- **Share data:** Provide raw and processed data alongside code to ensure transparency.

## References
- Luecken MD & Theis FJ (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol, 15(6):e8746.
- Andrews TS & Hemberg M (2019). Identifying cell populations with scRNA-seq data. Nat Biotechnol, 37, 1171–1177.
