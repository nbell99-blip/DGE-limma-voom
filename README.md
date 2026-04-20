# RNA-seq Differential Gene Expression Analysis
## Night Shift Gene Expression in Healthcare Workers (GSE282051)

---

## Study Overview

This repository contains the code and results for a differential gene expression (DGE) analysis of RNA-seq data from healthcare workers before and after a night shift.

**Publication:** Changes in Gene Expression in Healthcare Workers During Night Shifts: Implications for Immune Response and Health Risks  
**GEO Accession:** [GSE282051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282051)  
**PubMed ID:** 40069845  
**Organism:** Homo sapiens  
**Tissue:** Whole blood  
**Platform:** Illumina NovaSeq 6000  

---

## Experimental Design

- 4 medical doctors working night shifts in an emergency department
- Blood samples collected at 10am **before** and 10am **after** a night shift
- Paired design — each doctor sampled twice (n=4 subjects, 8 samples total)
- RNA sequencing performed on whole blood

| Sample | Condition | Subject |
|--------|-----------|---------|
| before_1 | Before night shift | Doctor 1 |
| after_1 | After night shift | Doctor 1 |
| before_2 | Before night shift | Doctor 2 |
| after_2 | After night shift | Doctor 2 |
| before_3 | Before night shift | Doctor 3 |
| after_3 | After night shift | Doctor 3 |
| before_4 | Before night shift | Doctor 4 |
| after_4 | After night shift | Doctor 4 |

---

## Repository Contents

```
├── GSE282051_DGE_SOP.R          # Full R script SOP
├── GSE282051_DGE_analysis.Rmd   # R Markdown analysis file
├── GSE282051_DGE_results.csv    # Final DGE results table
├── GSE282051_volcano_plot.png   # Volcano plot
├── GSE282051_heatmap.png        # Heatmap of top 50 DE genes
└── README.md                    # This file
```

---

## Requirements

### R Version
R 4.0 or higher recommended

### R Packages

**Bioconductor packages:**
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", 
                       "edgeR", 
                       "BiocParallel", 
                       "AnnotationDbi", 
                       "org.Hs.eg.db"))
```

**CRAN packages:**
```r
install.packages(c("tidyverse", 
                   "RColorBrewer", 
                   "gplots", 
                   "ggplot2"))
```

---

## How to Reproduce This Analysis

### Step 1 — Download raw count data from GEO
Go to [GSE282051](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282051) and download the 8 individual raw count files:
- GSM8634930_before_1_Raw_Count.txt.gz
- GSM8634931_after_1_Raw_Count.txt.gz
- GSM8634932_before_2_Raw_Count.txt.gz
- GSM8634933_after_2_Raw_Count.txt.gz
- GSM8634934_before_3_Raw_Count.txt.gz
- GSM8634935_after_3_Raw_Count.txt.gz
- GSM8634936_before_4_Raw_Count.txt.gz
- GSM8634937_after_4_Raw_Count.txt.gz

### Step 2 — Download HGNC complete set
The analysis will automatically download this file, or you can download it manually from:  
https://www.genenames.org/download/statistics-and-files/

### Step 3 — Run the analysis
Open `GSE282051_DGE_SOP.R` or `GSE282051_DGE_analysis.Rmd` in RStudio and run all sections in order.

---

## Analysis Pipeline

The analysis follows a standard limma-voom workflow:

1. **Data loading** — Load raw count matrix and sample metadata
2. **Quality control** — Library sizes, count distributions, density plots
3. **Gene ID validation** — Validate HGNC symbols against org.Hs.eg.db
4. **DGEList object** — Create edgeR DGEList object
5. **Filtering** — Remove lowly expressed genes using filterByExpr()
6. **TMM normalization** — Correct for compositional bias between samples
7. **Design matrix** — Paired design blocking by subject
8. **Voom transformation** — Convert counts to log2 CPM with precision weights
9. **Sample QC** — MDS plot, hierarchical clustering, correlation heatmap
10. **Linear model** — Fit linear model using lmFit()
11. **Contrasts** — Define after vs before contrast using makeContrasts()
12. **eBayes** — Apply empirical Bayes smoothing
13. **Results** — Extract DE results using topTable()
14. **HGNC annotation** — Join results with HGNC complete set
15. **Visualization** — Volcano plot and heatmap

---

## Results Summary

| Metric | Value |
|--------|-------|
| Total genes tested | 16,338 |
| Genes retained after HGNC join | 14,635 |
| Upregulated after night shift (P < 0.05, \|logFC\| > 0.5) | 381 |
| Downregulated after night shift (P < 0.05, \|logFC\| > 0.5) | 119 |
| Genes passing FDR < 0.05 | 0 |

**Note on statistical significance:** No genes passed FDR-adjusted p-value threshold (adj.P.Val < 0.05) due to the small sample size (n=4). Results are reported using raw p-values with |logFC| > 0.5 threshold, consistent with the original publication. Both raw and adjusted p-values are provided in the results table.

---

## Output Files

### GSE282051_DGE_results.csv
Final results table with 4 columns:

| Column | Description |
|--------|-------------|
| HGNC_Symbol | Official HGNC gene symbol |
| logFC | Log2 fold change (positive = higher after shift) |
| raw_P.value | Raw unadjusted p-value |
| adj_P.value | BH-adjusted p-value (FDR corrected) |

---

## Key QC Notes

- **after_4** was consistently flagged across all QC metrics (higher zero counts, lowest normalization factor, outlier in hierarchical clustering and correlation heatmap) but was retained due to small sample size (n=4)
- Cook's distance analysis was uninformative due to model saturation with n=4 paired samples
- Counts had decimal values due to featureCounts `-M -O --fraction` option and were rounded prior to analysis

---

## Citation

If you use this code, please cite the original study:

> Nukiwa R, Oda S, Matsumoto H, Motooka D. Changes in Gene Expression in Healthcare Workers During Night Shifts: Implications for Immune Response and Health Risks. 2025. PMID: 40069845

---

## Contact

For questions about this analysis, please open an issue in this repository.
