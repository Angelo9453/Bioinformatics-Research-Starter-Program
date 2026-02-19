# Microarray vs RNA-seq

## Introduction

Microarray and RNA-seq are both methods to measure gene expression.  

**Microarray** uses a chip containing probes for known genes. RNA from a sample is converted to cDNA, which hybridizes to these probes. The intensity of binding quantifies gene expression.  

**RNA-seq** sequences cDNA fragments and counts how often each gene appears. It is more comprehensive but requires mapping sequences to the genome.  

For GEO2R analysis, **Microarray datasets** are required because GEO2R works with probe-based arrays, not raw RNA-seq reads.  

---

## Mapping

A sequence is *mapped* when it is aligned to a specific location in the genome. This allows assignment of reads to their corresponding genes. Accurate mapping is crucial for interpreting which genes are up- or downregulated.  

---

## Dataset: GSE19080

**Title**: Gene expression profiling in patients infected with HTLV-1: Identification of ATL and HAM/TSP-specific genetic profiles  

**Summary**: CD4+ T-cells were isolated from 7 ATL patients, 12 HAM/TSP patients, and 11 asymptomatic carriers. Microarray analysis identified approximately 1039 immune-related genes with differential expression. A refined 122-gene signature distinguishes ATL-specific, HAM/TSP-specific, and common infection-related genes.  

**Design**: T lymphocytes from infected and asymptomatic individuals were analyzed using the human ImmuneArray cDNA array. Stratagene Universal Control (Cy5-labelled) was used in competitive hybridization.  

**Contributors**: Hernandez E, Oliere S  
**Submission**: Nov 18, 2009  
**Last Update**: Jan 17, 2013   

---

## Analysis Steps

Access [GEO2R](https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE19080). Define groups as Healthy and Infected. Click “Analyze” to generate results. Download tables for fold change, p-value, and log₂ expression values. These are used to generate the plots.  

---

## Results and Figures

### Figure 1: Volcano Plot

The Volcano plot displays the statistical significance versus magnitude of expression change. X-axis shows **log₂ fold change**, Y-axis shows **-log10(p-value)**.  

<img width="957" height="737" alt="Screenshot 2026-02-11 070250" src="https://github.com/user-attachments/assets/b50e349a-d0ca-47ca-8003-cc6341f7449c" />

Interpretation:

Genes far right and high are strongly upregulated in infected CD4+ T-cells. These include immune genes, cytokines, and stress response genes.

Genes far left and high are strongly downregulated, representing suppressed metabolism, housekeeping, and differentiation pathways.

The global spread of points shows HTLV-1 infection alters transcription widely, not locally.

Figure 2: MA Plot
The MA plot compares average gene expression (A) to log₂ fold change (M).

<img width="872" height="721" alt="Screenshot 2026-02-11 071805" src="https://github.com/user-attachments/assets/6a0b7334-9284-4745-82dd-f7669c1d5ff6" />

Interpretation:

Red points above zero are upregulated genes. These correspond to immune activation, inflammation, and viral response pathways.

Blue points below zero are downregulated genes, representing suppressed normal cellular functions.

The MA plot confirms system-wide transcriptome reprogramming by HTLV-1, including both upregulation of defense pathways and downregulation of housekeeping functions.

Figure 3: UMAP Plot
UMAP reduces high-dimensional gene expression data into 2D clusters. Healthy samples are green, infected are purple.

<img width="541" height="446" alt="Screenshot 2026-02-11 072710" src="https://github.com/user-attachments/assets/e3e08047-fcef-4042-b29a-e9eb192c79d6" />

Interpretation:

Healthy CD4+ T-cells form a distinct cluster separate from infected cells.

Infected cells split into two subclusters, likely corresponding to ATL and HAM/TSP subtypes.

This confirms that infection produces disease-specific transcriptional profiles, supporting the 122-gene signature classification.

Figure 4: Limma Differential Expression
Limma identifies differentially expressed genes (DEGs) across comparisons.

<img width="425" height="380" alt="Screenshot 2026-02-11 073211" src="https://github.com/user-attachments/assets/9a44cf7b-2c71-465d-85bf-4415cd9d09c2" />

Interpretation:

773 genes are shared across comparisons, representing common infection responses.

2229 genes are condition-specific, distinguishing ATL or HAM/TSP.

These data explain the finer structure in UMAP clusters and support functional interpretation of disease-specific gene expression changes.

Figure 5: Boxplot of GSE19080
This boxplot shows normalized expression for each sample. X-axis lists sample IDs, Y-axis shows standardized expression values.

<img width="808" height="483" alt="Screenshot 2026-02-11 073755" src="https://github.com/user-attachments/assets/8ae60adb-3c14-4033-9808-695e8951edaa" />

Interpretation:

Healthy samples (green) have similar median expression and low variance.

Infected samples (purple) show greater spread, indicating heterogeneity in gene expression.

This variation corresponds to ATL- and HAM/TSP-specific regulation and supports the 122-gene signature.

Conclusion
HTLV-1 infection reprograms CD4+ T-cells by activating immune, inflammatory, and viral response genes while suppressing normal cellular pathways. Volcano and MA plots demonstrate the magnitude and direction of these changes, UMAP confirms the transcriptional clustering of healthy versus infected cells, Limma identifies shared and disease-specific DEGs, and boxplots show sample-level variability. These analyses provide a molecular explanation for HTLV-1-driven pathology in ATL and HAM/TSP.
