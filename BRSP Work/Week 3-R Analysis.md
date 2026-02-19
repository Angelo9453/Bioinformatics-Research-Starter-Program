# Lung Cancer Gene Expression Analysis

Dataset: GSE10072
Disease: Lung Adenocarcinoma vs Normal Lung Tissue
Platform: Affymetrix Human Genome U133A (GPL96)
Method: Differential Expression Analysis using limma

## 1. Project Overview

This project performs gene expression analysis using microarray data from the public database Gene Expression Omnibus (GEO).

We analyze dataset GSE10072, which contains:

-Lung adenocarcinoma tumor samples

-Normal lung tissue samples

---

## Objectives

Identify Differentially Expressed Genes (DEGs)

Visualize gene expression differences

Perform GO enrichment analysis

Perform KEGG pathway enrichment analysis

Interpret the biological meaning of results

---

## 2. What is Gene Expression Analysis?

Gene expression analysis compares how strongly genes are expressed between two biological conditions.

In this study:

Condition 1: Lung adenocarcinoma (cancer)

Condition 2: Normal lung tissue

We use the limma method (Linear Models for Microarray Data), which is a widely used statistical framework for microarray analysis.

---

## 3. Complete R Code

The full reproducible R workflow is provided below.

## 3.1 Install and Load Packages
### Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

### Install Bioconductor packages
BiocManager::install(c("GEOquery", "limma", "hgu133a.db"), 
                     ask = FALSE, update = FALSE)

### Install CRAN packages
install.packages(c("pheatmap", "ggplot2", "dplyr"))

if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

### Load libraries
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

## 3.2 Download Data from GEO
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

ex <- exprs(gset)

## 3.3 Log2 Transformation
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

## 3.4 Define Sample Groups
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)

gset$group <- factor(groups)
group_levels <- levels(gset$group)

design <- model.matrix(~0 + gset$group)
colnames(design) <- group_levels

contrast_formula <- paste(group_levels[1], "-", group_levels[2])
print(group_levels)


<img width="658" height="35" alt="Screenshot 2026-02-19 091220" src="https://github.com/user-attachments/assets/daa642cf-4799-440e-b8cb-f4e3d9ea4538" />

**Figure 1**
This shows both Samples have identfied groups

## 3.5 Differential Expression Analysis (limma)
fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, 
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)
head(topTableResults)

<img width="929" height="223" alt="Screenshot 2026-02-19 091435" src="https://github.com/user-attachments/assets/026210a7-280d-4dd0-ba44-463400f61f77" />

**Figure 2** 
This is the output table from limma differential expression analysis

## 3.6 Gene Annotation
probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

<img width="861" height="210" alt="Screenshot 2026-02-19 091738" src="https://github.com/user-attachments/assets/cf854925-7ca2-40d7-99cf-a6164eba3897" />

**Figure 3** This shows the PROBEID, SYMBOL, and GENENAME

## 3.7 BoxPlot
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Expression Value Distribution per Sample",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

<img width="1919" height="957" alt="Screenshot 2026-02-19 111842" src="https://github.com/user-attachments/assets/e1bba13b-43d0-4b8b-b46b-4595d33bad6f" />

**Figure 4** 
This figure shows the distribution of gene expression values for each sample in the dataset.
The title, “Expression Value Distribution per Sample”, means that each boxplot represents how gene expression values are spread within one sample. The y-axis shows expression values on a log2 scale, which is a common transformation used in gene expression studies to reduce large differences and make the data easier to compare.
Each vertical box represents one sample (labeled along the x-axis, e.g., GSM IDs). Black boxes represent adenocarcinoma of the lung (cancer) samples. Pink boxes represent normal lung tissue samples.

For each boxplot:

-The middle horizontal line inside the box is the median (the middle expression value).

-The box shows the interquartile range (IQR), meaning the middle 50% of gene expression values.

-The whiskers extend to the lower and upper values, showing the overall spread of the data.

From the plot, we can see that:

-The medians are very similar across all samples, both cancer and normal. Most medians are around 7–7.5 (log2 scale).

-The spread of values (box size and whisker length) is also similar between cancer and normal samples.

-There are no obvious outlier samples with very different distributions.

This suggests that:

-The data appear to be well normalized, since the overall expression distributions are consistent across samples.

-There are no major technical biases or batch effects visible at the global level.

## 3.8 Density Plot
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Gene Expression Distribution",
    x = "Expression Value (log2)",
    y = "Density"
  )

<img width="746" height="656" alt="Screenshot 2026-02-19 112351" src="https://github.com/user-attachments/assets/0e84c630-8e7b-4b8e-bff0-0affdbb7663c" />

**Figure 5** This figure shows the overall distribution of gene expression values for two groups: lung adenocarcinoma samples and normal lung tissue samples.
The x-axis represents gene expression values on a log2 scale. This transformation is commonly used in gene expression analysis because it reduces large differences between values and makes the data easier to compare. The y-axis represents density, which shows how frequently certain expression values occur. Higher peaks indicate that many genes have expression levels around that value.

There are two curves in the plot:

-The red curve represents adenocarcinoma of the lung samples.

-The blue curve represents normal lung tissue samples.

Both curves have a very similar shape and almost completely overlap. This means that, at a global level, the overall distribution of gene expression values is very similar between cancer and normal samples.The peak of both curves is around 7–8 (log2 scale), indicating that most genes have expression values in this range. The distribution is slightly right-skewed, meaning there are fewer genes with very high expression values (above 10–12), but some genes do show higher expression levels, forming a long tail to the right.

The strong overlap between the two curves suggests that:

-The data are well normalized.

-There is no major global shift in overall gene expression between cancer and normal tissues.

## 3.9 Uniform Manifold Approximation and Projection (UMAP)
umap_input <- t(ex)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot of Samples Based on Gene Expression",
    x = "UMAP 1",
    y = "UMAP 2"
  )

<img width="749" height="653" alt="Screenshot 2026-02-19 112702" src="https://github.com/user-attachments/assets/1d781e23-e44a-4d6a-8050-31aa39555556" />

**Figure 6** This figure shows a UMAP plot of all samples based on their gene expression profiles. UMAP (Uniform Manifold Approximation and Projection) is a dimensionality reduction method that summarizes thousands of gene expression values into two components (UMAP 1 and UMAP 2) so that similarities between samples can be visualized in a two-dimensional plot. Each point represents one sample. The position of each point is determined by its overall gene expression pattern. Samples that are close to each other have similar gene expression profiles, while samples that are far apart are more different.

In the plot, two clear clusters can be seen:

-The red points represent lung adenocarcinoma samples.

-The blue points represent normal lung tissue samples.

Most cancer samples group together on the right side of the plot, while most normal samples cluster together on the left side. This clear separation indicates that gene expression patterns are substantially different between lung adenocarcinoma and normal lung tissue. This supports the presence of strong biological differences between the two groups.

However, there appears to be one adenocarcinoma sample located within or very close to the normal lung tissue cluster. This can happen for several possible reasons:

-Biological heterogeneity: Tumors are not all identical. Some cancer samples may have gene expression patterns that are more similar to normal tissue, especially if the tumor is early-stage, less aggressive, or contains a high proportion of normal cells.

-Tumor purity: The sample may contain a large amount of surrounding normal tissue mixed with tumor cells. This lowers the overall cancer-specific signal and makes the expression profile resemble normal tissue.

-Technical variation: Small technical differences during sample processing, RNA extraction, or sequencing can slightly affect gene expression patterns.

## 3.9 Volcano Plot
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Lung Cancer DEG Volcano Plot")

<img width="780" height="659" alt="Screenshot 2026-02-19 112934" src="https://github.com/user-attachments/assets/45cb9433-86c0-4090-9e97-e6817182b25c" />

**Figure 7** This figure shows a volcano plot of differentially expressed genes (DEGs) between lung adenocarcinoma and normal lung tissue samples.
Each point represents one gene. 

The position of each gene is determined by two values:

-The x-axis shows the log fold change (logFC). This indicates how much a gene’s expression differs between cancer and normal tissue.

-Genes on the right side (positive logFC) are upregulated in lung cancer, meaning they are expressed at higher levels in cancer compared to normal tissue.

-Genes on the left side (negative logFC) are downregulated in lung cancer, meaning they are expressed at lower levels in cancer.


The y-axis shows −log10(adjusted p-value). This represents the statistical significance of the difference in expression.
Higher values on the y-axis indicate more statistically significant results. Genes near the bottom have less significant differences. The dashed vertical lines indicate the fold change cutoff (for example, logFC = ±1). Genes beyond these lines show a biologically meaningful change in expression. The dashed horizontal line indicates the adjusted p-value cutoff (for example, 0.05). Genes above this line are statistically significant after correcting for multiple testing. 

The colors represent gene status:

-Red points (UP) are significantly upregulated genes in lung cancer.

-Blue points (DOWN) are significantly downregulated genes in lung cancer.

-Grey points (NO) are genes that are not significantly different.


From the plot, we can see that:

-There are many significantly differentially expressed genes between cancer and normal tissue.

-Both upregulated and downregulated genes are present.

-Some genes show very strong statistical significance (very high on the y-axis) and large fold changes, indicating strong differences between the two groups.

## 4.0 Heatmap of Top 50 Genes
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]
gene_label <- ifelse(is.na(top50$SYMBOL) | top50$SYMBOL == "", top50$PROBEID, top50$SYMBOL)
rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

<img width="884" height="794" alt="Screenshot 2026-02-19 113745" src="https://github.com/user-attachments/assets/aac7a0db-c2cf-42a0-b43c-9a537ba1df0f" />
**Figure 8**  This figure shows a heatmap of the top 50 differentially expressed genes between lung adenocarcinoma and normal lung tissue samples.
Each row represents one gene, and each column represents one sample. The genes shown are the 50 most statistically significant genes identified from the differential expression analysis. The expression values have been scaled by row, meaning each gene’s expression is standardized across samples. This allows easier comparison of relative expression differences between cancer and normal tissues.

The color scale represents relative expression levels:

-Red indicates higher expression.

-Blue indicates lower expression.


Yellow represents intermediate expression.


At the top of the heatmap, there is a sample annotation bar labeled “Group.” The group colors were manually assigned as:
 
 "Adenocarcinoma.of.the.Lung" = "pink"
 
 "Normal.Lung.Tissue" = "cyan"
 
This means that pink columns represent lung adenocarcinoma samples, and cyan columns represent normal lung tissue samples.
The heatmap shows a clear separation between most cancer and normal samples. The majority of adenocarcinoma samples cluster together and show similar expression patterns, while most normal samples also cluster together in a separate group. This indicates that the top 50 genes strongly distinguish between cancer and normal tissues.
In general, many genes appear highly expressed (red) in cancer samples but low in normal samples (blue), while others show the opposite pattern. This confirms strong biological differences between the two groups.

However, one normal lung tissue sample appears within the adenocarcinoma cluster. This can happen for several possible reasons:

-Biological variability: Not all normal tissues are identical. Some normal samples may naturally have gene expression patterns that partially resemble tumor tissue.

-Sample heterogeneity: The normal sample may contain a small proportion of abnormal or pre-cancerous cells, which can influence its expression profile.

## 4.1 Gene Ontology Analysis (KO) AND Kyoto Encyclopedia of Genes and Genomes (KEGG)

## 4.1.1 Differential Gene Filtering (Input for Enrichment)

**This section prepares the significant genes that will be used for GO and KEGG**
deg_sig <- topTableResults %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

gene_symbols <- deg_sig$SYMBOL
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

gene_entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

## 4.1.2 GO Analysis
go_results <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",   # BP = Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

dotplot(go_results, showCategory = 10) +
  ggtitle("GO Biological Process Enrichment")

<img width="612" height="460" alt="Screenshot 2026-02-19 114937" src="https://github.com/user-attachments/assets/7519261e-1a2c-49f1-bb80-d55914c3d64e" />
**Figure 9** To investigate the biological mechanisms underlying the differences between lung adenocarcinoma tumor samples and normal lung tissue, Gene Ontology (GO) Biological Process enrichment analysis was performed using significantly differentially expressed genes (FDR < 0.05 and |logFC| > 1).

The analysis identified several highly significant biological processes that are strongly enriched in tumor samples. The most enriched categories include extracellular matrix organization, extracellular structure organization, external encapsulating structure organization, wound healing, regulation of angiogenesis, regulation of vasculature development, and cell–substrate adhesion. These processes demonstrated extremely strong statistical significance (adjusted p-values < 10⁻¹⁵) and involved a large number of genes (approximately 44–52 genes per term), indicating robust and coordinated biological alterations in tumor tissue.

1.Extracellular Matrix Remodeling

The strongest enrichment was observed in extracellular matrix (ECM)-related processes. The extracellular matrix is the structural scaffold that surrounds cells and provides mechanical and biochemical support. In cancer, remodeling of the ECM is a central feature of tumor progression. The enrichment of ECM organization suggests that tumor tissue undergoes increased collagen production, structural reorganization of surrounding tissue, degradation of basement membrane barriers, and increased tissue stiffness. These structural alterations enable tumor cells to invade adjacent tissues, infiltrate surrounding structures, and ultimately metastasize. ECM remodeling is widely recognized as a hallmark of aggressive adenocarcinoma and reflects active modification of the tumor microenvironment.

2.Wound Healing–Like Processes

Tumors are often described as “wounds that do not heal.” The enrichment of wound healing pathways indicates activation of fibroblasts, growth factor signaling, inflammatory mediators, and tissue repair mechanisms. Under normal physiological conditions, these processes promote recovery after injury. However, in cancer, they are hijacked to create a microenvironment that supports tumor survival and expansion. The activation of wound healing programs suggests enhanced stromal activation, modulation of immune responses, and increased production of growth-promoting signals. This supports the concept that tumor progression is not driven solely by malignant cells, but also by dynamic interactions between cancer cells and their surrounding tissue environment.

3.Angiogenesis and Vasculature Development

Significant enrichment was also observed in biological processes related to angiogenesis and regulation of vasculature development. Angiogenesis refers to the formation of new blood vessels, a process that tumors rely on for sustained growth. As tumors expand, they require increased oxygen and nutrient supply while simultaneously removing metabolic waste. Activation of angiogenic pathways enables the development of new vascular networks that support rapid tumor proliferation. Additionally, enhanced vascularization increases the potential for tumor cells to enter circulation, thereby facilitating metastasis. The strong enrichment of angiogenesis-related processes indicates that the adenocarcinoma samples exhibit active vascular remodeling consistent with aggressive tumor behavior.

4.Cell–Substrate Adhesion

Cell–substrate adhesion processes were also significantly enriched. These processes regulate how cells attach to the extracellular matrix through molecules such as integrins and focal adhesion proteins. Alterations in cell adhesion are critical for cancer progression because they influence cell migration, invasion, and structural organization. Dysregulation of adhesion pathways allows tumor cells to detach from their original tissue architecture, gain motility, and invade surrounding tissues. This contributes directly to metastatic potential and loss of normal tissue integrity.

Interpretation of the GO Enrichment Plot

The GeneRatio values (approximately 0.07–0.08) indicate that roughly 7–8% of the significant differentially expressed genes are involved in each enriched biological process. The relatively large dot sizes observed in the enrichment plot reflect the high number of contributing genes (44–52 genes per term), further supporting the robustness of these findings. The extremely low adjusted p-values confirm that these enrichments are statistically strong and highly unlikely to have occurred by chance.

## 4.1.3 KEGG Analysis
kegg_results <- enrichKEGG(
  gene         = gene_entrez$ENTREZID,
  organism     = "hsa",   # human
  pvalueCutoff = 0.05
)

kegg_results <- setReadable(
  kegg_results,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

dotplot(kegg_results, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment")
  
<img width="535" height="456" alt="Screenshot 2026-02-19 115258" src="https://github.com/user-attachments/assets/c3a599f1-7ca9-4746-815b-5a622df39e04" />
**Figure 10** To further characterize functional pathway alterations between lung adenocarcinoma tumor samples and normal lung tissue, KEGG pathway enrichment analysis was performed using significantly differentially expressed genes (FDR < 0.05 and |logFC| > 1). The analysis identified several significantly enriched pathways, including cytoskeleton organization–related pathways, focal adhesion, integrin signaling, complement and coagulation cascades, relaxin signaling, AGE–RAGE signaling pathway in diabetic complications, protein digestion and absorption, ECM–receptor interaction, and the p53 signaling pathway. These pathways demonstrated strong statistical significance (adjusted p-values ranging from 10⁻³ to 10⁻⁶) and involved approximately 10–30 genes per pathway, indicating coordinated molecular alterations in tumor samples.

1.Focal Adhesion and Integrin Signaling

Focal adhesion and integrin signaling pathways regulate how cells attach to the extracellular matrix and transmit signals between the extracellular environment and intracellular signaling networks. These pathways control cell movement, survival, and structural organization. Their enrichment suggests increased adhesion dynamics and enhanced migration capacity in tumor cells. In adenocarcinoma, integrins and focal adhesion proteins facilitate detachment from normal tissue architecture, migration into surrounding tissue, and initiation of metastatic processes. Activation of these pathways is strongly associated with invasive tumor behavior.

2.ECM–Receptor Interaction

The ECM–receptor interaction pathway directly reflects communication between tumor cells and the extracellular matrix. Enrichment of this pathway indicates structural remodeling of tissue and altered cell–matrix signaling.This suggests that tumor cells are actively modifying their surrounding microenvironment to promote invasion and survival. Altered ECM signaling enhances motility, facilitates tissue infiltration, and supports metastatic potential. This finding reinforces the GO results highlighting extracellular matrix remodeling as a central biological theme.

3.Complement and Coagulation Cascades

The complement and coagulation cascades pathway is associated with immune response and blood clotting mechanisms. Its enrichment indicates activation of tumor-associated inflammatory processes and potential modulation of the immune microenvironment. In cancer, activation of coagulation pathways can create a pro-thrombotic state that supports tumor progression, enhances metastatic spread, and protects circulating tumor cells from immune detection. This suggests that adenocarcinoma may interact extensively with immune and vascular systems.

3.p53 Signaling Pathway

The p53 signaling pathway is a major tumor suppressor pathway that regulates cell cycle arrest, DNA repair, and apoptosis. Enrichment of this pathway suggests altered regulation of cellular stress responses and genomic stability mechanisms.In lung adenocarcinoma, disruption of p53 signaling is common and contributes to uncontrolled proliferation and impaired apoptosis. The identification of this pathway supports the presence of tumor suppressor dysregulation and DNA damage–associated signaling in tumor samples.

4.Cytoskeleton Organization–Related Pathways

Although labeled as “cytoskeleton in muscle cells,” this pathway reflects broader cytoskeletal organization and structural protein regulation. Cytoskeletal remodeling plays a critical role in changes in cell shape, polarity, and motility. Alterations in cytoskeletal dynamics are essential for tumor invasion and metastatic dissemination. Enrichment of cytoskeleton-related pathways suggests enhanced structural reorganization that enables tumor cells to migrate and invade surrounding tissue.

5.AGE–RAGE and Relaxin Signaling Pathways

The AGE–RAGE and relaxin signaling pathways are associated with inflammatory signaling, oxidative stress, and extracellular matrix remodeling. In cancer, activation of these pathways can promote tumor progression, fibrosis-like tissue changes, angiogenesis, and microenvironment activation. Their enrichment indicates that tumor samples may experience heightened oxidative stress and inflammatory signaling, further supporting aggressive tumor behavior.

Interpretation of the KEGG Enrichment Plot

The GeneRatio values (approximately 0.03–0.08) indicate that 3–8% of significant differentially expressed genes are involved in each enriched pathway. Larger dots in the enrichment plot correspond to pathways containing a greater number of contributing genes (up to ~30 genes). Red-colored dots represent stronger statistical significance, reflecting lower adjusted p-values. Cytoskeleton organization–related pathways appear among the most statistically significant, followed closely by focal adhesion and integrin signaling pathways, indicating that structural remodeling and cell–matrix communication are dominant molecular features of the tumor samples.

## 4.2 Save Results
write.csv(topTableResults, "GSE10072_DEG_Results.csv")

