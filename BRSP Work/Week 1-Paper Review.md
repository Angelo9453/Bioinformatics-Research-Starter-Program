# Bioinformatics-Research-Starter-Program-OmicsLite
Program design to work on projects, ameliorate research critical thinking

# Paper could be access here: https://atm.amegroups.org/article/view/95790/html 

# Transcriptomic Analysis of Sepsis-Associated Acute Kidney Injury (SA-AKI)

## Overview
Sepsis-associated acute kidney injury (SA-AKI) is a common and life-threatening complication of sepsis with limited targeted therapies. The molecular mechanisms driving renal injury remain incompletely understood.  

This project analyzes kidney transcriptomic changes in a murine sepsis model and integrates bioinformatics with laboratory validation to identify critical genes, pathways, and RNA regulatory mechanisms involved in SA-AKI.

---

## Objectives
- Identify differentially expressed genes (DEGs) in septic kidney tissue  
- Determine enriched biological pathways and molecular functions  
- Detect hub genes through protein–protein interaction (PPI) networks  
- Construct miRNA–mRNA–lncRNA ceRNA regulatory networks  
- Validate key targets experimentally  
- Investigate the role of m6A RNA methylation regulators  

---

## Experimental Design

### Animal Model
- Species: C57BL/6 mice  
- Sepsis induction: Cecal ligation and puncture (CLP)  
- Groups: CLP vs Control  
- Tissue: Kidney  

### Transcriptomic Analysis
- RNA sequencing (RNA-seq)
- Differential expression analysis
- Functional enrichment:
  - Gene Ontology (GO)
  - KEGG pathways
  - Gene Set Enrichment Analysis (GSEA)
- Protein–protein interaction (PPI) networks
- Hub gene identification

### Regulatory Network Analysis
- miRNA–mRNA–lncRNA ceRNA network construction
- m6A RNA methylation regulator profiling

### Experimental Validation
- qPCR
- Western blot
- ROS detection
- Malondialdehyde (MDA) assay
- Flow cytometry
- Cell model: TCMK-1 renal tubular epithelial cells (LPS stimulation)

---

## Key Results

### Differential Expression
- **4,754 DEGs** identified between CLP and control groups

### Enriched Pathways
Top pathways:
- Ferroptosis (most significant)
- Apoptosis
- PI3K–Akt signaling
- NF-κB signaling
- IL-17 signaling

### Hub Genes
Fifteen hub genes identified, including:
- Hmox1
- Lcn2
- Cxcl1
- Cxcl12
- Mapk14
- Socs3
- Spp1

### Ferroptosis Mechanism
- **Hmox1** identified as a central ferroptosis gene
- Regulated by **mmu-miR-7212-5p**

### Functional Validation
miR-7212-5p inhibition:
- ↑ Hmox1 expression
- ↓ oxidative stress
- ↓ ferroptosis
- ↓ inflammatory mediator release

### Epigenetic Regulation
m6A regulators:
- METTL3
- ALKBH5

These showed significant expression changes and correlation with hub genes.

---

## Conclusions
SA-AKI pathogenesis is strongly associated with:
- Immune activation
- Ferroptosis
- RNA regulatory mechanisms

Key regulatory axis:
miR-7212-5p → Hmox1 → ferroptosis modulation

m6A RNA methylation may further influence disease progression.

These pathways represent potential therapeutic targets.

---

## Methods & Results Diagram
<img width="416" height="592" alt="image" src="https://github.com/user-attachments/assets/f18a9d2e-88ef-4104-b024-2ddb03b6786b" />
