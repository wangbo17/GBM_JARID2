### Overview
This R-based workflow, developed following the methodology outlined in Tanner et al., 2024, provides a reproducible framework for analysing transcriptional changes between primary and recurrent IDHwt glioblastomas, enabling the stratification of patients into distinct response groups based on their transcriptional profiles.

<img src="https://github.com/user-attachments/assets/81fc50ac-8c25-42cc-8bdf-356fe1e1e488" alt="NES_vs_PC1_plot" width="90%"/>

---

### Analysis Workflow
#### **Step 0: Data Preparation**
- **Input**:  
  Gene annotation file (e.g., `.gtf`), metadata file, and raw gene expression data.  
- **Process**:  
  Gene lengths are calculated using genome annotation. Metadata is filtered to include relevant paired samples, and raw expression data is aligned with metadata and annotated genes.  
- **Output**:  
  Filtered metadata, processed expression data, and gene length information.  

#### **Step 1: Normalization**
- **Input**:  
  Processed expression data and gene length information.  
- **Process**:  
  Raw counts are converted to FPKM values using the edgeR package, accounting for sequencing depth through CPM normalization and adjusting for gene lengths to ensure comparability across genes.
- **Output**:  
  Normalized FPKM matrix.  

#### **Step 2: Log2 Fold Change Calculation**
- **Input**:  
  Normalized FPKM matrix and filtered metadata for paired samples.  
- **Process**:  
  Log2 fold changes (log2FC) are computed by comparing gene expression between paired conditions (e.g., primary vs. recurrent or untreated vs. treated). A small constant (+0.01) is added to avoid division by zero.  
- **Output**:  
  Log2FC values for each gene and patient.  

#### **Step 3: Gene Set Enrichment Analysis**
- **Input**:  
  Log2FC values and a gene set file (e.g., `.gmt`).  
- **Process**:  
  GSEA is performed for the JARID2-related gene set. Genes are ranked by log2FC values, and normalized enrichment scores (NES) are calculated using the `fgsea` package.  
- **Output**:  
  NES results for each sample pair.  

#### **Step 4: Principal Component Analysis**
- **Input**:  
  Log2FC values and precomputed PCA rotation data.  
- **Process**:  
  Genes overlapping between PCA rotation data and log2FC results are used to compute PC1 scores. Log2FC values are  mapped onto an existing PCA space, using the PCA rotation data to calculate PC1 scores. 
- **Output**:  
  PC1 scores for each sample pair.  

#### **Step 5: Visualization**
- **Input**:  
  NES results and PC1 scores.  
- **Process**:  
  NES and PC1 scores are combined to stratify samples into response groups. Patients are categorized as “Up” and “Down” responders based on NES values. A scatter plot is generated to illustrate the relationship between NES and PC1 scores.  
- **Output**:  
  High-resolution scatter plot (`NES_vs_PC1_plot.png`).  

---

### References
Tanner, G., Barrow, R., Ajaib, S., Al-Jabri, M., Ahmed, N., Pollock, S., Finetti, M., Rippaus, N., Bruns, A. F., Syed, K., Poulter, J. A., Matthews, L., Hughes, T., Wilson, E., Johnson, C., Varn, F. S., Brüning-Richardson, A., Hogg, C., Droop, A., Gusnanto, A., Care, M. A., Cutillo, L., Westhead, D. R., Short, S. C., Jenkinson, M. D., Brodbelt, A., Chakrabarty, A., Ismail, A., Verhaak, R. G. W. & Stead, L. F. (2024) 'IDHwt glioblastomas can be stratified by their transcriptional response to standard treatment, with implications for targeted therapy', *Genome Biology*, vol. 25, p. 45.

