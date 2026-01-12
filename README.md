# Impact of TMM and Median normalization method on RNA‑seq differential gene expression analysis 

Normalization is a pre-processing step in RNA-seq data analysis that transforms raw counts into values suitable for inter-sample comparisons. Different normalization methods, such as Trimmed Mean of M-values (TMM), Median, and others, have been developed to address biases in RNA-seq data. Despite their importance, the choice of normalization method is often missed, and its impact on differential expression analysis remains not completely investigated. Recent studies, including Tong et al., 2020, have highlighted the influence of normalization choices on computational speed, precision, and accuracy, underscoring the need to evaluate these methods systematically.

This study aims to investigate the impact of normalization methods on RNA-seq differential gene expression analysis. Specifically, it compares two widely used normalization methods “TMM and Median normalization” on their ability to identify differentially expressed (DE) genes in a publicly available RNA-seq dataset. By analyzing differences in DE gene detection, computational performance, and result reproducibility, this study seeks to address the following questions: How do the TMM and Median normalization methods differ in their impact on the identification of DE genes? Are findings from differential expression analyses consistent between the two methods? 

## Description of Data Set 
For this study, I analyzed a publicly available RNA-seq dataset obtained from the Gene Expression Omnibus (GEO) database at NCBI (Accession: GSE132698). The dataset was downloaded on December 12, 2024, directly from the GEO platform using download link from R utilities. The dataset originates from a study by Wang et al., 2019, which investigated the role of ASCL1, a basic helix-loop-helix transcription factor, in the oncogenic pathways of neuroblastoma a type of cancer associated with mutations in regulatory sequences of the LMO1 gene. The data include four raw count files representing two experimental groups: untreated cells as a control (shGFP) and knockdown (shASCL1), where a ASCL1 gene has been silenced to assess its impact on global gene expression, with two biological replicates for each group. These profiles are derived from human cell to investigate differential gene expression. This dataset is well-suited for the objective of evaluating the effects of normalization methods on RNA-seq differential expression analysis due to its clear experimental design and reproducibility. The dataset comprises raw read counts for 63,681 genes, represented as Ensembl gene IDs. After filtering out lowly expressed genes and technical artifacts, 15,697 genes with sufficient expression remained for analysis. 


 <img width="600" height="458" alt="image" src="https://github.com/user-attachments/assets/f7d714f3-0327-453b-97f6-3677927fe62b" />
## Figure 1: 
The density of log-CPM values for raw data (A) and filtered data (B) are shown for the  control samples (shGFP_1 and shGFP_2 ) and knockdown samples (shASCL1_1 and shGFP_2).

 <img width="623" height="476" alt="image" src="https://github.com/user-attachments/assets/8fd84ddc-58f7-45a5-8a31-c6cd5145adb4" />
## Figure 2: 
Boxplots of log-CPM values showing expression distributions for unnormalized data (Non) and TMM and Median normalization methods for the  control samples (shGFP_1 and shGFP_2 ) and knockdown samples (shASCL1_1 and shGFP_2).

<img width="932" height="537" alt="image" src="https://github.com/user-attachments/assets/992f8449-09ec-489a-88f0-3e1d5eaf1d5a" />
## Figure 3: 
Multidimensional scaling (MDS) plots of log-CPM values, illustrating the effect of normalization on separation and clustering of samples based on their gene expression profiles. (A, C) show MDS plots over dimensions 1 and 2, with samples colored and labeled by sample groups. (B, D) MDS plots over dimensions 2 and 3, with samples colored and labeled by biological replicates. 

<img width="881" height="508" alt="image" src="https://github.com/user-attachments/assets/e846fa16-820e-4811-a2a0-48247d21194a" />
## Figure 4: 
The left panels show Means (x-axis) and variances (y-axis) trend of each gene before applying voom precision weights, the right panels display differential expression plots for TMM (top) and Median (bottom) normalization methods.

<img width="898" height="494" alt="image" src="https://github.com/user-attachments/assets/d9c4b246-35bc-4611-b26a-6041b6bf2af0" />
## Figure 5:
Mean-Difference (MD) plots showing differential expression results for TMM (left) and Median (right) normalization methods.

<img width="905" height="522" alt="image" src="https://github.com/user-attachments/assets/7eb7395c-9a52-47e8-bf27-29def0da5f08" />
## Figure 6: 
Heatmap log-CPM values for top 50 genes DE in the  control samples (shGFP_1 and shGFP_2 ) and knockdown samples (shASCL1_1 and shGFP_2).
