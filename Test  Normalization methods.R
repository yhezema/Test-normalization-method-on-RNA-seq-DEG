#----------------------------------------------------------------------------------------------#
#Impact of TMM and Median normalization method on RNA‑seq differential gene expression analysis 
# Author: "Yasmine Hezema"
# Date: "12-12-2024"
#---------------------------------=------------------------------------------------------------#

# Load Required Libraries -----
library(limma) 
library(Glimma)
library(edgeR)
library(org.Hs.eg.db)
library(R.utils)
library(gplots)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
 
# Code Section 1 – Data Acquisition, Exploration, Filtering, and Quality Control ----

## Data Acquisition ----
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132698&format=file"
utils::download.file(url, destfile = "GSE132698_RAW.tar", mode = "wb") 
utils::untar("GSE132698_RAW.tar", exdir = ".")

files <- c("GSM3887869_EXP2_shGFP_1.txt", 
           "GSM3887870_EXP2_shGFP_2.txt", 
           "GSM3887871_EXP2_shASCL1_1.txt", 
           "GSM3887872_EXP2_shASCL1_2.txt")
# Decompress and clean files
for (file in paste0(files, ".gz")) {
  R.utils::gunzip(file, overwrite = TRUE)
}
# The data are derived from Wang et al. (2019), and include two treatment conditions; (control and knockdown) with two biological replicates each.

## Data Exploration ----
# Preview the structure of one file
read.delim(files[1], nrow = 5, header = FALSE)
# The file contains two columns the first is the EntrezIDs and the second is the count with no header.

# Combine raw count data into a matrix
x <- readDGE(files, columns = c(1, 2), header = FALSE)
dim(x)  # Dimensions of the DGEList object
# The resulting DGEList-object contains a matrix of 63681 rows (unique Entrez gene identifiers (IDs)) and four columns (samples). However, Meta tags were detected in the files e.g."__no_feature". Therefore, the files need to be cleaned.

# Clean raw data by removing meta tags (e.g., "__no_feature")
cleaned_files <- c()
for (file in files) {
  # Read the file
  data <- read.delim(file, header = FALSE)
  # Filter out rows that have meta tags (starting with "__")
  cleaned_data <- data[!grepl("^__", data$V1), ]
  # Add column names to the cleaned data
  colnames(cleaned_data) <- c("IDs", "Count")
  # Create a new file name for the cleaned file
  cleaned_file <- paste0("cleaned_", file)
  # Save the cleaned file with column names
  write.table(cleaned_data, file = cleaned_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  # Add the cleaned file name to the vector
  cleaned_files <- c(cleaned_files, cleaned_file)
}
# Confirm cleaned files and reload data
x_cleaned <- readDGE(cleaned_files, columns = c(1, 2), header = TRUE)
# Update sample names for clarity
samplenames <- sub(".*EXP2_", "", colnames(x_cleaned))
colnames(x_cleaned) <- samplenames # Rename Column Names
# Assign group information to samples
x_cleaned$samples$group <- factor(c("Control", "Control", "Knockdown", "Knockdown"))

# Update sample names to include only the part after "EXP2_"
samplenames <- sub(".*EXP2_", "", colnames(x_cleaned))
colnames(x_cleaned) <- samplenames
samplenames # check the names

## Annotation Preparation ----
# Retrieve Ensembl IDs from the rownames
geneid <- rownames(x_cleaned$counts)
head(geneid) # Check gene id's 

# From SRA link "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA548772&o=acc_s%3Aa" I have tried to find the run information, since there is only one run/sample and no indication of multiple lanes or batches, I did not specify lane information in the analysis.

# The IDs in the geneid are Ensembl Gene IDs (e.g., ENSG00000000003) which are not valid ENTREZID identifiers so we need to map the Ensembl Gene IDs to Entrez Gene IDs before querying org.Hs.eg.db.

# Map Ensembl IDs to Entrez IDs and SYMBOLs
ensembl_to_entrez <- AnnotationDbi::select(  x = org.Hs.eg.db, keys = geneid,  columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
head(ensembl_to_entrez)
# Remove duplicates and inspect results
ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez$ENSEMBL), ]

# Map Entrez IDs to chromosomes
txdb_chrom <- AnnotationDbi::select(x = TxDb.Hsapiens.UCSC.hg38.knownGene,  keys = ensembl_to_entrez$ENTREZID,  columns = c("TXCHROM"), keytype = "GENEID")
head(txdb_chrom)

# Remove duplicate entries
txdb_chrom <- txdb_chrom[!duplicated(txdb_chrom$GENEID), ]

# Combine mappings and remove unnecessary columns
annotations <- merge(ensembl_to_entrez, txdb_chrom, by.x = "ENTREZID", by.y = "GENEID", all.x = TRUE)
head(annotations)
# Drop the ENSEMBL column
annotations <- annotations[, !colnames(annotations) %in% "ENSEMBL"] 
head(annotations)

# Attach annotations to the DGEList object
x_cleaned$genes <- annotations
head(x_cleaned$genes)
x_cleaned

## Data processing ----
# Calculate counts per million (CPM) and log-CPM for raw data
cpm_raw <- cpm(x_cleaned)
lcpm_raw <- cpm(x_cleaned, log = TRUE)

# Summarize the raw CPM and log-CPM distributions
summary(cpm_raw)
summary(lcpm_raw)

# Library size summary (mean and median in millions)
L <- mean(x_cleaned$samples$lib.size) * 1e-6
M <- median(x_cleaned$samples$lib.size) * 1e-6
c(L, M)
# The library sizes (mean = 31.13M reads, median = 31.00M reads) indicate balanced sequencing depths across samples.

# Identify the number of non-expressed genes (zero counts across all samples)
table(rowSums(x_cleaned$counts == 0) == 4)
# The dataset contains 37,054 non-expressed genes (zero counts in all samples) and 26,623 expressed genes.

# Define group variable from x_cleaned samples
group <- x_cleaned$samples$group
group

### Data Filtering ----
# Filter out lowly expressed genes
keep.exprs <- filterByExpr(x_cleaned, group = group, min.count = 10)
filtered_x_cleaned <- x_cleaned[keep.exprs, , keep.lib.sizes = FALSE]

# Check the dimensions after filtering
dim(filtered_x_cleaned)
table(keep.exprs)
# Of 63,681 total genes, 15,697 genes with sufficient expression remain after filtering, while 47,980 genes with low counts are excluded.

# Calculate log-CPM for filtered data
lcpm_filtered <- cpm(filtered_x_cleaned, log = TRUE)

# Define cutoff for filtering
lcpm_cutoff <- log2(10 / M + 2 / L)

# Number of samples and assign colors for plotting
nsamples <- ncol(x_cleaned$counts)
col <- brewer.pal(nsamples, "Paired")

# Define plotting parameters
visualize_densities <- function(lcpm_data, title, cutoff, colors, samples, ylim = c(0, 0.26)) {
  plot(density(lcpm_data[, 1]), col = colors[1], lwd = 2, ylim = ylim, las = 2, main = "", xlab = "Log-cpm")
  title(main = title)
  abline(v = cutoff, lty = 3)
  for (i in 2:ncol(lcpm_data)) {
    lines(density(lcpm_data[, i]), col = colors[i], lwd = 2)
  }
  legend("topright", legend = samples, text.col = colors, bty = "n")
}
# Set up plotting area for side-by-side plots
par(mfrow = c(1, 2))
# Example usage for `lcpm_raw`
visualize_densities(lcpm_raw, "A. Raw Data", lcpm_cutoff, col, samplenames)
visualize_densities(lcpm_filtered, "B. Filtered Data", lcpm_cutoff, col, samplenames)

### Data Normalization ----
# Normalize filtered data using TMM (Trimmed Mean of M-values) method
filtered_x_cleaned <- calcNormFactors(filtered_x_cleaned, method = "TMM")

# Display normalization factors
filtered_x_cleaned$samples$norm.factors
# The calculated normalization factors (e.g., 1.0328, 1.0273, 0.9752, 0.9666); suggest that only minor scaling adjustments were necessary. This confirms that the library sizes across samples were already well-balanced.

# Log-CPM for CPM-normalized data
lcpm_tmm <- cpm(filtered_x_cleaned, log = TRUE)

# Normalize filtered data using Median method
median_factors <- apply(filtered_x_cleaned$counts, 2, median)
median_normalized_counts <- sweep(filtered_x_cleaned$counts, 2, median_factors, FUN = "/")

# Create a new DGEList object for median-normalized data
filtered_x_cleaned_median <- filtered_x_cleaned
filtered_x_cleaned_median$counts <- median_normalized_counts

# Log-CPM for Median-normalized data
lcpm_median <- log2(median_normalized_counts + 1)

### Create Boxplots Comparing Normalization Methods -----
# Combine data into a long-format dataframe for easier plotting
boxplot_data <- data.frame(
  Method = rep(c("Non", "TMM", "Median"), each = nrow(lcpm_filtered) * ncol(lcpm_raw)),
  Sample = rep(rep(colnames(filtered_x_cleaned$counts), each = nrow(lcpm_filtered)), 3),
  LogCPM = c(as.vector(lcpm_filtered), as.vector(lcpm_tmm), as.vector(lcpm_median))
)

ggplot(boxplot_data, aes(x = Sample, y = LogCPM, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  labs(
    title = "Comparison of Normalization Methods",
    x = " ",
    y = "Log-CPM",
    fill = "Methods"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )


### Multidimensional Scaling (MDS)----
# The MDS plot shows the leading log-fold change dimensions for the most variable genes.

# Define a function to plot MDS
plot_mds <- function(data, title, dim = c(1, 2), xlim = c(-1.2, 1.2), ylim = c(-0.8, 0.8), col, labels) {
  plotMDS(data, labels = labels, col = col, dim = dim, cex = 1.5, xlim = xlim, ylim = ylim)
  title(main = title, cex.main = 2)
}

# Add biological replicates to sample metadata
filtered_x_cleaned$samples$replicate <- c("Rep1", "Rep2", "Rep1", "Rep2")

# Number of samples and colors
nsamples <- ncol(filtered_x_cleaned$counts)
col <- brewer.pal(nsamples, "Paired")
labels <- colnames(filtered_x_cleaned$counts)

# Set up plot area for six plots
par(mfrow = c(2, 2), mar = c(5, 6, 4, 2))

# Plot MDS for TMM Normalized Data
plot_mds(lcpm_tmm, "A. Sample Groups (TMM)", col = col, labels = labels)
plot_mds(lcpm_tmm, "B. Biological Replicates (TMM)", dim = c(2, 3), col = col, labels = labels)

# Plot MDS for Median Normalized Data
plot_mds(lcpm_median, "C. Sample Groups (Median)", col = col, labels = labels)
plot_mds(lcpm_median, "D. Biological Replicates (Median)", dim = c(2, 3), col = col, labels = labels)

# TMM normalization captured variance of 67 % that likely maintains more biological signal and variation, which is crucial for downstream analyses of identifying differentially expressed genes. While Median normalization captured variance of 93% that might simplify the data more aggressively, which might reduce noise but risks losing important biological distinctions. In conclusion, the MDS plots suggested that median normalization may introduce more uniformity across samples by dampening biological variability, whereas TMM normalization better retains the biological signal while addressing technical differences.

# Code Section 2 – Main Analysis ------
## Examining Differentially Expressed Genes ----
# Function for Differential Expression Analysis
run_limma_pipeline <- function(norm_data, group, contrast, design_matrix) {
  # Voom transformation
  voom_data <- voom(norm_data, design_matrix, plot = TRUE)
  fit <- lmFit(voom_data, design_matrix) # Linear model fitting
  fit <- contrasts.fit(fit, contrast) # Apply contrasts
  efit <- eBayes(fit) # # Empirical Bayes moderation
  # Return the efit object for further analysis
  return(efit)
}

# Define the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

# Define the contrast matrix
contrast_matrix <- makeContrasts(ControlvsKnockdown = Control - Knockdown, levels = colnames(design))

# Function to analyze and plot results
analyze_and_plot <- function(efit, method_name) {
  # Summarize results
  de_results <- summary(decideTests(efit))
  print(paste(method_name, "Normalization Results:"))
  print(de_results)
  
  # Create a mean-difference (MD) plot with method name in the title
  plotMD(efit, column = 1, xlim = c(-8, 13),
    main = paste("Differential Expression Plot:", method_name, "Normalization"),
    ylab = "Log Fold Change", legend = TRUE)
}

# Set up plot area for two plots
par(mfrow = c(2, 2))

# Run the pipeline and plot for TMM normalization
efit_tmm <- run_limma_pipeline(filtered_x_cleaned, group, contrast_matrix, design)
analyze_and_plot(efit_tmm, "TMM")

# Run the pipeline and plot for Median normalization
efit_median <- run_limma_pipeline(filtered_x_cleaned_median, group, contrast_matrix, design)
analyze_and_plot(efit_median, "Median")
# TMM and median methods identified, in Knockdown vs. Control, 3738 and 1965 significantly upregulated , 3973 and 5457 significantly downregulated in and 7986 and 8275 genes not significantly differentially expressed, respectively.

## Applying a Criterion for Log-Fold Changes ----
# Function to apply log-fold change criterion
apply_lfc_threshold <- function(efit, lfc_threshold) {
  # Apply treat function with the log-fold change threshold
  tfit <- treat(efit, lfc = lfc_threshold)
  dt <- decideTests(tfit) # Decide tests
  
  # Summarize results
  logfc_results <- summary(dt)
  return(list(tfit = tfit, dt = dt, logfc_results = logfc_results))
}

# Apply the criterion for TMM-Normalized Data
lfc_threshold <- 1
results_tmm <- apply_lfc_threshold(efit_tmm, lfc_threshold)
print("Log-Fold Change Results for TMM Normalization:")
print(results_tmm$logfc_results)

# Apply the criterion for Median-Normalized Data
results_median <- apply_lfc_threshold(efit_median, lfc_threshold)
print("Log-Fold Change Results for Median Normalization:")
print(results_median$logfc_results)
# However, Median showed higher differential expression results based on statistical testing (P value), after applying log fold filter, TMM normalization resulted in more differentially expressed genes in both upregulated (TMM 73  vs Meduian 9) and downregulated (TMM 104  vs Median 32) compared to the Median normalization method.

## Identify DE genes for a given method----
# Define a function to identify DE genes for a given method
identify_de_genes <- function(tfit, dt) {
  # Identify DE genes for the single contrast
  de.common <- which(dt[, 1] != 0)
  num_de_genes <- length(de.common) # Number of DE genes
  # Retrieve the gene symbols of the DE genes
  de_gene_symbols <- head(tfit$genes$SYMBOL[de.common], n = 20)
  # Return results as a list
  list(num_de_genes = num_de_genes, de_gene_symbols = de_gene_symbols)
}

# Apply the function to TMM normalization
de_tmm <- identify_de_genes(results_tmm$tfit, results_tmm$dt)
print(paste("Number of DE genes (TMM):", de_tmm$num_de_genes))
print("Top DE genes (TMM):")
print(de_tmm$de_gene_symbols)

# Apply the function to Median normalization
de_median <- identify_de_genes(results_median$tfit, results_median$dt)
print(paste("Number of DE genes (Median):", de_median$num_de_genes))
print("Top DE genes (Median):")
print(de_median$de_gene_symbols)
# Number of DE genes using TMM are 177 and 41 genes by using Median and the DEGs in each methods are different.
#If the study focuses on detecting large biological changes, LFC results may be more informative. However, if the goal is to identify reliable genes for further validation, DE results are more robust and statistically defensible.

## Examining Individual Gene Results and Plotting Differential Gene Expression ----
# Function to analyze and plot differential gene expression for a given normalization method
analyze_and_plot <- function(tfit, dt, method_name) {
  # Extract top significant results
  gene_results <- topTreat(tfit, coef = 1, n = Inf)
  
  # Display the top significant results
  print(paste("Top significant results for", method_name, "normalization:"))
  print(head(gene_results))
  
  # Plot Mean-Difference (MD) Plot
  plotMD(tfit, column = 1,
    status = dt[, 1],  
    xlim = c(-8, 13), main = paste("Differential Expression Plot:", method_name, "Normalization"),
    ylab = "Log Fold Change", legend = TRUE)
}

# Graphical representations of differential expression results
# Applying a Criterion for Log-Fold Changes
apply_lfc_threshold <- function(efit, lfc_threshold) {
  tfit <- treat(efit, lfc = lfc_threshold)
  dt <- decideTests(tfit)
  return(list(tfit = tfit, dt = dt))
}

# Define the log-fold change threshold
lfc_threshold <- 1

# Apply the criterion for TMM-normalized data
tmm_results <- apply_lfc_threshold(efit_tmm, lfc_threshold)
tfit_tmm <- tmm_results$tfit
dt_tmm <- tmm_results$dt

# Apply the criterion for Median-normalized data
median_results <- apply_lfc_threshold(efit_median, lfc_threshold)
tfit_median <- median_results$tfit
dt_median <- median_results$dt

# Set up plot area for two plots
par(mfrow = c(1, 2))

# Apply the function for TMM-normalized data
analyze_and_plot(tfit_tmm, dt_tmm, "TMM")

# Apply the function for Median-normalized data
analyze_and_plot(tfit_median, dt_median, "Median")

## Heatmap ------
# Function to update row names to gene symbols
update_row_names <- function(data_matrix, gene_symbols) {
  rownames(data_matrix) <- gene_symbols
  return(data_matrix)
}

# Check and update row names with gene symbols for both TMM and Median
if (sum(is.na(filtered_x_cleaned$genes$SYMBOL)) > 0) {
  filtered_x_cleaned$genes$SYMBOL[is.na(filtered_x_cleaned$genes$SYMBOL)] <- rownames(filtered_x_cleaned$genes)[is.na(filtered_x_cleaned$genes$SYMBOL)]
}
if (sum(is.na(filtered_x_cleaned_median$genes$SYMBOL)) > 0) {
  filtered_x_cleaned_median$genes$SYMBOL[is.na(filtered_x_cleaned_median$genes$SYMBOL)] <- rownames(filtered_x_cleaned_median$genes)[is.na(filtered_x_cleaned_median$genes$SYMBOL)]
}

# Select the top 50 genes by variance for each normalization method
top_genes_tmm <- head(order(apply(lcpm_tmm, 1, var), decreasing = TRUE), 50)
top_genes_median <- head(order(apply(lcpm_median, 1, var), decreasing = TRUE), 50)

# Update row names to gene symbols
lcpm_tmm <- update_row_names(lcpm_tmm, filtered_x_cleaned$genes$SYMBOL)
lcpm_median <- update_row_names(lcpm_median, filtered_x_cleaned_median$genes$SYMBOL)

# Define a color palette for the heatmap
mycol <- colorRampPalette(c("blue", "white", "red"))(100)

# Reorder samples so that treatments are grouped together
reordered_samples <- c("shGFP_1", "shGFP_2", "shASCL1_1", "shASCL1_2")
lcpm_tmm <- lcpm_tmm[, reordered_samples]
lcpm_median <- lcpm_median[, reordered_samples]

# Create heatmaps for TMM and Median normalized data
ht1 <- Heatmap(
  lcpm_tmm[top_genes_tmm, ],
  name = "Expression (TMM)",
  col = mycol,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),  
  column_names_gp = gpar(fontsize = 10),  
  column_title = "TMM Normalization",
  column_title_gp = gpar(fontface = "bold"),
  heatmap_legend_param = list(title = "Expression (TMM)")
)

ht2 <- Heatmap(
  lcpm_median[top_genes_median, ],
  name = "Expression (Median)",
  col = mycol,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),  # Adjust font size for gene names
  column_names_gp = gpar(fontsize = 10),  # Adjust font size for sample names
  column_title = "Median Normalization",
  column_title_gp = gpar(fontface = "bold"),
  heatmap_legend_param = list(title = "Expression (Median)")
)

# Combine heatmaps into a single figure with one legend
ht_list <- ht1 + ht2

# Draw the combined heatmap
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)

## Reproducibility Analysis ----
# Define `logfc_comparison` if it is not defined
logfc_comparison <- data.frame(
  Gene = rownames(efit_tmm$coefficients),
  LogFC_TMM = efit_tmm$coefficients[, "ControlvsKnockdown"],
  LogFC_Median = efit_median$coefficients[, "ControlvsKnockdown"]
)

# Function to capture differential expression genes and runtime
analyze_de_genes <- function(fit, method_name) {
  # Measure runtime
  runtime <- system.time({
    de_genes <- rownames(logfc_comparison)[decideTests(fit)[, "ControlvsKnockdown"] != 0]
  })
  
  # Return results and runtime
  list(method = method_name, de_genes = de_genes, runtime = runtime["elapsed"])
}

# Analyze DE genes for TMM-normalized data
tmm_results <- analyze_de_genes(efit_tmm, "TMM")
# Analyze DE genes for Median-normalized data
median_results <- analyze_de_genes(efit_median, "Median")
# Calculate overlap of DE genes
common_genes <- intersect(tmm_results$de_genes, median_results$de_genes)
unique_tmm_genes <- setdiff(tmm_results$de_genes, median_results$de_genes)
unique_median_genes <- setdiff(median_results$de_genes, tmm_results$de_genes)
# Calculate percentage overlap
overlap_percentage <- length(common_genes) / length(unique(c(tmm_results$de_genes, median_results$de_genes))) * 100
# Print results
cat("Reproducibility Analysis:\n")
cat("TMM DE Genes:", length(tmm_results$de_genes), "(", tmm_results$runtime, "seconds runtime)\n") # 7711 ( 0.05 seconds runtime)
cat("Median DE Genes:", length(median_results$de_genes), "(", median_results$runtime, "seconds runtime)\n") # 7422 ( 0 seconds runtime)
cat("Common DE Genes:", length(common_genes), "\n") # 5926 
cat("Unique to TMM:", length(unique_tmm_genes), "\n") # 1785
cat("Unique to Median:", length(unique_median_genes), "\n") # 1496 
cat("Percentage overlap of DE genes:", round(overlap_percentage, 2), "%\n") # 64.36 %


