# =============================================================================
# SOP: RNA-seq Differential Gene Expression Analysis using limma-voom
# =============================================================================
# Dataset:    GSE282051
# Title:      Changes in Gene Expression in Healthcare Workers During Night 
#             Shifts: Implications for Immune Response and Health Risks
# Organism:   Homo sapiens
# Tissue:     Whole blood
# Design:     Paired - 4 doctors sampled before and after a night shift
# Platform:   Illumina NovaSeq 6000
# Author:     
# Date:       
# =============================================================================
 
 
# =============================================================================
# SECTION 1: LOAD PACKAGES
# =============================================================================
 
library(limma)          # linear modeling and differential expression
library(edgeR)          # DGEList, normalization, filtering
library(BiocParallel)   # parallel processing
library(tidyverse)      # data manipulation
library(AnnotationDbi)  # gene annotation
library(org.Hs.eg.db)   # human gene annotation database
library(RColorBrewer)   # color palettes
library(gplots)         # heatmap.2
library(ggplot2)        # visualization


# =============================================================================
# SECTION 2: LOAD COUNT MATRIX
# =============================================================================

# Set working directory to folder containing count files
#setwd("path/to/your/files")  # replace with your actual path
setwd("C:/Users/nicki/xRStudiox/DEG SOP/SOP/SaimaSP2026")

# List all raw count files
files <- list.files(pattern = "*_Raw_Count.txt.gz")
print(files)  # verify order before loading
 
# Load and combine all 8 files into one count matrix
counts_list <- lapply(files, function(f) {
  df <- read.delim(gzfile(f), header = TRUE)
  rownames(df) <- make.unique(as.character(df$Gene_symbol))  # handle duplicates
  df <- df[, 2, drop = FALSE]  # keep only count column
  return(df)
})
 
counts <- do.call(cbind, counts_list)

# Rename columns to sample names
# NOTE: verify file order matches sample order before renaming
colnames(counts) <- c("before_1", "after_1", "before_2", "after_2",
                      "before_3", "after_3", "before_4", "after_4")
                      
# Round fractional counts (featureCounts -M -O --fraction option)
counts <- round(counts)
 
# Verify count matrix
dim(counts)       # should be 26364 x 8
head(counts)
all(counts >= 0)  # should be TRUE

# =============================================================================
# SECTION 3: LOAD SAMPLE METADATA
# =============================================================================
 
# Create metadata data frame
metadata <- data.frame(
  sample = c("before_1", "after_1", "before_2", "after_2",
             "before_3", "after_3", "before_4", "after_4"),
  condition = factor(c("before", "after", "before", "after",
                       "before", "after", "before", "after")),
  subject = factor(c("1", "1", "2", "2", "3", "3", "4", "4"))
)
rownames(metadata) <- metadata$sample
 
# Verify row names match column names of count matrix
all(rownames(metadata) == colnames(counts))  # should be TRUE

# Check metadata structure
str(metadata)
summary(metadata)

# =============================================================================
# SECTION 4: INITIAL DATA INSPECTION
# =============================================================================
 
dim(counts)
head(counts)
summary(counts)
 
dim(metadata)
head(metadata)
summary(metadata)


# =============================================================================
# SECTION 5: INITIAL QC PLOTS
# =============================================================================
 
# --- Library Size Distribution ---
barplot(colSums(counts),
        names.arg = colnames(counts),
        main = "Library Sizes",
        ylab = "Total Counts",
        xlab = "Sample",
        col = c(rep(c("steelblue", "tomato"), 4)),
        las = 2)
legend("topright", legend = c("Before", "After"),
       fill = c("steelblue", "tomato"))
       
# --- Log2 Count Distribution (Boxplot) ---
log_counts <- log2(counts + 1)
 
boxplot(log_counts,
        main = "Log2 Count Distribution per Sample",
        ylab = "Log2(counts + 1)",
        xlab = "Sample",
        col = c(rep(c("steelblue", "tomato"), 4)),
        las = 2)
        
# --- Density Plot ---
plot(density(log_counts[, 1]),
     main = "Density of Log2 Counts",
     xlab = "Log2(counts + 1)",
     ylim = c(0, 0.25),
     col = "steelblue")
     
for (i in 2:ncol(log_counts)) {
  lines(density(log_counts[, i]),
        col = ifelse(grepl("before", colnames(log_counts)[i]), 
                     "steelblue", "tomato"))
}
 
legend("topright", legend = c("Before", "After"),
       col = c("steelblue", "tomato"), lty = 1)
 
# --- Check library sizes and zero counts ---
colSums(counts)         # total reads per sample
colSums(counts == 0)    # number of zero count genes per sample


# =============================================================================
# SECTION 6: SPECIES CHECK & GENE ID VALIDATION
# =============================================================================
# Species: Homo sapiens -> org.Hs.eg.db
# Gene IDs: Already HGNC symbols (no conversion needed for this dataset)
 
# Validate gene symbols against official HGNC database
hgnc_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(counts),
                       column = "SYMBOL",
                       keytype = "SYMBOL",
                       multiVals = "first")
                       
# Check for NAs
sum(!is.na(hgnc_symbols))  # number mapped
sum(is.na(hgnc_symbols))   # number failed
 
# Remove NAs if any (not needed for this dataset - 0 NAs)
# counts <- counts[!is.na(hgnc_symbols), ]
# hgnc_symbols <- hgnc_symbols[!is.na(hgnc_symbols)]
 
# Confirm symbols match rownames
all(hgnc_symbols == rownames(counts))  # should be TRUE


# =============================================================================
# SECTION 7: CREATE DGEList OBJECT
# =============================================================================
 
dge <- DGEList(counts = counts,
               samples = metadata,
               group = metadata$condition)
 
# Verify
dim(dge)        # should be 26364 x 8
dge$samples     # check sample info and library sizes


# =============================================================================
# SECTION 8: FILTER LOWLY EXPRESSED GENES
# =============================================================================
 
# Filter using edgeR's recommended method
keep <- filterByExpr(dge, group = metadata$condition)
table(keep)  # see how many genes kept vs removed
 
# Apply filter
dge <- dge[keep, , keep.lib.sizes = FALSE]
 
# Verify dimensions after filtering
dim(dge)  # should be 16338 x 8


# =============================================================================
# SECTION 9: TMM NORMALIZATION
# =============================================================================
# TMM = Trimmed Mean of M-values
# Corrects for compositional bias between samples
 
dge <- calcNormFactors(dge, method = "TMM")
 
# Check normalization factors (should be close to 1)
dge$samples$norm.factors


# =============================================================================
# SECTION 10: CREATE DESIGN MATRIX
# =============================================================================
# Paired design: blocking by subject to account for individual differences
# Formula: ~ subject + condition
# Reference: subject 1, condition = after
 
design2 <- model.matrix(~ 0 + condition + subject, data = metadata)
 
# Clean up column names
colnames(design2) <- gsub("condition", "", colnames(design2))
colnames(design2)
 
# Verify design matrix
dim(design2)   # should be 8 x 5
design2


# =============================================================================
# SECTION 11: VOOM TRANSFORMATION
# =============================================================================
# Converts counts to log2 CPM and estimates mean-variance relationship
# Assigns precision weights to each observation
 
v <- voom(dge, design2, plot = TRUE)
 
# Verify voom output
dim(v$E)   # should be 16338 x 8
head(v$E)  # log2 CPM values


# =============================================================================
# SECTION 12: SAMPLE-LEVEL QC (POST-NORMALIZATION)
# =============================================================================
 
# --- MDS Plot ---
plotMDS(v,
        col = ifelse(metadata$condition == "before", "steelblue", "tomato"),
        pch = as.numeric(metadata$subject),
        main = "MDS Plot")
 
legend("topright",
       legend = c("Before", "After"),
       col = c("steelblue", "tomato"),
       pch = 1)
       
# --- Hierarchical Clustering ---
dist_mat <- dist(t(v$E))
hclust_res <- hclust(dist_mat)
 
plot(hclust_res,
     main = "Hierarchical Clustering of Samples",
     xlab = "Samples",
     sub = "",
     cex = 0.9)
     
# --- Sample-to-Sample Correlation Heatmap ---
cor_mat <- cor(v$E, method = "pearson")
 
heatmap.2(cor_mat,
          trace = "none",
          col = colorRampPalette(brewer.pal(9, "Blues"))(100),
          margins = c(8, 8),
          main = "Sample-to-Sample Correlation",
          key.title = "Correlation",
          dendrogram = "both",
          ColSideColors = ifelse(metadata$condition == "before",
                                 "steelblue", "tomato"),
          RowSideColors = ifelse(metadata$condition == "before",
                                 "steelblue", "tomato"))
 
legend("topright", legend = c("Before", "After"),
       fill = c("steelblue", "tomato"), cex = 0.8)     

# IMPORTANT: Always reset graphics after heatmap.2 to avoid corrupted 
# graphics state when running ggplot afterwards
graphics.off()

# --- Cook's Distance ---
# NOTE: Cook's distance was uninformative for this dataset due to model
# saturation with n=4 paired samples. Sample QC was assessed via MDS,
# hierarchical clustering, and correlation heatmap instead.
# after_4 was consistently flagged as different across all QC metrics
# but retained due to small sample size (n=4).


# =============================================================================
# SECTION 13: FIT LINEAR MODEL & DEFINE CONTRASTS
# =============================================================================
 
# Fit linear model
fit2 <- lmFit(v, design2)

# Define contrasts (after vs before night shift)
contrast_matrix <- makeContrasts(
  AfterVsBefore = after - before,
  levels = design2
)
 
contrast_matrix  # verify contrast

# Fit contrasts
fit2 <- contrasts.fit(fit2, contrast_matrix)
 
# Apply empirical Bayes smoothing
# Shrinks gene-wise variances toward a common value - improves power
fit2 <- eBayes(fit2)
 
# Check fit object
names(fit2)
dim(fit2$coefficients)  # should be 16338 x 1


# =============================================================================
# SECTION 14: EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# =============================================================================
 
results2 <- topTable(fit2,
                     coef = "AfterVsBefore",
                     number = Inf,
                     adjust.method = "BH",
                     sort.by = "P")
 
# Check results
head(results2)
dim(results2)

# Summary of significant genes
# NOTE: No genes pass adj.P.Val < 0.05 due to small sample size (n=4)
sum(results2$adj.P.Val < 0.05)   # adjusted p-value threshold
sum(results2$P.Value < 0.05)     # raw p-value threshold
 
# Genes passing paper's threshold (P < 0.05 & |logFC| > 0.5)
sig_genes2 <- results2[results2$P.Value < 0.05 &
                         abs(results2$logFC) > 0.5, ]
nrow(sig_genes2)
sum(sig_genes2$logFC > 0)  # upregulated after shift
sum(sig_genes2$logFC < 0)  # downregulated after shift


# =============================================================================
# SECTION 15: HGNC ANNOTATION & JOIN
# =============================================================================
 
# Download HGNC complete set file
download.file(
  url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
  destfile = "hgnc_complete_set.txt",
  mode = "wb"
)
 
# Load HGNC complete set
hgnc_complete_set_df <- read_tsv("hgnc_complete_set.txt")

# Rename symbol column
colnames(hgnc_complete_set_df)[colnames(hgnc_complete_set_df) == "symbol"] <- "HGNC_Symbol"

# Convert results to tibble with gene symbols as column
res_df <- results2 %>%
  as_tibble(rownames = "HGNC_Symbol")
  
# Perform inner join with HGNC database
joined_results <- res_df %>%
  inner_join(dplyr::select(hgnc_complete_set_df, HGNC_Symbol),
             by = "HGNC_Symbol")
             
# Check how many genes retained
dim(joined_results)
nrow(results2) - nrow(joined_results)  # genes lost in join


# =============================================================================
# SECTION 16: HANDLE DUPLICATES & FINALIZE RESULTS TABLE
# =============================================================================
 
# Handle duplicate gene symbols (keep lowest adj.P.Val)
# Remove NAs and select final columns
final_results <- joined_results %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::distinct(HGNC_Symbol, .keep_all = TRUE) %>%
  dplyr::select(HGNC_Symbol, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(raw_P.value = P.Value,
                adj_P.value = adj.P.Val)
                
# Verify final table
colnames(final_results)
dim(final_results)
head(final_results)
sum(is.na(final_results$adj_P.value))  # should be 0
 
# Save final results table
write.csv(final_results,
          "GSE282051_DGE_results.csv",
          row.names = FALSE)
          

# =============================================================================
# SECTION 17: VISUALIZATION
# =============================================================================
 
# --- Volcano Plot ---
final_results <- final_results %>%
  dplyr::mutate(significance = case_when(
    raw_P.value < 0.05 & logFC > 0.5  ~ "Up after shift",
    raw_P.value < 0.05 & logFC < -0.5 ~ "Down after shift",
    TRUE ~ "Not significant"
  ))
  
ggplot(final_results, aes(x = logFC, y = -log10(raw_P.value),
                           color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up after shift" = "tomato",
                                "Down after shift" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "black") +
  labs(title = "Volcano Plot - Night Shift Gene Expression",
       subtitle = "After vs Before Night Shift",
       x = "Log2 Fold Change",
       y = "-Log10(P-value)",
       color = "Direction") +
  theme_bw() +
  theme(legend.position = "right")
 
ggsave("GSE282051_volcano_plot.png", width = 8, height = 6, dpi = 300)  

# --- Heatmap of Top 50 DE Genes ---
 
# Get top 50 significant genes
top_genes <- final_results %>%
  dplyr::filter(raw_P.value < 0.05 & abs(logFC) > 0.5) %>%
  dplyr::arrange(raw_P.value) %>%
  head(50)
  
# Extract and scale expression data
top_expr <- v$E[top_genes$HGNC_Symbol, ]
top_expr_scaled <- t(scale(t(top_expr)))
 
# Order columns by condition
col_order <- c("before_1", "before_2", "before_3", "before_4",
               "after_1", "after_2", "after_3", "after_4")
top_expr_ordered <- top_expr_scaled[, col_order]
 
# Average before and after samples
before_avg <- rowMeans(top_expr_ordered[, grepl("before", colnames(top_expr_ordered))])
after_avg <- rowMeans(top_expr_ordered[, grepl("after", colnames(top_expr_ordered))])
 
top_expr_avg <- data.frame(Before = before_avg, After = after_avg)
 
# Prepare for ggplot
heatmap_data_avg <- top_expr_avg %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(cols = -Gene,
                      names_to = "Condition",
                      values_to = "Expression")
 
# Plot heatmap
ggplot(heatmap_data_avg, aes(x = Condition, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "#053061",
                       mid = "white",
                       high = "#67001F",
                       midpoint = 0,
                       limits = c(-2, 2),
                       name = "Z-score") +
  labs(title = "Top 50 DE Genes - Night Shift",
       x = "Condition", y = "Gene") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(10, 10, 10, 60),
        panel.border = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank()) +
  scale_y_discrete(expand = expansion(add = 0.5))
 
ggsave("GSE282051_heatmap.png", width = 6, height = 14, dpi = 300)

# =============================================================================
# SECTION 18: SAVE WORKSPACE
# =============================================================================
 
save.image("GSE282051_analysis.RData")
 
# To reload in a future session:
# load("GSE282051_analysis.RData")
# Then reload all libraries listed in Section 1
 
# =============================================================================
# END OF SOP
# =============================================================================