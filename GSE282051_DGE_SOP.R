library(limma)          # linear modeling and differential expression
> library(edgeR)          # DGEList, normalization, filtering
> library(BiocParallel)   # parallel processing
> library(tidyverse)      # data manipulation
── Attaching core tidyverse packages ──────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.2.0     ✔ readr     2.2.0
✔ forcats   1.0.1     ✔ stringr   1.6.0
✔ ggplot2   4.0.2     ✔ tibble    3.3.1
✔ lubridate 1.9.5     ✔ tidyr     1.3.2
✔ purrr     1.2.1     
── Conflicts ────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package to force all conflicts to become errors
> library(AnnotationDbi)  # gene annotation
Loading required package: stats4
Loading required package: BiocGenerics
Loading required package: generics

Attaching package: ‘generics’

The following object is masked from ‘package:lubridate’:
  
  as.difftime

The following object is masked from ‘package:dplyr’:
  
  explain

The following objects are masked from ‘package:base’:
  
  as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
setequal, union


Attaching package: ‘BiocGenerics’

The following object is masked from ‘package:dplyr’:
  
  combine

The following object is masked from ‘package:limma’:
  
  plotMA

The following objects are masked from ‘package:stats’:
  
  IQR, mad, sd, var, xtabs

The following objects are masked from ‘package:base’:
  
  anyDuplicated, aperm, append, as.data.frame, basename, cbind,
colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
Reduce, rownames, sapply, saveRDS, table, tapply, unique, unsplit,
which.max, which.min

Loading required package: Biobase
Welcome to Bioconductor

Vignettes contain introductory material; view with
'browseVignettes()'. To cite Bioconductor, see 'citation("Biobase")',
and for packages 'citation("pkgname")'.

Loading required package: IRanges
Loading required package: S4Vectors

Attaching package: ‘S4Vectors’

The following objects are masked from ‘package:lubridate’:
  
  second, second<-
  
  The following objects are masked from ‘package:dplyr’:
  
  first, rename

The following object is masked from ‘package:tidyr’:
  
  expand

The following object is masked from ‘package:utils’:
  
  findMatches

The following objects are masked from ‘package:base’:
  
  expand.grid, I, unname


Attaching package: ‘IRanges’

The following object is masked from ‘package:lubridate’:
  
  %within%
  
  The following objects are masked from ‘package:dplyr’:
  
  collapse, desc, slice

The following object is masked from ‘package:purrr’:
  
  reduce

The following object is masked from ‘package:grDevices’:
  
  windows


Attaching package: ‘AnnotationDbi’

The following object is masked from ‘package:dplyr’:
  
  select
> library(org.Hs.eg.db)   # human gene annotation database

> library(RColorBrewer)   # color palettes
> library(gplots)         # heatmap.2

---------------------
  gplots 3.3.0 loaded:
  * Use citation('gplots') for citation info.
* Homepage: https://talgalili.github.io/gplots/
  * Report issues: https://github.com/talgalili/gplots/issues
* Ask questions: https://stackoverflow.com/questions/tagged/gplots
* Suppress this message with: suppressPackageStartupMessages(library(gplots))
---------------------
  
  
  Attaching package: ‘gplots’

The following object is masked from ‘package:IRanges’:
  
  space

The following object is masked from ‘package:S4Vectors’:
  
  space

The following object is masked from ‘package:stats’:
  
  lowess
> library(ggplot2)        # visualization
> getwd()
[1] "C:/Users/nicki/xRStudiox/DEG SOP/SOP"
> # Set working directory to folder containing count files
  > #setwd("path/to/your/files")  # replace with your actual path
  > setwd("C:/Users/nicki/xRStudiox/DEG SOP/SOP/SaimaSP2026")
> # List all raw count files
  > files <- list.files(pattern = "*_Raw_Count.txt.gz")
> print(files)  # verify order before loading
[1] "GSM8634930_before_1_Raw_Count.txt.gz" "GSM8634931_after_1_Raw_Count.txt.gz" 
[3] "GSM8634932_before_2_Raw_Count.txt.gz" "GSM8634933_after_2_Raw_Count.txt.gz" 
[5] "GSM8634934_before_3_Raw_Count.txt.gz" "GSM8634935_after_3_Raw_Count.txt.gz" 
[7] "GSM8634936_before_4_Raw_Count.txt.gz" "GSM8634937_after_4_Raw_Count.txt.gz" 
>  
  > # Load and combine all 8 files into one count matrix
  > counts_list <- lapply(files, function(f) {
    +   df <- read.delim(gzfile(f), header = TRUE)
    +   rownames(df) <- make.unique(as.character(df$Gene_symbol))  # handle duplicates
    +   df <- df[, 2, drop = FALSE]  # keep only count column
    +   return(df)
    + })
>  
  > counts <- do.call(cbind, counts_list)
> # Rename columns to sample names
  > # NOTE: verify file order matches sample order before renaming
  > colnames(counts) <- c("before_1", "after_1", "before_2", "after_2",
                          +                       "before_3", "after_3", "before_4", "after_4")
> # Round fractional counts (featureCounts -M -O --fraction option)
  > counts <- round(counts)
>  
  > # Verify count matrix
  > dim(counts)       # should be 26364 x 8
[1] 26364     8
> head(counts)
before_1 after_1 before_2 after_2 before_3 after_3 before_4 after_4
DDX11L1        325     241      756     986      178     154      149     115
WASH7P        2536    2330     3050    4518     2438    1878     3075    1458
MIR6859-3       54      48       73      89       54      30       74      33
MIR6859-2       54      48       73      89       54      30       74      33
MIR6859-4       54      48       73      89       54      30       74      33
MIR6859-1       54      48       73      89       54      30       74      33
> all(counts >= 0)  # should be TRUE
[1] TRUE
> # =============================================================================
> # SECTION 3: LOAD SAMPLE METADATA
  > # =============================================================================
>  
  > # Create metadata data frame
  > metadata <- data.frame(
    +   sample = c("before_1", "after_1", "before_2", "after_2",
                   +              "before_3", "after_3", "before_4", "after_4"),
    +   condition = factor(c("before", "after", "before", "after",
                             +                        "before", "after", "before", "after")),
    +   subject = factor(c("1", "1", "2", "2", "3", "3", "4", "4"))
    + )
> rownames(metadata) <- metadata$sample
>  
  > # Verify row names match column names of count matrix
  > all(rownames(metadata) == colnames(counts))  # should be TRUE
[1] TRUE
> # Check metadata structure
  > str(metadata)
'data.frame':	8 obs. of  3 variables:
  $ sample   : chr  "before_1" "after_1" "before_2" "after_2" ...
$ condition: Factor w/ 2 levels "after","before": 2 1 2 1 2 1 2 1
$ subject  : Factor w/ 4 levels "1","2","3","4": 1 1 2 2 3 3 4 4
> summary(metadata)
sample           condition subject
Length:8           after :4   1:2    
Class :character   before:4   2:2    
Mode  :character              3:2    
4:2    
> # =============================================================================
> # SECTION 4: INITIAL DATA INSPECTION
  > # =============================================================================
>  
  > dim(counts)
[1] 26364     8
> head(counts)
before_1 after_1 before_2 after_2 before_3 after_3 before_4 after_4
DDX11L1        325     241      756     986      178     154      149     115
WASH7P        2536    2330     3050    4518     2438    1878     3075    1458
MIR6859-3       54      48       73      89       54      30       74      33
MIR6859-2       54      48       73      89       54      30       74      33
MIR6859-4       54      48       73      89       54      30       74      33
MIR6859-1       54      48       73      89       54      30       74      33
> summary(counts)
before_1           after_1           before_2          after_2        
Min.   :       0   Min.   :      0   Min.   :      0   Min.   :       0  
1st Qu.:       0   1st Qu.:      0   1st Qu.:      0   1st Qu.:       0  
Median :      90   Median :     93   Median :     84   Median :     109  
Mean   :    3392   Mean   :   2806   Mean   :   3064   Mean   :    3723  
3rd Qu.:    1514   3rd Qu.:   1491   3rd Qu.:   1369   3rd Qu.:    1698  
Max.   :11965563   Max.   :5774089   Max.   :8101009   Max.   :11313581  
before_3           after_3            before_4           after_4        
Min.   :       0   Min.   :       0   Min.   :       0   Min.   :       0  
1st Qu.:       0   1st Qu.:       0   1st Qu.:       0   1st Qu.:       0  
Median :      80   Median :      69   Median :      83   Median :      38  
Mean   :    3030   Mean   :    2949   Mean   :    4134   Mean   :    2880  
3rd Qu.:    1251   3rd Qu.:    1076   3rd Qu.:    1376   3rd Qu.:     672  
Max.   :11721612   Max.   :11562323   Max.   :20375207   Max.   :15000541  
>  
  dim(metadata)
[1] 8 3
> head(metadata)
sample condition subject
before_1 before_1    before       1
after_1   after_1     after       1
before_2 before_2    before       2
after_2   after_2     after       2
before_3 before_3    before       3
after_3   after_3     after       3
> summary(metadata)
sample           condition subject
Length:8           after :4   1:2    
Class :character   before:4   2:2    
Mode  :character              3:2    
4:2    
> # =============================================================================
> # SECTION 5: INITIAL QC PLOTS
  > # =============================================================================
>  
  > # --- Library Size Distribution ---
  > barplot(colSums(counts),
            +         names.arg = colnames(counts),
            +         main = "Library Sizes",
            +         ylab = "Total Counts",
            +         xlab = "Sample",
            +         col = c(rep(c("steelblue", "tomato"), 4)),
            +         las = 2)
> legend("topright", legend = c("Before", "After"),
         +        fill = c("steelblue", "tomato"))
> # --- Log2 Count Distribution (Boxplot) ---
  > log_counts <- log2(counts + 1)
>  
  > boxplot(log_counts,
            +         main = "Log2 Count Distribution per Sample",
            +         ylab = "Log2(counts + 1)",
            +         xlab = "Sample",
            +         col = c(rep(c("steelblue", "tomato"), 4)),
            +         las = 2)
> # --- Density Plot ---
  > plot(density(log_counts[, 1]),
         +      main = "Density of Log2 Counts",
         +      xlab = "Log2(counts + 1)",
         +      ylim = c(0, 0.25),
         +      col = "steelblue")
> for (i in 2:ncol(log_counts)) {
  +   lines(density(log_counts[, i]),
            +         col = ifelse(grepl("before", colnames(log_counts)[i]), 
                                   +                      "steelblue", "tomato"))
  + }
>  
  > legend("topright", legend = c("Before", "After"),
           +        col = c("steelblue", "tomato"), lty = 1)
>  
  > # --- Check library sizes and zero counts ---
  > colSums(counts)         # total reads per sample
before_1   after_1  before_2   after_2  before_3   after_3  before_4   after_4 
89425943  73965378  80788168  98145901  79868756  77755385 108984658  75923157 
> colSums(counts == 0)    # number of zero count genes per sample
before_1  after_1 before_2  after_2 before_3  after_3 before_4  after_4 
7180     7155     7186     6896     7329     7429     7303     8361 
> # =============================================================================
> # SECTION 6: SPECIES CHECK & GENE ID VALIDATION
  > # =============================================================================
> # Species: Homo sapiens -> org.Hs.eg.db
  > # Gene IDs: Already HGNC symbols (no conversion needed for this dataset)
  >  
  > # Validate gene symbols against official HGNC database
  > hgnc_symbols <- mapIds(org.Hs.eg.db,
                           +                        keys = rownames(counts),
                           +                        column = "SYMBOL",
                           +                        keytype = "SYMBOL",
                           +                        multiVals = "first")
> # Check for NAs
  > sum(!is.na(hgnc_symbols))  # number mapped
[1] 26364
> sum(is.na(hgnc_symbols))   # number failed
[1] 0
>  
  > # Remove NAs if any (not needed for this dataset - 0 NAs)
  > # counts <- counts[!is.na(hgnc_symbols), ]
  > # hgnc_symbols <- hgnc_symbols[!is.na(hgnc_symbols)]
  >  
  > # Confirm symbols match rownames
  > all(hgnc_symbols == rownames(counts))  # should be TRUE
[1] TRUE
> # =============================================================================
> # SECTION 7: CREATE DGEList OBJECT
  > # =============================================================================
>  
  > dge <- DGEList(counts = counts,
                   +                samples = metadata,
                   +                group = metadata$condition)
>  
  > # Verify
  > dim(dge)        # should be 26364 x 8
[1] 26364     8
> dge$samples     # check sample info and library sizes
group  lib.size norm.factors   sample condition subject
before_1 before  89425943            1 before_1    before       1
after_1   after  73965378            1  after_1     after       1
before_2 before  80788168            1 before_2    before       2
after_2   after  98145901            1  after_2     after       2
before_3 before  79868756            1 before_3    before       3
after_3   after  77755385            1  after_3     after       3
before_4 before 108984658            1 before_4    before       4
after_4   after  75923157            1  after_4     after       4
> # =============================================================================
> # SECTION 8: FILTER LOWLY EXPRESSED GENES
  > # =============================================================================
>  
  > # Filter using edgeR's recommended method
  > keep <- filterByExpr(dge, group = metadata$condition)
> table(keep)  # see how many genes kept vs removed
keep
FALSE  TRUE 
10030 16334 
>  
  > # Apply filter
  > dge <- dge[keep, , keep.lib.sizes = FALSE]
>  
  > # Verify dimensions after filtering
  > dim(dge)  # should be 16338 x 8
[1] 16334     8
> # =============================================================================
> # SECTION 9: TMM NORMALIZATION
  > # =============================================================================
> # TMM = Trimmed Mean of M-values
  > # Corrects for compositional bias between samples
  >  
  > dge <- calcNormFactors(dge, method = "TMM")
>  
  > # Check normalization factors (should be close to 1)
  > dge$samples$norm.factors
[1] 1.1413050 1.3622102 1.1402850 1.1673923 1.0335507 0.9280813 0.8384839 0.6007744
> # =============================================================================
> # SECTION 10: CREATE DESIGN MATRIX
  > # =============================================================================
> # Paired design: blocking by subject to account for individual differences
  > # Formula: ~ subject + condition
  > # Reference: subject 1, condition = after
  >  
  > design2 <- model.matrix(~ 0 + condition + subject, data = metadata)
>  
  > # Clean up column names
  > colnames(design2) <- gsub("condition", "", colnames(design2))
> colnames(design2)
[1] "after"    "before"   "subject2" "subject3" "subject4"
>  
  > # Verify design matrix
  > dim(design2)   # should be 8 x 5
[1] 8 5
> design2
after before subject2 subject3 subject4
before_1     0      1        0        0        0
after_1      1      0        0        0        0
before_2     0      1        1        0        0
after_2      1      0        1        0        0
before_3     0      1        0        1        0
after_3      1      0        0        1        0
before_4     0      1        0        0        1
after_4      1      0        0        0        1
attr(,"assign")
[1] 1 1 2 2 2
attr(,"contrasts")
attr(,"contrasts")$condition
[1] "contr.treatment"

attr(,"contrasts")$subject
[1] "contr.treatment"

> # =============================================================================
> # SECTION 11: VOOM TRANSFORMATION
  > # =============================================================================
> # Converts counts to log2 CPM and estimates mean-variance relationship
  > # Assigns precision weights to each observation
  >  
  > v <- voom(dge, design2, plot = TRUE)
>  
  > # Verify voom output
  > dim(v$E)   # should be 16338 x 8
[1] 16334     8
> head(v$E)  # log2 CPM values
before_1   after_1   before_2    after_2   before_3   after_3
DDX11L1    1.6734782  1.261496  3.0380087  3.1063412  1.1128817  1.098518
WASH7P     4.6355880  4.532041  5.0496425  5.3017942  4.8848797  4.702420
MIR6859-3 -0.9048512 -1.054470 -0.3255152 -0.3560182 -0.5987142 -1.242208
MIR6859-2 -0.9048512 -1.054470 -0.3255152 -0.3560182 -0.5987142 -1.242208
MIR6859-4 -0.9048512 -1.054470 -0.3255152 -0.3560182 -0.5987142 -1.242208
MIR6859-1 -0.9048512 -1.054470 -0.3255152 -0.3560182 -0.5987142 -1.242208
before_4    after_4
DDX11L1    0.7103535  1.3405170
WASH7P     5.0729571  4.9990376
MIR6859-3 -0.2944797 -0.4451429
MIR6859-2 -0.2944797 -0.4451429
MIR6859-4 -0.2944797 -0.4451429
MIR6859-1 -0.2944797 -0.4451429
> # =============================================================================
> # SECTION 12: SAMPLE-LEVEL QC (POST-NORMALIZATION)
  > # =============================================================================
>  
  > # --- MDS Plot ---
  > plotMDS(v,
            +         col = ifelse(metadata$condition == "before", "steelblue", "tomato"),
            +         pch = as.numeric(metadata$subject),
            +         main = "MDS Plot")
>  
  > legend("topright",
           +        legend = c("Before", "After"),
           +        col = c("steelblue", "tomato"),
           +        pch = 1)
> # --- Hierarchical Clustering ---
  > dist_mat <- dist(t(v$E))
> hclust_res <- hclust(dist_mat)
>  
  > plot(hclust_res,
         +      main = "Hierarchical Clustering of Samples",
         +      xlab = "Samples",
         +      sub = "",
         +      cex = 0.9)
> # --- Sample-to-Sample Correlation Heatmap ---
  > cor_mat <- cor(v$E, method = "pearson")

heatmap.2(cor_mat,
          +           trace = "none",
          +           col = colorRampPalette(brewer.pal(9, "Blues"))(100),
          +           margins = c(8, 8),
          +           main = "Sample-to-Sample Correlation",
          +           key.title = "Correlation",
          +           dendrogram = "both",
          +           ColSideColors = ifelse(metadata$condition == "before",
                                             +                                  "steelblue", "tomato"),
          +           RowSideColors = ifelse(metadata$condition == "before",
                                             +                                  "steelblue", "tomato"))
Error in plot.new() : figure margins too large

> # --- Cook's Distance ---
  > # NOTE: Cook's distance was uninformative for this dataset due to model
  > # saturation with n=4 paired samples. Sample QC was assessed via MDS,
  > # hierarchical clustering, and correlation heatmap instead.
  > # after_4 was consistently flagged as different across all QC metrics
  > # but retained due to small sample size (n=4).
  > # =============================================================================
> # SECTION 13: FIT LINEAR MODEL & DEFINE CONTRASTS
  > # =============================================================================
>  
  > # Fit linear model
  > fit2 <- lmFit(v, design2)
> # Define contrasts (after vs before night shift)
  > contrast_matrix <- makeContrasts(
    +   AfterVsBefore = after - before,
    +   levels = design2
    + )
>  
  > contrast_matrix  # verify contrast
Contrasts
Levels     AfterVsBefore
after                1
before              -1
subject2             0
subject3             0
subject4             0
> # Fit contrasts
  > fit2 <- contrasts.fit(fit2, contrast_matrix)
>  
  > # Apply empirical Bayes smoothing
  > # Shrinks gene-wise variances toward a common value - improves power
  > fit2 <- eBayes(fit2)
>  
  > # Check fit object
  > names(fit2)
[1] "coefficients"     "stdev.unscaled"   "sigma"            "df.residual"     
[5] "cov.coefficients" "pivot"            "rank"             "Amean"           
[9] "method"           "design"           "contrasts"        "df.prior"        
[13] "s2.prior"         "var.prior"        "proportion"       "s2.post"         
[17] "t"                "df.total"         "p.value"          "lods"            
[21] "F"                "F.p.value"       
> dim(fit2$coefficients)  # should be 16338 x 1
[1] 16334     1
> # =============================================================================
> # SECTION 14: EXTRACT DIFFERENTIAL EXPRESSION RESULTS
  > # =============================================================================
>  
  > results2 <- topTable(fit2,
                         +                      coef = "AfterVsBefore",
                         +                      number = Inf,
                         +                      adjust.method = "BH",
                         +                      sort.by = "P")
>  
  > # Check results
  > head(results2)
logFC    AveExpr         t      P.Value adj.P.Val         B
EIF1AY  0.6658283  3.3863661  5.403400 0.0003609424 0.5774714 -2.893620
MRC1    1.2224130 -0.6724927  5.289821 0.0004216238 0.5774714 -4.108330
CEBPE  -0.7297401  1.9688379 -5.240555 0.0004512804 0.5774714 -3.501153
FADS2  -0.7492539  2.8566965 -5.079934 0.0005645969 0.5774714 -3.077388
RSPH14  1.1590552 -0.3628088  4.963750 0.0006654545 0.5774714 -4.123997
KRCC1   0.4988557  3.9866808  4.955482 0.0006733330 0.5774714 -2.828447
> dim(results2)
[1] 16334     6
> # =============================================================================
> # SECTION 14: EXTRACT DIFFERENTIAL EXPRESSION RESULTS
  > # =============================================================================
>  
  > results2 <- topTable(fit2,
                         +                      coef = "AfterVsBefore",
                         +                      number = Inf,
                         +                      adjust.method = "BH",
                         +                      sort.by = "P")
>  
  > # Check results
  > head(results2)
logFC    AveExpr         t      P.Value adj.P.Val         B
EIF1AY  0.6658283  3.3863661  5.403400 0.0003609424 0.5774714 -2.893620
MRC1    1.2224130 -0.6724927  5.289821 0.0004216238 0.5774714 -4.108330
CEBPE  -0.7297401  1.9688379 -5.240555 0.0004512804 0.5774714 -3.501153
FADS2  -0.7492539  2.8566965 -5.079934 0.0005645969 0.5774714 -3.077388
RSPH14  1.1590552 -0.3628088  4.963750 0.0006654545 0.5774714 -4.123997
KRCC1   0.4988557  3.9866808  4.955482 0.0006733330 0.5774714 -2.828447
> dim(results2)
[1] 16334     6
> # Summary of significant genes
  > # NOTE: No genes pass adj.P.Val < 0.05 due to small sample size (n=4)
  > sum(results2$adj.P.Val < 0.05)   # adjusted p-value threshold
[1] 0
> sum(results2$P.Value < 0.05)     # raw p-value threshold
[1] 1383
>  
  > # Genes passing paper's threshold (P < 0.05 & |logFC| > 0.5)
  > sig_genes2 <- results2[results2$P.Value < 0.05 &
                             +                          abs(results2$logFC) > 0.5, ]
> nrow(sig_genes2)
[1] 501
> sum(sig_genes2$logFC > 0)  # upregulated after shift
[1] 381
> sum(sig_genes2$logFC < 0)  # downregulated after shift
[1] 120

# =============================================================================
> # SECTION 15: HGNC ANNOTATION & JOIN
  > # =============================================================================
>  
  > # Download HGNC complete set file
  > download.file(
    +   url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    +   destfile = "hgnc_complete_set.txt",
    +   mode = "wb"
    + )
trying URL 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
Content type 'text/plain;charset=utf-8' length 17033504 bytes (16.2 MB)
downloaded 16.2 MB

>  
  > # Load HGNC complete set
  > hgnc_complete_set_df <- read_tsv("hgnc_complete_set.txt")
Rows: 44982 Columns: 54                                                             
── Column specification ────────────────────────────────────────────────────────────
Delimiter: "\t"
chr  (43): hgnc_id, symbol, name, locus_group, locus_type, status, location, loc...
dbl   (4): entrez_id, omim_id, homeodb, orphanet
lgl   (3): kznf_gene_catalog, mamit-trnadb, intermediate_filament_db
date  (4): date_approved_reserved, date_symbol_changed, date_name_changed, date_...

ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
Warning message:
  One or more parsing issues, call `problems()` on your data frame for details, e.g.:
  dat <- vroom(...)
problems(dat) 

> # Check what the parsing issues are
  > problems(hgnc_complete_set_df)
# A tibble: 33 × 5
row   col expected           actual               file                         
<int> <int> <chr>              <chr>                <chr>                        
  1   823    32 a double           312095|465000        C:/Users/nicki/xRStudiox/DEG…
2  5515    32 a double           146910|147010|147070 C:/Users/nicki/xRStudiox/DEG…
3  6012    32 a double           308385|430000        C:/Users/nicki/xRStudiox/DEG…
4 11543    32 a double           300015|402500        C:/Users/nicki/xRStudiox/DEG…
5 11544    32 a double           300162|400011        C:/Users/nicki/xRStudiox/DEG…
6 15197    44 1/0/T/F/TRUE/FALSE 7                    C:/Users/nicki/xRStudiox/DEG…
7 15199    44 1/0/T/F/TRUE/FALSE 12                   C:/Users/nicki/xRStudiox/DEG…
8 15200    44 1/0/T/F/TRUE/FALSE 11                   C:/Users/nicki/xRStudiox/DEG…
9 15201    44 1/0/T/F/TRUE/FALSE 13                   C:/Users/nicki/xRStudiox/DEG…
10 15203    44 1/0/T/F/TRUE/FALSE 20                   C:/Users/nicki/xRStudiox/DEG…
# ℹ 23 more rows
# ℹ Use `print(n = ...)` to see more rows
> # Check dimensions
  > dim(hgnc_complete_set_df)  # should be 44982 x 54
[1] 44982    54
> 
  > # Check it still has the symbol column
  > "symbol" %in% colnames(hgnc_complete_set_df)
[1] TRUE
> # Rename symbol column
  > colnames(hgnc_complete_set_df)[colnames(hgnc_complete_set_df) == "symbol"] <- "HGNC_Symbol"
> # Convert results to tibble with gene symbols as column
  > res_df <- results2 %>%
  +   as_tibble(rownames = "HGNC_Symbol")
> # Perform inner join with HGNC database
  > joined_results <- res_df %>%
  +   inner_join(dplyr::select(hgnc_complete_set_df, HGNC_Symbol),
                 +              by = "HGNC_Symbol")
> # =============================================================================
> # SECTION 16: HANDLE DUPLICATES & FINALIZE RESULTS TABLE
  > # =============================================================================
>  
  > # Handle duplicate gene symbols (keep lowest adj.P.Val)
  > # Remove NAs and select final columns
  > final_results <- joined_results %>%
  +   dplyr::arrange(adj.P.Val) %>%
  +   dplyr::distinct(HGNC_Symbol, .keep_all = TRUE) %>%
  +   dplyr::select(HGNC_Symbol, logFC, P.Value, adj.P.Val) %>%
  +   dplyr::rename(raw_P.value = P.Value,
                    +                 adj_P.value = adj.P.Val)
> # Verify final table
  > colnames(final_results)
[1] "HGNC_Symbol" "logFC"       "raw_P.value" "adj_P.value"
> dim(final_results)
[1] 14628     4
> head(final_results)
# A tibble: 6 × 4
HGNC_Symbol  logFC raw_P.value adj_P.value
<chr>        <dbl>       <dbl>       <dbl>
  1 EIF1AY       0.666    0.000361       0.577
2 MRC1         1.22     0.000422       0.577
3 CEBPE       -0.730    0.000451       0.577
4 FADS2       -0.749    0.000565       0.577
5 RSPH14       1.16     0.000665       0.577
6 KRCC1        0.499    0.000673       0.577
> sum(is.na(final_results$adj_P.value))  # should be 0
[1] 0
>  
  > # Save final results table
  > write.csv(final_results,
              +           "GSE282051_DGE_results.csv",
              +           row.names = FALSE)
> # =============================================================================
> # SECTION 17: VISUALIZATION
  > # =============================================================================
>  
  > # --- Volcano Plot ---
  > final_results <- final_results %>%
  +   dplyr::mutate(significance = case_when(
    +     raw_P.value < 0.05 & logFC > 0.5  ~ "Up after shift",
    +     raw_P.value < 0.05 & logFC < -0.5 ~ "Down after shift",
    +     TRUE ~ "Not significant"
    +   ))
> ggplot(final_results, aes(x = logFC, y = -log10(raw_P.value),
                            +                            color = significance)) +
  +   geom_point(alpha = 0.6, size = 1.5) +
  +   scale_color_manual(values = c("Up after shift" = "tomato",
                                    +                                 "Down after shift" = "steelblue",
                                    +                                 "Not significant" = "grey70")) +
  +   geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 +              color = "black") +
  +   geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
                 +              color = "black") +
  +   labs(title = "Volcano Plot - Night Shift Gene Expression",
           +        subtitle = "After vs Before Night Shift",
           +        x = "Log2 Fold Change",
           +        y = "-Log10(P-value)",
           +        color = "Direction") +
  +   theme_bw() +
  +   theme(legend.position = "right")
Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
  invalid graphics state

> # Reset graphics device
  > graphics.off()
> 
  > # Then rerun the volcano plot
  > final_results <- final_results %>%
  +     dplyr::mutate(significance = case_when(
    +         raw_P.value < 0.05 & logFC > 0.5  ~ "Up after shift",
    +         raw_P.value < 0.05 & logFC < -0.5 ~ "Down after shift",
    +         TRUE ~ "Not significant"
    +     ))
> 
  > ggplot(final_results, aes(x = logFC, y = -log10(raw_P.value),
                              +                           color = significance)) +
  +     geom_point(alpha = 0.6, size = 1.5) +
  +     scale_color_manual(values = c("Up after shift" = "tomato",
                                      +                                   "Down after shift" = "steelblue",
                                      +                                   "Not significant" = "grey70")) +
  +     geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   +                color = "black") +
  +     geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
                   +                color = "black") +
  +     labs(title = "Volcano Plot - Night Shift Gene Expression",
             +          subtitle = "After vs Before Night Shift",
             +          x = "Log2 Fold Change",
             +          y = "-Log10(P-value)",
             +          color = "Direction") +
  +     theme_bw() +
  +     theme(legend.position = "right")
> # --- Heatmap of Top 50 DE Genes ---
  >  
  > # Get top 50 significant genes
  > top_genes <- final_results %>%
  +   dplyr::filter(raw_P.value < 0.05 & abs(logFC) > 0.5) %>%
  +   dplyr::arrange(raw_P.value) %>%
  +   head(50)
> # Extract and scale expression data
  > top_expr <- v$E[top_genes$HGNC_Symbol, ]
> top_expr_scaled <- t(scale(t(top_expr)))
>  
  > # Order columns by condition
  > col_order <- c("before_1", "before_2", "before_3", "before_4",
                   +                "after_1", "after_2", "after_3", "after_4")
> top_expr_ordered <- top_expr_scaled[, col_order]
>  
  > # Average before and after samples
  > before_avg <- rowMeans(top_expr_ordered[, grepl("before", colnames(top_expr_ordered))])
> after_avg <- rowMeans(top_expr_ordered[, grepl("after", colnames(top_expr_ordered))])
>  
  > top_expr_avg <- data.frame(Before = before_avg, After = after_avg)
>  
  > # Prepare for ggplot
  > heatmap_data_avg <- top_expr_avg %>%
  +   tibble::rownames_to_column("Gene") %>%
  +   tidyr::pivot_longer(cols = -Gene,
                          +                       names_to = "Condition",
                          +                       values_to = "Expression")
>  
  > # Plot heatmap
  > ggplot(heatmap_data_avg, aes(x = Condition, y = Gene, fill = Expression)) +
  +   geom_tile() +
  +   scale_fill_gradient2(low = "#053061",
                           +                        mid = "white",
                           +                        high = "#67001F",
                           +                        midpoint = 0,
                           +                        limits = c(-2, 2),
                           +                        name = "Z-score") +
  +   labs(title = "Top 50 DE Genes - Night Shift",
           +        x = "Condition", y = "Gene") +
  +   theme_bw() +
  +   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            +         axis.text.y = element_text(size = 9),
            +         strip.text = element_text(face = "bold"),
            +         plot.margin = margin(10, 10, 10, 60),
            +         panel.border = element_blank(),
            +         axis.line = element_blank(),
            +         panel.grid = element_blank()) +
  +   scale_y_discrete(expand = expansion(add = 0.5))
>  
  > ggsave("GSE282051_heatmap.png", width = 6, height = 14, dpi = 300)
> # =============================================================================
> # SECTION 18: SAVE WORKSPACE
  > # =============================================================================
>  
  > save.image("GSE282051_analysis.RData")
>  
  > # To reload in a future session:
  > # load("GSE282051_analysis.RData")
  > # Then reload all libraries listed in Section 1
  >  
  > # =============================================================================
> # END OF SOP
  > # =============================================================================
> 