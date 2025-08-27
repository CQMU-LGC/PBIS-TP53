###Figure 6A(Cycling cell is shown here as an example; the same applies to other cell types.)
library(SCPA)
library(msigdbr)
library(Seurat)
library(dplyr)
library(ggplot2)
Idents(scobj)<-"Cell_Type"
seurat_object <- subset(scobj, idents=c('Cycling cell'))
low<-seurat_extract(seurat_object ,
                      meta1="group",value_meta1="Low")
high<-seurat_extract(seurat_object ,
                    meta1="group",value_meta1="High")
pathways<-msigdbr("Homo sapiens","H")%>%
    format_pathways()
scpa_out<-compare_pathways(samples=list(low,high),
                             pathways=pathways)
plot_rank(scpa_out=scpa_out,
          pathway=c("Allograft Rejection"),
          base_point_size=3,
          highlight_point_size=5)

###Figure 6B
library(Seurat)
library(SeuratData)
library(RcppML)
library(irGSEA)
scobj<-irGSEA.score(object=scobj,assay="RNA",
                             slot="data",seeds=123,ncores=1,
                             min.cells=3,min.feature=0,
                             custom=F,geneset=NULL,msigdb=T,
                             species="Homo sapiens",category="H",
                             subcategory=NULL,geneid="symbol",
                             method=c("ssgsea"),
                             aucell.MaxRank=NULL,ucell.MaxRank=NULL,
                             kcdf='Gaussian')
result.dge<-irGSEA.integrate(object=scobj,
                               group.by="cell_groups",
                               method=c("ssgsea"))
irGSEA.heatmap.plot<-irGSEA.heatmap(object=result.dge,
                                      method="ssgsea",
                                      top=50,
                                      show.geneset=NULL)
###Figure 6C
irGSEA.upset.plot<-irGSEA.upset(object=result.dge,
                                  method="ssgsea",
                                  mode="intersect",
                                  upset.width=20,
                                  upset.height=10,
                                  set.degree=2,
                                  pt_size=grid::unit(2,"mm"))

###Figure 6D
library(stringr)
library(Seurat)
library(patchwork)
library(SummarizedExperiment)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(SCENIC)
library(qs)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads=8)
scobj <- readRDS("scobj.rds")
saveRDS(Epithelial,file="Epithelial.rds")
scobj <- readRDS("Epithelial.rds")
scRNA <- SetupForWGCNA(
  scobj,
  gene_select ="fraction",
  fraction =0.05,
  wgcna_name = "KICH_KIRC"
)
scRNA<- MetacellsByGroups(
  seurat_obj = scRNA,k=25,
  max_shared = 15,
  reduction ='harmony',
  group.by = c("hdWGCNA",'orig.ident'), 
  ident.group = 'hdWGCNA' 
)
scRNA <- NormalizeMetacells(scRNA)
seurat_obj  <- SetDatExpr(
  scRNA,assay = 'RNA',slot = 'data',
  group_name = "KICH_KIRC", 
  group.by='hdWGCNA'
)
seurat_obj<-TestSoftPowers(
  seurat_obj,
  networkType="signed"
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)

###Figure 6E
softpower=6
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power=softpower,
  setDatExpr=FALSE,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  detectCutHeight = 0.995,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  tom_outdir = "TOM", #
  tom_name = 'KICH_KIRC' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main='KICH_KIRC hdWGCNA Dendrogram')

###Figure 6F
library(fmsb)
library(scales)
analyze_correlation_and_plot <- function(data_file,
                                         target_col = "PBIS.TP53",
                                         output_table = "correlation_results.txt",
                                         output_plot = "correlation_radar.pdf",
                                         color = "#1B9E77") {
  # --- Step 1: Load data
  sample_df <- read.table(data_file, header = TRUE, sep = "\t")
  # --- Step 2: variable names (exclude target)
  var_names <- setdiff(colnames(sample_df), target_col)
  # --- Step 3: Safe Shapiro test
  shapiro_safe <- function(x, max_n = 5000) {
    x <- na.omit(x)
    n <- length(x)
    if (n < 3) {
      return(list(p.value = 0))  # 不满足正态
    } else if (n > max_n) {
      x <- sample(x, max_n)
    }
    return(shapiro.test(x))
  }
  # --- Step 4: Compute correlations
  cor_vals <- numeric(length(var_names))
  p_vals <- numeric(length(var_names))
  methods_used <- character(length(var_names))

  for (i in seq_along(var_names)) {
    x <- as.numeric(sample_df[[target_col]])
    y <- as.numeric(sample_df[[var_names[i]]])
    shapiro_x <- shapiro_safe(x)
    shapiro_y <- shapiro_safe(y)
    method <- if (shapiro_x$p.value > 0.05 & shapiro_y$p.value > 0.05) "pearson" else "spearman"
    methods_used[i] <- method
    cor_test <- cor.test(x, y, method = method)
    cor_vals[i] <- cor_test$estimate
    p_vals[i] <- cor_test$p.value
  }
  # --- Step 5: Significance stars
  get_stars <- function(p) {
    if (p <= 0.001) return("***")
    else if (p <= 0.01) return("**")
    else if (p <= 0.05) return("*")
    else return("")
  }
  stars <- sapply(p_vals, get_stars)
  # --- Step 6: Save results table
  cor_results <- data.frame(
    Module = var_names,
    Correlation = cor_vals,
    P_value = p_vals,
    Method = methods_used,
    Significance = stars
  )
  write.table(cor_results, file = output_table, sep = "\t", row.names = FALSE, quote = FALSE)
  # --- Step 7: Prepare data for radar chart
  scaled_vals <- (cor_vals + 1) / 2
  radar_data <- as.data.frame(rbind(rep(1, length(var_names)),
                                    rep(0, length(var_names)),
                                    scaled_vals))
  colnames(radar_data) <- var_names
  rownames(radar_data) <- c("Max", "Min", "Cor")
  # --- Step 8: Plot radar chart
  plot_correlation_radarchart <- function(data, stars, color) {
    radarchart(data,
               axistype = 1,
               pcol = color,
               pfcol = alpha(color, 0.5),
               plwd = 2,
               cglcol = "grey",
               cglty = 1,
               cglwd = 0.8,
               axislabcol = "black",
               caxislabels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
               vlcex = 0.8,
               title = "Correlation Radar Chart")
    # Add significance stars
    n_vars <- ncol(data)
    for (i in 1:n_vars) {
      angle <- (2 * pi / n_vars) * (i - 1)
      r <- 1.05
      x <- r * sin(angle)
      y <- r * cos(angle)
      text(x, y, labels = stars[i], cex = 1.4, col = "red")
    }
  }
  pdf(file = output_plot, width = 5, height = 5)
  par(mar = c(1, 1, 2, 1))
  plot_correlation_radarchart(radar_data, stars, color)
  dev.off()
  message("✅ Correlation analysis finished.")
  message("Results table saved to: ", output_table)
  message("Radar chart saved to: ", output_plot)
  return(cor_results)
}
results <- analyze_correlation_and_plot(
  data_file = "Module.txt",
  target_col = "PBIS.TP53",
  output_table = "PBIS.TP53_vs_modules.txt",
  output_plot = "Figure 6F.pdf",
  color = "#1B9E77"
)

###Figure 6G
library(tidyverse)
library(corrplot)
library(RColorBrewer)
analyze_celltype_correlation <- function(data_file,
                                         target_col = "PBIS.TP53",
                                         celltype_col = "Cell_Type",
                                         exclude_cols = c("Cell_Type", "PBIS.TP53"),
                                         output_pdf = "CellType_Module_Correlation.pdf",
                                         max_n = 5000,
                                         width = 5,
                                         height = 5) {
  # --- Step 1: Load data
  sample_df <- read.table(data_file, header = TRUE, sep = "\t")
  # --- Step 2: Extract module names
  module_names <- colnames(sample_df)[!(colnames(sample_df) %in% exclude_cols)]
  # --- Step 3: Extract cell types
  cell_types <- unique(sample_df[[celltype_col]])
  # --- Step 4: Initialize matrices
  cor_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(module_names),
                       dimnames = list(cell_types, module_names))
  pval_matrix <- cor_matrix
  # --- Step 5: Loop through each cell type
  for (cell in cell_types) {
    sub_df <- sample_df[sample_df[[celltype_col]] == cell, ]
    x <- sub_df[[target_col]]
    for (mod in module_names) {
      y <- sub_df[[mod]]
      # downsample if too large
      n_used <- min(length(x), max_n)
      sample_x <- sample(x, n_used)
      sample_y <- sample(y, n_used)
      # test for normality
      if (shapiro.test(sample_x)$p.value > 0.05 &
          shapiro.test(sample_y)$p.value > 0.05) {
        method <- "pearson"
      } else {
        method <- "spearman"
      }
      cor_result <- suppressWarnings(cor.test(x, y, method = method))
      cor_matrix[cell, mod] <- cor_result$estimate
      pval_matrix[cell, mod] <- cor_result$p.value
    }
  }
  # --- Step 6: Convert to dataframe
  cor_df <- as.data.frame(cor_matrix)
  pval_df <- as.data.frame(pval_matrix)
  # --- Step 7: Plot correlation heatmap
  col <- colorRampPalette(brewer.pal(11,"RdBu")[3:9])(100)
  pdf(file = output_pdf, width = width, height = height)
  corrplot(as.matrix(cor_df),
           p.mat = as.matrix(pval_df),
           method = "square",
           insig = "label_sig",
           sig.level = c(0.001, 0.01, 0.05),
           pch.cex = 0.9,
           pch.col = "white",
           is.corr = FALSE,
           cl.length = 6,
           col = col,
           tl.col = "black",
           tl.cex = 0.7)
  dev.off()
  message("✅ Correlation analysis finished.")
  message("Correlation matrix and p-values computed for ", length(cell_types), " cell types.")
  message("PDF saved to: ", output_pdf)
  return(list(cor_matrix = cor_df, pval_matrix = pval_df))
}
results <- analyze_celltype_correlation(
  data_file = "module.txt",
  target_col = "PBIS.TP53",
  celltype_col = "Cell_Type",
  exclude_cols = c("Cell_Type", "PBIS.TP53"),
  output_pdf = "Figure 6G.pdf"
)

###Figure6 (H-j)(Through String database)






