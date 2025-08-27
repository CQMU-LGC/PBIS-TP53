# 🧬 PBIS-TP53: Pathway-Based Integrated Scoring for TP53
This repository provides an R-based pipeline to compute the **PBIS-TP53** score, a transcriptional signature that quantifies the downstream impact of TP53 mutations using weighted ssGSEA.
---
![PBIS-TP53 Overview](https://raw.githubusercontent.com/CQMU-LGC/PBIS-TP53/main/PBIS-TP53.png?raw=true&v=2)
## 📦 Requirements
This script requires the following R packages:
- `GSVA`
- `Matrix`
```r
install.packages("BiocManager")
BiocManager::install("GSVA")
install.packages("Matrix")
# Calculation of PBIS-TP53 in Bulk data
library(Matrix)
library(GSVA)
# Function definition
run_weighted_gsva <- function(expr_file,
                              group_rds,
                              coef_file,
                              pos_file,
                              neg_file,
                              output_file = "Bulk_TP53.txt") {
  # --- Step 1: Load expression data ---
  exprSet <- read.table(file = expr_file, sep = "\t", header = TRUE, row.names = 1)
  X <- t(exprSet)
  # Ensure proper column names
  if (length(colnames(X)) != ncol(X)) {
    if (!is.null(colnames(X))) {
      colnames(X) <- colnames(X)[1:ncol(X)]
    } else {
      colnames(X) <- paste0("V", 1:ncol(X))
    }
  }
  # --- Step 2: incidenceMatrix helper ---
  incidenceMatrix <- function(X, group) {
    n <- nrow(X); p <- ncol(X)
    if (!is.list(group)) stop("Argument 'group' must be a list!")
    J <- length(group)
    grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE)
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:p, sep = "")
    if (is.null(names(group))) names(group) <- paste("grp", 1:J, sep = "")
    if (is.numeric(group[[1]])) {
      for (i in 1:J) {
        ind <- group[[i]]
        grp.mat[i, ind] <- 1
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    } else {
      for (i in 1:J) {
        grp.i <- as.character(group[[i]])
        ind <- colnames(X) %in% grp.i
        grp.mat[i, ] <- 1 * ind
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    }
    rownames(grp.mat) <- names(group)
    if (all(grp.mat == 0)) stop("The names in X do not match the group list!")
    grp.mat
  }
  # --- Step 3: expandX helper ---
  expandX <- function(X, group) {
    incidence.mat <- incidenceMatrix(X, group)
    over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE)
    grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
    X.latent <- NULL
    names <- NULL
    for (i in 1:nrow(incidence.mat)) {
      idx <- incidence.mat[i, ] == 1
      X.latent <- cbind(X.latent, X[, idx, drop = FALSE])
      names <- c(names, colnames(incidence.mat)[idx])
    }
    colnames(X.latent) <- paste("grp", grp.vec, "_", names, sep = "")
    X.latent
  }
  # --- Step 4: Build latent expression ---
  group <- readRDS(group_rds)
  expression_data <- expandX(X, group)
  # --- Step 5: Read SGL outputs ---
  if (!file.exists(coef_file) || !file.exists(pos_file) || !file.exists(neg_file)) {
    stop("One or more SGL input files are missing!")
  }
  weighted_expression_data <- read.table(coef_file, header = TRUE, sep = "\t")
  weighted_expression_data$Weights <- abs(weighted_expression_data$coefficients) + 1
  Pos <- read.table(pos_file, header = TRUE, sep = "\t")
  Neg <- read.table(neg_file, header = TRUE, sep = "\t")
  gene_weights <- setNames(weighted_expression_data$Weights, rownames(weighted_expression_data))
  if (!setequal(colnames(expression_data), names(gene_weights))) {
    stop("Mismatch between gene weights and expression data columns!")
  }
  gene_weights <- gene_weights[colnames(expression_data)]
  # --- Step 6: Apply weights ---
  weighted_expression_matrix <- expression_data
  for (gene in colnames(expression_data)) {
    weighted_expression_matrix[, gene] <- expression_data[, gene] * gene_weights[gene]
  }
  # --- Step 7: GSVA enrichment ---
  cc <- t(weighted_expression_matrix)
  Enrichment_score <- gsva(expr = as.matrix(cc),
                           list(Pos = Pos$x, Neg = Neg$x),
                           kcdf = "Gaussian", method = "ssgsea")

  score <- Enrichment_score["Pos", ] - Enrichment_score["Neg", ]
  Enrichment_score <- as.data.frame(t(Enrichment_score))
  Enrichment_score$Score <- score
  Enrichment_score$Sample <- rownames(Enrichment_score)
  # --- Step 8: Save output ---
  write.table(Enrichment_score, file = output_file, sep = "\t", quote = FALSE)
  return(Enrichment_score)
}
res <- run_weighted_gsva(
  expr_file = "Expression matrix.txt",
  group_rds = "filtered_gene_sets.rds",
  coef_file = "coefficients_df(0.65_0.02).txt",
  pos_file = "SGL_pos(0.65_0.02).txt",
  neg_file = "SGL_neg(0.65_0.02).txt",
  output_file = "Bulk_TP53.txt"
)

# Calculation of PBIS-TP53 in scRNA-seq data
library(Seurat)
library(Matrix)
library(GSVA)
# Function definition
run_sc_weighted_gsva <- function(seurat_file,
                                 gene_length_file,
                                 group_rds,
                                 coef_file,
                                 pos_file,
                                 neg_file,
                                 output_file = "PBIS-TP53.txt") {
  # --- Step 1: Load Seurat object and extract counts ---
  scobj <- readRDS(seurat_file)
  counts <- GetAssayData(object = scobj, slot = "counts")
  # --- Step 2: TPM normalization ---
  gene_lengths <- read.table(gene_length_file, header = TRUE, sep = "\t")
  genes_in_counts <- rownames(counts)
  matched_gene_lengths <- gene_lengths[gene_lengths$genesymbol %in% genes_in_counts, ]
  counts <- counts[matched_gene_lengths$genesymbol, ]
  gene_lengths_vector <- matched_gene_lengths$length
  RPK <- counts / (gene_lengths_vector / 1000)
  TPM <- sweep(RPK, 2, colSums(RPK), FUN = "/") * 1e6
  TPM <- log2(TPM + 1)
  rownames(TPM) <- gsub("\\.", "-", rownames(TPM))
  # --- Step 3: transpose expression for downstream use ---
  X <- t(as.matrix(TPM))
  if (length(colnames(X)) != ncol(X)) {
    if (!is.null(colnames(X))) {
      colnames(X) <- colnames(X)[1:ncol(X)]
    } else {
      colnames(X) <- paste0("V", 1:ncol(X))
    }
  }
  # --- Step 4: helper functions for grouping ---
  incidenceMatrix <- function(X, group) {
    n <- nrow(X); p <- ncol(X)
    if (!is.list(group)) stop("Argument 'group' must be a list!")
    J <- length(group)
    grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE)
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:p, sep = "")
    if (is.null(names(group))) names(group) <- paste("grp", 1:J, sep = "")
    if (is.numeric(group[[1]])) {
      for (i in 1:J) {
        ind <- group[[i]]
        grp.mat[i, ind] <- 1
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    } else {
      for (i in 1:J) {
        grp.i <- as.character(group[[i]])
        ind <- colnames(X) %in% grp.i
        grp.mat[i, ] <- 1 * ind
        colnames(grp.mat)[ind] <- colnames(X)[ind]
      }
    }
    rownames(grp.mat) <- as.character(names(group))
    if (all(grp.mat == 0)) stop("The names in X do not match the group list!")
    grp.mat
  }
  expandX <- function(X, group) {
    incidence.mat <- incidenceMatrix(X, group)
    over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE)
    grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
    X.latent <- NULL
    names <- NULL
    for (i in 1:nrow(incidence.mat)) {
      idx <- incidence.mat[i, ] == 1
      X.latent <- cbind(X.latent, X[, idx, drop = FALSE])
      names <- c(names, colnames(incidence.mat)[idx])
    }
    colnames(X.latent) <- paste("grp", grp.vec, "_", names, sep = "")
    X.latent
  }
  # --- Step 5: expand expression by groups ---
  group <- readRDS(group_rds)
  expression_data <- expandX(X, group)
  # --- Step 6: read coefficient/pos/neg files ---
  if (!file.exists(coef_file) || !file.exists(pos_file) || !file.exists(neg_file)) {
    stop("One or more input files are missing!")
  }
  weighted_expression_data <- read.table(coef_file, header = TRUE, sep = "\t")
  weighted_expression_data$Weights <- abs(weighted_expression_data$coefficients) + 1
  Pos <- read.table(pos_file, header = TRUE, sep = "\t")
  Neg <- read.table(neg_file, header = TRUE, sep = "\t")
  gene_weights <- setNames(weighted_expression_data$Weights, rownames(weighted_expression_data))
  if (!setequal(colnames(expression_data), names(gene_weights))) {
    stop("Mismatch between gene weights and expression data columns!")
  }
  gene_weights <- gene_weights[colnames(expression_data)]
  # --- Step 7: apply weights ---
  weighted_expression_matrix <- expression_data
  for (gene in colnames(expression_data)) {
    weighted_expression_matrix[, gene] <- expression_data[, gene] * gene_weights[gene]
  } 
  # --- Step 8: run GSVA (ssGSEA) ---
  cc <- t(weighted_expression_matrix)
  Enrichment_score <- gsva(expr = as.matrix(cc),
                           list(Pos = Pos$x, Neg = Neg$x),
                           kcdf = "Gaussian", method = "ssgsea")  
  score <- Enrichment_score["Pos", ] - Enrichment_score["Neg", ]
  Enrichment_score <- as.data.frame(t(Enrichment_score))
  Enrichment_score$Score <- score
  Enrichment_score$Sample <- rownames(Enrichment_score) 
  # --- Step 9: save results ---
  write.table(Enrichment_score, file = output_file, sep = "\t", quote = FALSE)
  return(Enrichment_score)
}
res <- run_sc_weighted_gsva(
  seurat_file = "scRNA-seq.rds",
  gene_length_file = "gencode.v22.annotation.txt",
  group_rds = "filtered_gene_sets.rds",
  coef_file = "coefficients_df(0.65_0.02).txt",
  pos_file = "SGL_pos(0.65_0.02).txt",
  neg_file = "SGL_neg(0.65_0.02).txt",
  output_file = "PBIS-TP53.txt"
)
