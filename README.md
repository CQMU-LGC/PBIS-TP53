# 🧬 PBIS-TP53: Pathway-Based Integrated Scoring for TP53
This repository provides an R-based pipeline to compute the **PBIS-TP53** score, a transcriptional signature that quantifies the downstream impact of TP53 mutations using weighted ssGSEA.
---
## 📦 Requirements
This script requires the following R packages:
- `GSVA`
- `Matrix`
```r
install.packages("BiocManager")
BiocManager::install("GSVA")
install.packages("Matrix")
calculate_total_score <- function(expr_file, geneSet, coef_file, pos_file, neg_file) {
  library(GSVA)
  library(Matrix)
  exprSet <- read.table(file = expr_file, sep = "\t", header = TRUE, row.names = 1)
  X <- t(exprSet)
  if (length(colnames(X)) != ncol(X)) {
    if (!is.null(colnames(X))) {
      colnames(X) <- colnames(X)[1:ncol(X)]
    } else {
      colnames(X) <- paste0("V", 1:ncol(X))
    }
  }
  incidenceMatrix <- function(X, group) {
    n <- nrow(X)
    p <- ncol(X)
    if (!is.list(group)) stop("Argument 'group' must be a list!")
    J <- length(group)
    grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE)
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:ncol(X), sep = "")
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
    if (all(grp.mat == 0)) stop("Group names don't match expression matrix!")
    return(grp.mat)
  }
  expandX <- function(X, group) {
    incidence.mat <- incidenceMatrix(X, group)
    incidence.mat <- as(incidence.mat, "CsparseMatrix") 
    over.mat <- incidence.mat %*% t(incidence.mat)
    grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
    X.latent <- NULL
    names <- NULL
    for (i in 1:nrow(incidence.mat)) {
      idx <- incidence.mat[i, ] == 1
      X.latent <- cbind(X.latent, X[, idx, drop = FALSE])
      names <- c(names, colnames(incidence.mat)[idx])
    }
    colnames(X.latent) <- paste("grp", grp.vec, "_", names, sep = "")
    return(as.matrix(X.latent))
  }
  X.latent <- expandX(X, geneSet)
  weighted_expression_data <- coef_file
  weighted_expression_data$Weights <- abs(weighted_expression_data$coefficients) + 1
  expression_data <- X.latent
  gene_weights <- setNames(weighted_expression_data$Weights, rownames(weighted_expression_data))
  if (!setequal(colnames(expression_data), names(gene_weights))) {
    stop("Mismatch between gene weights and expression matrix columns.")
  }
  gene_weights <- gene_weights[colnames(expression_data)]
  weighted_expression_matrix <- expression_data
  for (gene in colnames(expression_data)) {
    weighted_expression_matrix[, gene] <- expression_data[, gene] * gene_weights[gene]
  }
  gene_sets <- list(Pos = pos_file$x, Neg = neg_file$x)
  cc <- t(weighted_expression_matrix)
  Enrichment_score <- gsva(expr = as.matrix(cc), gene_sets, kcdf = "Gaussian", method = "ssgsea", abs.ranking = TRUE)
  score <- Enrichment_score["Pos", ] - Enrichment_score["Neg", ]
  score <- as.data.frame(score)
  Enrichment_score <- t(Enrichment_score)
  if (!all(rownames(score) == rownames(Enrichment_score))) {
    stop("Score sample names do not match enrichment results!")
  }
  Enrichment_score <- as.data.frame(Enrichment_score)
  Enrichment_score$Score <- score$score
  Enrichment_score$Sample <- rownames(Enrichment_score)
  return(Enrichment_score)
}
filtered_gene_sets <- readRDS("filtered_gene_sets.rds")
coef_file <- readRDS("coef_file.rds")
pos_file <- readRDS("pos_file.rds")
neg_file <- readRDS("neg_file.rds")
expr_file <- "TPM.txt"
result <- calculate_total_score(expr_file, filtered_gene_sets, coef_file, pos_file, neg_file)
write.table(result, file = "PBIS-TP53.txt", sep = "\t", quote = FALSE)
