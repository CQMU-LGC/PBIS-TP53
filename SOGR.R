###
library(GSA)
exprSet <- read.table("intersection_genes.txt", sep = "\t", header = TRUE)
gmt_file <- "ReactomePathways.gmt"  
gmt_data <- GSA.read.gmt(gmt_file)
pathway_list <- setNames(gmt_data$genesets, gmt_data$geneset.names)
pathway_list <- lapply(pathway_list, function(x) x[nchar(x) > 0])  
l <- Filter(function(x) length(x) > 0, pathway_list)  
gene_sets <- lapply(l, unique)  
filtered_gene_sets <- list()
for (gene_set_name in names(gene_sets)) {
    genes <- gene_sets[[gene_set_name]]
    filtered_genes <- genes[genes %in% exprSet$genesymbol]  
    if (length(filtered_genes) >= 10 && length(filtered_genes) <= 500) {
        filtered_gene_sets[[gene_set_name]] <- filtered_genes
    }
}
if (length(filtered_gene_sets) > 0) {
    max_length <- max(lengths(filtered_gene_sets))
    gene_df <- do.call(cbind, lapply(filtered_gene_sets, function(x) {
        c(x, rep(NA, max_length - length(x)))
    }))
    colnames(gene_df) <- names(filtered_gene_sets)
    print(gene_df)
    write.csv(gene_df, "filtered_gene_sets.csv", row.names = FALSE, na = "")
} else {
    print("No gene sets passed the filtering criteria.")
}

###
library(clusterProfiler)
library(org.Hs.eg.db)  
library(msigdbr)      
library(openxlsx)      
library(readxl)     
gmt_file <- "ReactomePathways.gmt"  
term2gene <- read.gmt(gmt_file)
file_path <- "filtered_gene_sets.xlsx"
gene_sets_data <- read_excel(file_path)
results_list <- list()
for (gene_set_name in colnames(gene_sets_data)) {
    gene_list <- na.omit(gene_sets_data[[gene_set_name]]) 
    if (length(gene_list) > 0) {  
        enrich_res <- enricher(
            gene = gene_list, 
            TERM2GENE = term2gene, 
            pvalueCutoff = 0.05,
            minGSSize = 10,    
            maxGSSize = 500 
        )
        if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
            significant_pathways <- enrich_res@result[enrich_res@result$p.adjust < 0.05, ]
            if (nrow(significant_pathways) > 0) {
                top_pathway <- significant_pathways[order(significant_pathways$p.adjust), ][1, ]   
                results_list[[gene_set_name]] <- data.frame(
                    Name = gene_set_name,
                    Pathway = top_pathway$Description,
                    adj.p = top_pathway$p.adjust,
                    Gene_Num = length(gene_list)
                )
            }
        }
    }
}
final_results <- do.call(rbind, results_list)
write.csv(final_results, "gene_set_enrichment_results.csv", row.names = FALSE)

###
library(dplyr)
library(readxl)
df <- read_excel("gene_set_enrichment_results.xlsx")
df_filtered <- df %>% filter(Name == Pathway)
write.csv(df_filtered, "df_filtered.csv", row.names = FALSE)

###
library(readxl)  
library(GSA)
exprSet <- read.table("intersection_genes.txt", sep = "\t", header = TRUE)
gmt_file <- "ReactomePathways.gmt"  
gmt_data <- GSA.read.gmt(gmt_file)
pathway_list <- setNames(gmt_data$genesets, gmt_data$geneset.names)
pathway_list <- lapply(pathway_list, function(x) x[nchar(x) > 0])  
l <- Filter(function(x) length(x) > 0, pathway_list)  
gene_sets <- lapply(l, unique)  
filtered_gene_sets <- list()
for (gene_set_name in names(gene_sets)) {
    genes <- gene_sets[[gene_set_name]]
    filtered_genes <- genes[genes %in% exprSet$genesymbol]    
    if (length(filtered_genes) >= 10 && length(filtered_genes) <= 500) {
        filtered_gene_sets[[gene_set_name]] <- filtered_genes
    }
}
df <- read_excel("df_filtered.xlsx")
selected_names <- df$Pathway 
selected_gene_sets <- filtered_gene_sets[names(filtered_gene_sets) %in% selected_names]
saveRDS(selected_gene_sets,file="filtered_gene_sets.rds")

###SOGR
library(Matrix)
library(asgl)
exprSet <- read.table(file = "TCGA-TPM(normalize).txt", sep = "\t", header = T, row.names = 1)
X=t(exprSet)
y <- read.table("phenotype.txt", header = TRUE, sep = "\t")
identical(y$sample, colnames(exprSet))
y <- y$group
if (length(colnames(X)) != ncol(X)) {
  if (!is.null(colnames(X))) {
    colnames(X) <- colnames(X)[1:ncol(X)]
  } else {
    colnames(X) <- paste0("V", 1:ncol(X))
  }
}
incidenceMatrix <- function (X, group) {
  n <- nrow(X)
  p <- ncol(X)
  if (!is.list(group)) {
    stop("Argument 'group' must be a list of integer indices or character names of variables!")
  }
  J <- length(group)
  grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE, dimnames = list(as.character(rep(NA, J)), as.character(rep(NA, p))))
  if (is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep = "")
  }
  if (is.null(names(group))) {
    names(group) <- paste("grp", 1:J, sep = "")
  }
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
  if (all(grp.mat == 0)) {
    stop("The names of variables in X don't match with names in group!")
  }
  print("incidenceMatrix dimensions:")
  print(dim(grp.mat))  
  print("incidenceMatrix:")
  print(grp.mat)
  grp.mat
}
expandX <- function (X, group) {
  incidence.mat <- incidenceMatrix(X, group)
  over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE)
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
  X.latent <- NULL
  names <- NULL
  for (i in 1:nrow(incidence.mat)) {
    idx <- incidence.mat[i, ] == 1
    print(paste("Processing group:", i))
    print("Index of columns being added:")
    print(which(idx))
    print("Number of columns being added:")
    print(sum(idx))
    print("Dimensions of X.latent before binding:")
    print(dim(X.latent))
    X.latent <- cbind(X.latent, X[, idx, drop = FALSE])
    print("Dimensions of X.latent after binding:")
    print(dim(X.latent))
    names <- c(names, colnames(incidence.mat)[idx])
    print("Current length of names:")
    print(length(names))
  }
  print("Length of final column names before assignment:")
  print(length(names))
  print("Dimensions of X.latent before final assignment:")
  print(dim(X.latent))
  colnames(X.latent) <- paste("grp", grp.vec, "_", names, sep = "")
  print("Final dimensions of X.latent:")
  print(dim(X.latent))
  print("Length of final column names:")
  print(length(colnames(X.latent)))
  X.latent
}
group <- readRDS("filtered_gene_sets.rds")
X.latent <- expandX(X, group)
incid.mat <- incidenceMatrix(X, group)
over.mat <- Matrix(incid.mat %*% t(incid.mat))
grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat))
run_asgl_pipeline <- function(X.latent, y, grp.vec, 
                              alphas = seq(0.05, 0.95, by = 0.05),
                              lambda_min = 0.02, tolerance = 1e-6, nfolds = 10) {
  # Ensure reproducibility
  set.seed(123456)
  # Loop through all alpha values
  for (alpha in alphas) {
    message("Running alpha = ", alpha, ", lambda_min = ", lambda_min)
    # Step 1: Fit SGL model
    fit <- asgl(X.latent, y, grp.vec, 
                family = "binomial", 
                alpha = alpha, 
                lambda_min = lambda_min)
    # Step 2: Cross-validation
    source("cv_asgl.R")
    cv_asgl_res <- cv_asgl(
      X.latent, y, grp.vec,
      family = "binomial",
      alpha = alpha,
      nfolds = nfolds,
      lambda = fit$lambda,
      lambda_min = lambda_min
    )
    # Step 3: Extract the lambda that minimizes CV error
    lambda.min <- cv_asgl_res$lambda
    lambda_min_index <- which(abs(fit$lambda - lambda.min) < tolerance)
    # Step 4: Extract coefficients
    coefficients <- fit$beta[, lambda_min_index, drop = FALSE]
    coefficients_df <- data.frame(coefficients)
    rownames(coefficients_df) <- colnames(X.latent)
    # Step 5: Identify positive and negative coefficients
    SGL_pos <- rownames(coefficients_df)[coefficients_df$coefficients > 0]
    SGL_neg <- rownames(coefficients_df)[coefficients_df$coefficients < 0]
    # Step 6: Save results
    file_stub <- paste0("(", alpha, "_", lambda_min, ")")
    write.table(SGL_pos, file = paste0("SGL_pos", file_stub, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(SGL_neg, file = paste0("SGL_neg", file_stub, ".txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(coefficients_df, file = paste0("coefficients_df", file_stub, ".txt"),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    # Step 7: Sanity check
    if (nrow(coefficients_df) != ncol(X.latent)) {
      stop(paste("Mismatch: coefficient rows != X.latent columns! Alpha =", alpha))
    }
  }
  message("Pipeline completed for all alpha values ✅")
}
run_asgl_pipeline(X.latent, y, grp.vec)


