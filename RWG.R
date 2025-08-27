###RWG
library(glmnet) 
library(Matrix)
exprSet <- read.table(file = "TCGA-TPM(normalize).txt", sep = "\t", header = T, row.names = 1)
X = t(exprSet)
filtered_gene_sets <- readRDS("filtered_gene_sets.rds")
genes_in_group <- unique(unlist(filtered_gene_sets))
common_genes <- intersect(colnames(X), genes_in_group)  
X_selected <- X[, common_genes]  
print(dim(X_selected))
y <- read.table("phenotype.txt", header = TRUE, sep = "\t")
Y <- as.numeric(as.character(y$group))  
alphas <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
for (alpha in alphas) {
  set.seed(123456)  
  fit <- glmnet(X_selected, Y, alpha = alpha, family = "binomial")  
  cv.fit <- cv.glmnet(X_selected, Y, alpha = alpha, family = "binomial", nfolds = 10)
  # Extract coefficients
  coefficients <- coef(cv.fit, s = "lambda.min")
  coefficients <- as.matrix(coefficients)   # <-- important fix
  coefficients_no_intercept <- coefficients[-1, , drop = FALSE]
  # Positive and negative features
  selected_pos <- rownames(coefficients_no_intercept)[coefficients_no_intercept > 0]
  selected_neg <- rownames(coefficients_no_intercept)[coefficients_no_intercept < 0]
  # Save results
  write.table(selected_pos, file = paste0("selected_pos(", alpha, ").txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(selected_neg, file = paste0("selected_neg(", alpha, ").txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(coefficients_no_intercept, file = paste0("coefficients_df(", alpha, ").txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  print(paste("Completed: alpha =", alpha))
}
