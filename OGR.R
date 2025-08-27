###OGR
library(grpreg)
library(Matrix)
exprSet <- read.table(file = "TCGA-TPM(normalize).txt", sep = "\t", header = T, row.names = 1)
X=t(exprSet)
y <- read.table("phenotype.txt", header = TRUE, sep = "\t")
Y <-as.numeric(as.character(y$group))
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
alphas <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
for (alpha in alphas) {
  set.seed(123456)
  fit <- grpreg(
    X.latent, Y, group = grp.vec,
    penalty = "grLasso", family = "binomial", alpha = alpha,lambda.min=0.02
  )
  cv.fit <- cv.grpreg(
    X.latent, Y, group = grp.vec, lambda = fit$lambda,
    penalty = "grLasso", 
    family = "binomial", alpha = alpha, nfolds = 10, seed = 123456,lambda.min=0.02
  )  
  error <- cv.fit$cve
  lambda_opt_index <- which.min(error)
  coefficients <- fit$beta[-1, lambda_opt_index] 
  coefficients_df <- as.data.frame(coefficients)
  selected_pos <- colnames(X.latent)[which(coefficients_df$coefficients > 0)]
  selected_neg <- colnames(X.latent)[which(coefficients_df$coefficients < 0)]
  write.table(selected_pos, file = paste0("selected_pos(", alpha, ").txt"), sep = "\t", quote = FALSE)
  write.table(selected_neg, file = paste0("selected_neg(", alpha, ").txt"), sep = "\t", quote = FALSE)
  write.table(coefficients_df, file = paste0("coefficients_df(", alpha, ").txt"), sep = "\t", quote = FALSE, row.names = TRUE)
  print(paste("Completed: alpha =", alpha))
}