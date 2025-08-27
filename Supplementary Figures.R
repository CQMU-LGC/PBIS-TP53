###Figure S1
library(ggplot2)
library(plyr)
library(ggord)
# Make sure you have sourced geom_ord_ellipse.R before running
# source("geom_ord_ellipse.R")
plot_pca_with_groups <- function(expr_file,
                                 pheno_file,
                                 sample_col = "sample",
                                 group_col = "Tumor_type",
                                 output_file = "PCA_plot.pdf",
                                 colors = NULL,
                                 width = 6,
                                 height = 6) {
  # --- Step 1: Load data
  expr_data <- read.table(expr_file, header = TRUE, sep = "\t")
  pheno <- read.table(pheno_file, header = TRUE, sep = "\t")
  # --- Step 2: Match sample order
  matched_indices <- match(pheno[[sample_col]], colnames(expr_data))
  expr_data <- expr_data[, matched_indices]
  if (!identical(pheno[[sample_col]], colnames(expr_data))) {
    stop("❌ Sample IDs in phenotype file and expression matrix do not match.")
  }
  # --- Step 3: PCA
  expr_data_t <- t(expr_data)
  pca.results <- prcomp(expr_data_t, center = TRUE, scale. = FALSE)
  # --- Step 4: Define colors
  if (is.null(colors)) {
    default_colors <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D",
                        "#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D",
                        "#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
    colors <- default_colors[1:length(unique(pheno[[group_col]]))]
  }
  # --- Step 5: Plot
  pdf(file = output_file, width = width, height = height, bg = "white")
  p <- ggord(pca.results,
             grp_in = pheno[[group_col]],
             repel = TRUE,
             ellipse = FALSE,
             size = 2,
             alpha = 0.5,
             cols = colors,
             arrow = NULL, txt = NULL) +
    theme(panel.grid = element_blank()) +
    geom_ord_ellipse(ellipse_pro = 0.95,
                     size = 1.5,
                     lty = 1)
  print(p)
  dev.off()
  message("✅ PCA plot saved to: ", output_file)
  return(p)
}
# Before batch correction
plot_pca_with_groups(expr_file = "TCGA-TPM.txt",
                     pheno_file = "phenotype.txt",
                     output_file = "Figure S1A.pdf")
# After batch correction
plot_pca_with_groups(expr_file = "TCGA-TPM(normalize).txt",
                     pheno_file = "phenotype.txt",
                     output_file = "Figure S1B.pdf")


###Figure S2 (Figures S2B–L are interpreted in the same manner)
library(meta)
run_forest_meta <- function(data_file,
                            output_file = "forest_plot.pdf",
                            hr_col = "HR",
                            lci_col = "LCI",
                            uci_col = "UCI",
                            dataset_col = "Dataset",
                            width = 8,
                            height = 5) {
  # --- Step 1: Load data
  exprSet <- read.table(file = data_file, sep = "\t", header = TRUE)
  # --- Step 2: Meta-analysis
  mg <- metagen(log(exprSet[[hr_col]]),
                lower = log(exprSet[[lci_col]]),
                upper = log(exprSet[[uci_col]]),
                data = exprSet,
                sm = "HR",
                studlab = paste(exprSet[[dataset_col]]))
  # --- Step 3: Print summary in console
  print(summary(mg))
  # --- Step 4: Save forest plot
  pdf(output_file, width = width, height = height)
  forest(mg,
         common = FALSE,   # show only random-effects model
         layout = "JAMA")  # journal style
  dev.off()
  message("✅ Forest plot saved to: ", output_file)
  return(mg)
}
meta_result <- run_forest_meta(
  data_file = "ACC.txt",
  output_file = "Figure S2A.pdf"
)

###Figures S3A–C were generated using code analogous to that for Figure 8I.

###Figure S3C
run_meta_forest <- function(data_file,
                            output_file = "forest_plot.pdf",
                            hr_col = "HR",
                            lci_col = "LCI",
                            uci_col = "UCI",
                            dataset_col = "Dataset",
                            width = 8,
                            height = 5,
                            layout = "RevMan5",
                            show_common = FALSE) {
  # --- Step 1: Load data
  exprSet <- read.table(file = data_file, sep = "\t", header = TRUE)
  # --- Step 2: Meta-analysis
  mg <- metagen(log(exprSet[[hr_col]]),
                lower = log(exprSet[[lci_col]]),
                upper = log(exprSet[[uci_col]]),
                data = exprSet,
                sm = "HR",
                studlab = paste(exprSet[[dataset_col]]))
  # --- Step 3: Summary
  print(summary(mg))
  # --- Step 4: Forest plot
  pdf(output_file, width = width, height = height)
  forest(mg,
         common = show_common,   # FALSE = hide fixed-effect, TRUE = show
         layout = layout)
  dev.off()
  message("✅ Forest plot saved to: ", output_file)
  return(mg)
}
meta_result <- run_meta_forest(
  data_file = "ACC.txt",
  output_file = "Figure S3C.pdf",
  layout = "RevMan5",
  show_common = FALSE
)
###Figures S4–14 were generated using code analogous to that for Figure S3.

###Figures S15–18 were generated using the following code as reference.
library(ComplexHeatmap)
library(circlize)
library(grid)
plot_weighted_heatmap <- function(input_file,
                                  output_file = "heatmap.pdf",
                                  bar_title = "Average C-index",
                                  heatmap_title = "Weight Average C-index",
                                  bar_width = 4,
                                  pdf_width = 14,
                                  pdf_height = 30) {
  # --- Step 1: Load data
  statical_mat <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1)
  # --- Step 2: Extract Weight_average
  avg_statical <- statical_mat$Weight_average
  names(avg_statical) <- rownames(statical_mat)
  # --- Step 3: Custom sort (descending, "Ours" at top if present)
  statical_mat <- statical_mat[order(-avg_statical, rownames(statical_mat) != "Ours"), ]
  avg_statical <- avg_statical[order(avg_statical, decreasing = TRUE)]
  # --- Step 4: Color mapping for barplot
  bar_col_fun <- colorRamp2(
    c(min(avg_statical), mean(avg_statical), max(avg_statical)),
    c("blue", "white", "red")
  )
  # --- Step 5: Barplot annotation
  right_anno <- rowAnnotation(
    "{bar_title}" := anno_barplot(
      avg_statical,
      gp = gpar(fill = bar_col_fun(avg_statical), col = "grey"),
      border = TRUE,
      width = unit(bar_width, "cm"),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 10, col = "black")
    )
  )
  # --- Step 6: Heatmap body (exclude Weight_average col)
  heatmap_data <- statical_mat[, -ncol(statical_mat)]
  hm <- Heatmap(
    heatmap_data,
    name = heatmap_title,
    col = colorRamp2(c(0.45, 0.6, 0.75), c("#66C2A5", "grey95", "#FF6EC7")),
    rect_gp = gpar(col = "grey"),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(col = "black"),
    width = unit(ncol(heatmap_data) + 3, "cm"),
    height = unit(nrow(heatmap_data) / 1.5, "cm"),
    heatmap_legend_param = list(title = "C-index"),
    right_annotation = right_anno,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 10))
    }
  )
  # --- Step 7: Save PDF
  pdf(output_file, width = pdf_width, height = pdf_height)
  draw(hm)
  dev.off()
  message("✅ Heatmap saved to: ", output_file)
  return(hm)
}
plot_weighted_heatmap(
  input_file = "LUAD_C-index.txt",
  output_file = "LUAD_heatmap_WeightAverageC-index.pdf",
  bar_title = "Average C-index",
  heatmap_title = "Weight Average C-index",
  pdf_width = 14,
  pdf_height = 30
)

###Figures S19–28 were generated using the following code as reference.
run_cox_nomogram_analysis <- function(data,
                                      time_col = "OS.time",
                                      event_col = "OS",
                                      variables,
                                      ref_levels = list(
                                        Type = "Wild",
                                        gender.demographic = "male",
                                        ajcc_pathologic_stage.diagnoses = "I-II"
                                      ),
                                      output_prefix = "TCGA-HCC") {
  # ---- Load required packages ----
  require(dplyr)
  require(survival)
  require(rms)
  require(regplot)
  require(ggplot2)
  # ---- Relevel factors according to reference ----
  for (var in names(ref_levels)) {
    if (var %in% names(data)) {
      data[[var]] <- relevel(as.factor(data[[var]]), ref = ref_levels[[var]])
    }
  }
  # ---- Helper function: univariate Cox ----
  perform_univariate_cox <- function(data, variable) {
    formula <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~", variable))
    cox_model <- coxph(formula, data = data)
    summary_cox <- summary(cox_model)
    data.frame(
      Variable = variable,
      PValue   = summary_cox$coefficients[,"Pr(>|z|)"],
      HR       = summary_cox$coefficients[,"exp(coef)"],
      Lower95  = summary_cox$conf.int[,"lower .95"],
      Upper95  = summary_cox$conf.int[,"upper .95"]
    )
  }
  # ---- Run univariate Cox for each variable ----
  cox_results <- lapply(variables, function(var) perform_univariate_cox(data, var))
  cox_results_df <- do.call(rbind, cox_results)
  print(cox_results_df)
  # ---- Forest plot (univariate Cox) ----
  dat <- cox_results_df
  dat$color <- ifelse(dat$HR > 1 & dat$PValue < 0.05, "pink",
                      ifelse(dat$HR > 1, "grey",
                             ifelse(dat$PValue < 0.05, "green", "grey")))
  dat$Variable <- factor(dat$Variable, levels = rev(dat$Variable))
  p <- ggplot(dat, aes(x = Variable, y = HR, ymin = Lower95, ymax = Upper95)) +
    geom_pointrange(fatten = 1.5, size = 1, aes(col = color)) +
    geom_hline(yintercept = 1, lty = 2, size = 1, color = "black") +
    coord_flip() +
    scale_color_manual(values = c("pink" = "pink", "grey" = "grey", "green" = "green")) +
    ylab("Hazard Ratio (HR)") +
    scale_x_discrete(expand = expansion(add = 1.5)) +
    theme_bw(base_size = 16) +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 16))

  pdf(paste0("UnivariateCox-ForestPlot(", output_prefix, ").pdf"), width = 6, height = 5)
  print(p)
  dev.off()
  # ---- Select significant variables ----
  significant_vars <- dat$Variable[dat$PValue < 0.05]
  selected_cols <- c(time_col, event_col, as.character(significant_vars))
  filtered_data <- na.omit(data[, selected_cols, drop = FALSE])
  # ---- Multivariate Cox ----
  formula <- as.formula(paste("Surv(", time_col, ",", event_col, ") ~", 
                              paste(significant_vars, collapse = " + ")))
  res.cox <- coxph(formula, data = filtered_data)
  print(summary(res.cox))
  # ---- Nomogram ----
  pdf(paste0("Nomogram(", output_prefix, ").pdf"), width = 7, height = 5)
  regplot(res.cox,
          plots = c("spikes", "spikes"),
          spkcol = "#F781BF",
          clickable = FALSE,
          title = "Nomogram",
          points = TRUE,
          droplines = TRUE,
          observation = filtered_data[1, ],
          rank = "sd",
          failtime = c(365, 1095, 1825),
          prfail = FALSE)
  dev.off()
  # ---- Calibration curves ----
  filtered_data$risk <- predict(res.cox, newdata = filtered_data, type = "risk")
  pdf(paste0("Calibration(", output_prefix, ").pdf"), width = 5, height = 5)
  times <- c(365, 1095, 1825)
  cols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")
  for (i in seq_along(times)) {
    f <- cph(Surv(get(time_col), get(event_col)) ~ risk, 
             x = TRUE, y = TRUE, surv = TRUE, 
             data = filtered_data, time.inc = times[i])
    cal <- calibrate(f, cmethod = "KM", method = "boot", u = times[i], 
                     m = (nrow(filtered_data)/3), B = 500)
    plot(cal, xlim = c(0, 1), ylim = c(0, 1),
         xlab = ifelse(i == 1, "Nomogram-predicted OS (%)", ""), 
         ylab = ifelse(i == 1, "Observed OS (%)", ""), 
         lwd = 3, col = cols[i], sub = FALSE, add = i != 1)
  }
  legend("bottomright", c("1-year", "3-year", "5-year"), col = cols, lwd = 3, bty = "n")
  dev.off()
  return(list(univariate = cox_results_df, multivariate = summary(res.cox)))
}
results <- run_cox_nomogram_analysis(
  data = ii,
  variables = c("Score", "Type", "gender.demographic", "age_at_index.demographic", "ajcc_pathologic_stage.diagnoses"),
  output_prefix = "TCGA-HCC"
)

###Figure S29A
VlnPlot(merged_seurat, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()

###Figure S29B
VlnPlot(merged_seurat_filtered, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()

###Figure S30A
plot_umap_with_ellipses <- function(scobj,
                                    reduction = "pca",
                                    dims = 1:20,
                                    reduction_name = "umap_custom",
                                    palette = c4a("bold", 9),
                                    point_size = 0.4,
                                    point_alpha = 0.8,
                                    label_size = 4,
                                    legend_point_size = 5,
                                    ellipse_alpha = 0.1,
                                    arrow_length = 0.2) {
  # ---- Run UMAP ----
  scobj <- RunUMAP(scobj,
                   reduction = reduction,
                   dims = dims,
                   reduction.name = reduction_name)
  umap_df <- as.data.frame(scobj@reductions[[reduction_name]]@cell.embeddings)
  umap_df$orig.ident <- scobj@meta.data$orig.ident
  # ---- Base theme ----
  mytheme <- theme_void() +
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5))
  # ---- Base scatter ----
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = orig.ident), size = point_size, alpha = point_alpha)
  # ---- Confidence ellipses ----
  p <- p +
    stat_ellipse(aes(color = orig.ident, fill = orig.ident),
                 level = 0.95,
                 linetype = 1,
                 show.legend = FALSE,
                 geom = "polygon",
                 alpha = ellipse_alpha) +
    mytheme
  # ---- Axis arrows ----
  p <- p +
    theme_dr(xlength = arrow_length,
             ylength = arrow_length,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = "last", type = "closed")) +
    theme(panel.grid = element_blank())
  # ---- Median labels ----
  label_df <- umap_df %>%
    group_by(orig.ident) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
  p <- p +
    geom_text(data = label_df,
              aes(x = UMAP_1, y = UMAP_2, label = orig.ident),
              fontface = "bold",
              color = "black",
              size = label_size)
  # ---- Legend formatting ----
  p <- p +
    guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette)
  return(p)
}
library(Seurat)
library(ggplot2)
library(cols4all)
library(tidydr)
library(dplyr)
# Example call
p <- plot_umap_with_ellipses(scobj,
                             reduction = "pca",
                             dims = 1:20,
                             reduction_name = "umap_naive",
                             palette = c4a("pastel", 9))
print(p)

###Figure S30B
plot_umap_harmony <- function(scobj,
                              reduction = "harmony",
                              dims = 1:20,
                              reduction_name = "umap",
                              palette = c4a("bold", 9),
                              point_size = 0.4,
                              point_alpha = 0.8,
                              label_size = 4,
                              legend_point_size = 5,
                              ellipse_alpha = 0.1,
                              arrow_length = 0.2) {
  # ---- Run UMAP using harmony reduction ----
  scobj <- RunUMAP(scobj,
                   reduction = reduction,
                   dims = dims,
                   reduction.name = reduction_name)
  umap_df <- as.data.frame(scobj@reductions[[reduction_name]]@cell.embeddings)
  umap_df$orig.ident <- scobj@meta.data$orig.ident
  # ---- Custom theme ----
  mytheme <- theme_void() +
    theme(plot.margin = margin(5.5, 15, 5.5, 5.5))
  # ---- Base scatter ----
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = orig.ident), size = point_size, alpha = point_alpha)
  # ---- Confidence ellipses ----
  p <- p +
    stat_ellipse(aes(color = orig.ident, fill = orig.ident),
                 level = 0.95,
                 linetype = 1,
                 show.legend = FALSE,
                 geom = "polygon",
                 alpha = ellipse_alpha) +
    mytheme
  # ---- Axis arrows ----
  p <- p +
    theme_dr(xlength = arrow_length,
             ylength = arrow_length,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = "last", type = "closed")) +
    theme(panel.grid = element_blank())
  # ---- Median labels ----
  label_df <- umap_df %>%
    group_by(orig.ident) %>%
    summarise(UMAP_1 = median(UMAP_1),
              UMAP_2 = median(UMAP_2), .groups = "drop")
  p <- p +
    geom_text(data = label_df,
              aes(x = UMAP_1, y = UMAP_2, label = orig.ident),
              fontface = "bold",
              color = "black",
              size = label_size)
  # ---- Legend & colors ----
  p <- p +
    guides(color = guide_legend(override.aes = list(size = legend_point_size))) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette)
  return(p)
}
library(Seurat)
library(ggplot2)
library(cols4all)
library(tidydr)
library(dplyr)
# Run and plot UMAP using harmony reduction
p <- plot_umap_harmony(scobj,
                       reduction = "harmony",
                       dims = 1:20,
                       reduction_name = "umap",
                       palette = c4a("pastel", 9))
print(p)

###The code for Figure S31 can be referenced from the corresponding part of the code for Figure 5.


