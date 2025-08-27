###Figure 4A, C, E
library(ggplot2)
library(patchwork)
library(grid)
library(ggtext)
library(ggsignif)
library(stringr)
library(dplyr)
# ---------------------------
# Main function
# ---------------------------
plot_ICI_response <- function(data_file,
                              params,
                              output_file = "ICI_response_vs_PBIS_TP53.pdf") {
  # Load data
  data <- read.table(data_file, header = TRUE, sep = "\t")
  # ---------------------------
  # helper: automatically choose test
  # ---------------------------
  get_p_value <- function(df) {
    x <- na.omit(df$PBIS.TP53[df$Response == "CR/PR"])
    y <- na.omit(df$PBIS.TP53[df$Response == "SD/PD"])
    safe_wilcox <- function(x, y) {
      tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    }
    if (length(x) < 3 || length(y) < 3 || sd(x) == 0 || sd(y) == 0) {
      return(safe_wilcox(x, y))
    }
    p_norm_x <- tryCatch(shapiro.test(x)$p.value, error = function(e) 0)
    p_norm_y <- tryCatch(shapiro.test(y)$p.value, error = function(e) 0)
    if (is.finite(p_norm_x) && is.finite(p_norm_y) &&
        p_norm_x > 0.05 && p_norm_y > 0.05) {
      p_var <- tryCatch(var.test(x, y)$p.value, error = function(e) 0)
      if (is.finite(p_var) && p_var > 0.05) {
        return(t.test(x, y, var.equal = TRUE)$p.value)
      } else {
        return(t.test(x, y, var.equal = FALSE)$p.value)
      }
    } else {
      return(safe_wilcox(x, y))
    }
  }
  # ---------------------------
  # helper: single violin plot
  # ---------------------------
  create_violin_plot <- function(data, age_group, title,
                                 y_limits = NULL, y_breaks = NULL,
                                 p_position = NULL,
                                 fill_color = "#A184BC",
                                 bg_colors = c("white", "#EFEFEF", "#DDDDDD")) {
    subset_data <- subset(data, Cohort == age_group)
    subset_data$Response <- factor(subset_data$Response,
                                   levels = c("CR/PR", "SD/PD"))
    p_val <- get_p_value(subset_data)
    p_label <- paste0("italic(P)==", formatC(p_val, format = "e", digits = 2))
    # auto y-axis
    if (is.null(y_limits)) {
      rng <- range(subset_data$PBIS.TP53, na.rm = TRUE)
      pad <- diff(rng) * 0.08
      if (!is.finite(pad) || pad == 0) pad <- 0.5
      y_limits <- c(rng[1] - pad, rng[2] + pad)
    }
    if (is.null(y_breaks)) {
      y_breaks <- pretty(y_limits, n = 10)
    }
    if (is.null(p_position)) {
      p_position <- max(y_limits) - 0.1 * diff(range(y_limits))
    }
    gradient_grob <- rasterGrob(
      colorRampPalette(bg_colors)(256),
      width = unit(1, "npc"), height = unit(1, "npc"),
      interpolate = TRUE
    )
    ggplot(subset_data, aes(x = Response, y = PBIS.TP53, fill = Response)) +
      annotation_custom(gradient_grob,
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      geom_violin(trim = FALSE) +
      stat_summary(fun.data = "mean_sdl",
                   fun.args = list(mult = 1),
                   geom = "crossbar", width = 0.3, size = 0.3) +
      geom_signif(comparisons = list(c("CR/PR", "SD/PD")),
                  annotations = p_label, y_position = p_position,
                  textsize = 4, tip_length = 0, parse = TRUE) +
      scale_fill_manual(values = c("#9C9C9C", fill_color)) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) +
      labs(title = title, x = NULL, y = "PBIS.TP53") +
      theme_classic(base_size = 12) +
      theme(
        plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_textbox_simple(size = 11, color = "white",
                                            halign = 0.5, fill = fill_color,
                                            width = 1.2,
                                            padding = margin(3, 0, 3, 0),
                                            margin = margin(0, 0, 10, 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none"
      )
  }
  # ---------------------------
  # generate plots
  # ---------------------------
  plots <- lapply(params, function(p) {
    create_violin_plot(
      data = data,
      age_group = p$age_group,
      title = p$title,
      y_limits = p$y_limits,
      y_breaks = p$y_breaks,
      p_position = p$p_position,
      fill_color = p$fill_color,
      bg_colors = p$bg_colors
    )
  })
  # auto layout for multiple plots
  n_plot <- length(plots)
  ncol_auto <- if (n_plot <= 3) n_plot else if (n_plot <= 6) 2 else 3
  nrow_auto <- ceiling(n_plot / ncol_auto)
  pdf_width <- max(4, 4 * ncol_auto)
  pdf_height <- max(4, 4.8 * nrow_auto)
  pdf(file = output_file, width = pdf_width, height = pdf_height)
  combined_plot <- wrap_plots(lapply(plots, wrap_elements), ncol = ncol_auto)
  print(combined_plot)
  dev.off()
  return(combined_plot)
}
params <- list(
  list(age_group = "GSE78220", title = "SKCM<br>(anti-PD-1)",
       y_limits = c(-1.5, 0.9), y_breaks = seq(-1.5, 0.9, 0.2),
       p_position = 0.6, fill_color = "#A184BC",
       bg_colors = c("white", "#E8DFF0", "#BFA0CC")),
  list(age_group = "IMvigor210", title = "UC<br>(anti-PD-L1)",
       y_limits = c(-1.3, 1.0), y_breaks = seq(-1.3, 1.0, 0.2),
       p_position = 0.7, fill_color = "#3E885B",
       bg_colors = c("white", "#E2F0E7", "#A6D8B2")),
  list(age_group = "GSE91061", title = "SKCM<br>(anti-PD-1, GSE91061)",
       y_limits = c(-1.3, 0.5), y_breaks = seq(-1.3, 0.5, 0.2),
       p_position = 0.3, fill_color = "#D48265",
       bg_colors = c("white", "#F8EAE3", "#E9CABA"))
)
plot_ICI_response(
  data_file = "ICIs应答与PBIS.TP53的关系汇总.txt",
  params = params,
  output_file = "Figure 4A,C,E.pdf"
)

###Figure 4B
library(pROC)
plot_multi_roc <- function(data,
                           response_col = "Response",
                           response_levels = c("CR/PR", "SD/PD"),
                           predictors,
                           colors,
                           output_file = "multi_ROC.pdf",
                           width = 5,
                           height = 5) {
  # Convert response to binary (0/1)
  data$Group <- ifelse(data[[response_col]] == response_levels[1], 0, 1)
  # Compute ROC curves for each predictor
  rocs <- lapply(predictors, function(m) pROC::roc(data$Group, data[[m]], na.rm = TRUE))
  names(rocs) <- predictors
  # Decide which predictors are suitable for smoothing (>=5 unique values)
  can_smooth_vec <- sapply(predictors, function(m) {
    x <- data[[m]]
    ux <- unique(x[!is.na(x)])
    length(ux) >= 5
  })
  # Pick base ROC curve (prefer the first smoothable, else the first one)
  base_idx <- if (any(can_smooth_vec)) which(can_smooth_vec)[1] else 1
  # Open PDF device
  pdf(file = output_file, width = width, height = height)
  # Plot base curve
  if (can_smooth_vec[base_idx]) {
    plot(smooth(rocs[[base_idx]]), col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } else {
    plot(rocs[[base_idx]], col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  }
  # Overlay the rest
  for (i in seq_along(rocs)) {
    if (i == base_idx) next
    if (can_smooth_vec[i]) {
      ok <- try({
        plot(smooth(rocs[[i]]), col = colors[i], lwd = 2, add = TRUE)
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
      }
    } else {
      plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
    }
  }
  # Legend with AUC values
  auc_labels <- sprintf("%s (AUC: %.2f)", predictors,
                        sapply(rocs, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", cex = 1.0, legend = auc_labels,
         col = colors, lty = 1, lwd = 2)
  dev.off()
  message("ROC curves saved to: ", output_file)
}
data <- read.table("GSE78220_PBIS-TP53.txt", header = TRUE, sep = "\t")
measures <- c(
  "PBIS.TP53", "IMPRES", "IPS",
  "CTLA4", "ImmuCellAI", "TMB",
  "Age", "TIDE", "CD274", "MIAS",
  "PDCD1", "M.Stage", "Gender", "GEP"
)
cols <- c(
  "#0072B2", "#009E73", "#56B4E9", "#FF8000", "#7030A0", "#E69F00",
  "#008000", "#00B050", "#C00000", "#F0E442", "#999999", "#7F7F7F", "#8B0000", "#D55E00"
)
plot_multi_roc(
  data = data,
  response_col = "Response",
  response_levels = c("CR/PR", "SD/PD"),
  predictors = measures,
  colors = cols,
  output_file = "Figure 4B.pdf"
)

###Figure 4D
library(pROC)
plot_multi_roc <- function(data,
                           response_col = "Response",
                           response_levels = c("CR/PR", "SD/PD"),
                           predictors,
                           colors,
                           output_file = "multi_ROC.pdf",
                           width = 5,
                           height = 5) {
  # Convert response into binary (0/1)
  data$Group <- ifelse(data[[response_col]] == response_levels[1], 0, 1)
  # Compute ROC curves
  rocs <- lapply(predictors, function(m) pROC::roc(data$Group, data[[m]], na.rm = TRUE))
  names(rocs) <- predictors
  # Decide which predictors can be smoothed (>=5 unique values)
  can_smooth_vec <- sapply(predictors, function(m) {
    x <- data[[m]]
    ux <- unique(x[!is.na(x)])
    length(ux) >= 5
  })
  # Pick a base index for the first plot
  base_idx <- if (any(can_smooth_vec)) which(can_smooth_vec)[1] else 1
  # Open PDF
  pdf(file = output_file, width = width, height = height)
  # Draw base ROC
  if (can_smooth_vec[base_idx]) {
    plot(smooth(rocs[[base_idx]]), col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } else {
    plot(rocs[[base_idx]], col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } 
  # Add the rest
  for (i in seq_along(rocs)) {
    if (i == base_idx) next
    if (can_smooth_vec[i]) {
      ok <- try({
        plot(smooth(rocs[[i]]), col = colors[i], lwd = 2, add = TRUE)
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
      }
    } else {
      plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
    }
  }
  # Add legend with AUC values
  auc_labels <- sprintf("%s (AUC: %.2f)", predictors,
                        sapply(rocs, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", cex = 1.0, legend = auc_labels,
         col = colors, lty = 1, lwd = 2)
  dev.off()
  message("ROC curves saved to: ", output_file)
  invisible(rocs)
}
data <- read.table("IMvigor210_PBIS-TP53.txt", header = TRUE, sep = "\t")
measures <- c(
  "PBIS.TP53","ImmuCellAI","IPS", "CD274", "GEP", "PDCD1",
  "TIDE", "MIAS","CTLA4","IMPRES","Gender"
)
cols <- c(
  "#0072B2", "#7030A0", "#56B4E9", "#C00000","#D55E00",
  "#999999","#00B050", "#F0E442", "#FF8000", "#009E73", "#8B0000"
)
plot_multi_roc(
  data = data,
  response_col = "Response",
  response_levels = c("CR/PR", "SD/PD"),
  predictors = measures,
  colors = cols,
  output_file = "Figure 4D.pdf"
)

###Figure 4F
library(pROC)
plot_multi_roc <- function(data,
                           response_col = "Response",
                           response_levels = c("CR/PR", "SD/PD"),
                           predictors,
                           colors,
                           output_file = "multi_ROC.pdf",
                           width = 5,
                           height = 5) {
  # Convert response into binary (0/1)
  data$Group <- ifelse(data[[response_col]] == response_levels[1], 0, 1)
  # Compute ROC curves
  rocs <- lapply(predictors, function(m) pROC::roc(data$Group, data[[m]], na.rm = TRUE))
  names(rocs) <- predictors
  # Decide which predictors can be smoothed (>=5 unique values)
  can_smooth_vec <- sapply(predictors, function(m) {
    x <- data[[m]]
    ux <- unique(x[!is.na(x)])
    length(ux) >= 5
  })
  # Pick a base index for the first plot
  base_idx <- if (any(can_smooth_vec)) which(can_smooth_vec)[1] else 1
  # Open PDF
  pdf(file = output_file, width = width, height = height)
  # Draw base ROC
  if (can_smooth_vec[base_idx]) {
    plot(smooth(rocs[[base_idx]]), col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } else {
    plot(rocs[[base_idx]], col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  }
  # Overlay remaining curves
  for (i in seq_along(rocs)) {
    if (i == base_idx) next
    if (can_smooth_vec[i]) {
      ok <- try({
        plot(smooth(rocs[[i]]), col = colors[i], lwd = 2, add = TRUE)
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
      }
    } else {
      plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
    }
  }
  # Add legend with AUC values
  auc_labels <- sprintf("%s (AUC: %.2f)", predictors,
                        sapply(rocs, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", cex = 1.0, legend = auc_labels,
         col = colors, lty = 1, lwd = 2)
  dev.off()
  message("ROC curves saved to: ", output_file)
  invisible(rocs)
}
data <- read.table("GSE91061_PBIS-TP53.txt", header = TRUE, sep = "\t")
measures <- c(
  "PBIS.TP53","TIDE","PDCD1","IMPRES","ImmuCellAI","CTLA4",
  "MIAS","GEP","CD274","TMB","IPS","M.Stage"
)
cols <- c(
  "#0072B2","#00B050","#999999","#009E73","#7030A0",
  "#FF8000","#C00000","#D55E00","#F0E442","#E69F00",
  "#56B4E9","#7F7F7F"
)
plot_multi_roc(
  data = data,
  response_col = "Response",
  response_levels = c("CR/PR", "SD/PD"),
  predictors = measures,
  colors = cols,
  output_file = "Figure 4F.pdf"
)

###Figure 4G, I
library(ggplot2)
library(ggpubr)
library(dplyr)
plot_group_box <- function(data_file,
                           x_col = "XG",
                           y_col = "value",
                           group_col = "group",
                           levels = c("G1", "G2"),
                           colors = c("G1" = "#47cf73", "G2" = "#C00000"),
                           output_file = "group_comparison.pdf",
                           width = 6,
                           height = 4) {
  # Load data
  df <- read.table(data_file, header = TRUE, sep = "\t")
  # Auto test function
  get_p_value <- function(x, y) {
    p_norm_x <- shapiro.test(x)$p.value
    p_norm_y <- shapiro.test(y)$p.value
    p_var <- var.test(x, y)$p.value
    if (p_norm_x > 0.05 && p_norm_y > 0.05) {
      if (p_var > 0.05) {
        test <- t.test(x, y, var.equal = TRUE)
      } else {
        test <- t.test(x, y, var.equal = FALSE)
      }
    } else {
      test <- wilcox.test(x, y)
    }
    return(test$p.value)
  }
  # Compute p-values for each group
  pval_df <- df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      p = get_p_value(
        .data[[y_col]][.data[[x_col]] == levels[1]],
        .data[[y_col]][.data[[x_col]] == levels[2]]
      ),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("p = ", signif(p, 3)))
  # Merge back for annotation
  df_plot <- left_join(df, pval_df, by = setNames("group", group_col))
  # Open PDF
  pdf(file = output_file, width = width, height = height)
  # Plot
  p <- ggplot(df_plot, aes_string(x = x_col, y = y_col, fill = x_col)) +
    geom_boxplot(linewidth = 0.6) +
    geom_point(shape = 21, size = 2, position = position_jitter(width = 0.15)) +
    facet_grid(reformulate(group_col, response = ".")) +
    geom_text(
      data = distinct(df_plot, !!sym(group_col), label),
      aes(x = 1.5, y = max(df[[y_col]], na.rm = TRUE) + 0.15, label = label),
      inherit.aes = FALSE,
      size = 3.5,
      color = "black"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10, color = "#204056"),
      axis.text.x = element_text(size = 10, color = "#204056"),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_fill_manual(values = colors)
  print(p)
  dev.off()
  message("Plot saved to: ", output_file)
  return(p)
}
plot_group_box(
  data_file = "chemo.txt",
  x_col = "XG",
  y_col = "value",
  group_col = "group",
  levels = c("G1", "G2"),
  colors = c("G1" = "#47cf73", "G2" = "#C00000"),
  output_file = "Figure 4G,I.pdf"
)

###Figure 4H
library(pROC)
plot_multi_roc <- function(data_file,
                           response_col = "Response",
                           response_levels = c("CR/PR", "SD/PD"),
                           predictors,
                           colors,
                           output_file = "multi_ROC.pdf",
                           width = 5,
                           height = 5) {
  # Load data
  data <- read.table(data_file, header = TRUE, sep = "\t")
  # Convert response into binary (0/1)
  data$Group <- ifelse(data[[response_col]] == response_levels[1], 0, 1)
  # Compute ROC curves
  rocs <- lapply(predictors, function(m) pROC::roc(data$Group, data[[m]], na.rm = TRUE))
  names(rocs) <- predictors
  # Decide which predictors can be smoothed (>=5 unique values)
  can_smooth_vec <- sapply(predictors, function(m) {
    x <- data[[m]]
    ux <- unique(x[!is.na(x)])
    length(ux) >= 5
  })
  # Pick a base index (first smoothable, else first predictor)
  base_idx <- if (any(can_smooth_vec)) which(can_smooth_vec)[1] else 1
  # Open PDF
  pdf(file = output_file, width = width, height = height)
  # Draw base ROC curve
  if (can_smooth_vec[base_idx]) {
    plot(smooth(rocs[[base_idx]]), col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } else {
    plot(rocs[[base_idx]], col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  }
  # Overlay the rest
  for (i in seq_along(rocs)) {
    if (i == base_idx) next
    if (can_smooth_vec[i]) {
      ok <- try({
        plot(smooth(rocs[[i]]), col = colors[i], lwd = 2, add = TRUE)
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
      }
    } else {
      plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
    }
  }
  # Add legend with AUC values
  auc_labels <- sprintf("%s (AUC: %.2f)", predictors,
                        sapply(rocs, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", cex = 1.0, legend = auc_labels,
         col = colors, lty = 1, lwd = 2)
  dev.off()
  message("ROC curves saved to: ", output_file)
  invisible(rocs)
}
measures <- c("PBIS.TP53","age","stage","tumor.grage")
cols <- c("#0072B2","#008000","#C00000","#00B050")
plot_multi_roc(
  data_file   = "GSE163882_PBIS-TP53.txt",
  response_col = "Response",
  response_levels = c("CR/PR", "SD/PD"),
  predictors  = measures,
  colors      = cols,
  output_file = "Figure 4H.pdf"
)

###Figure 4J
library(pROC)
plot_multi_roc <- function(data_file,
                           response_col = "Response",
                           response_levels = c("CR/PR", "SD/PD"),
                           predictors,
                           colors,
                           output_file = "multi_ROC.pdf",
                           width = 5,
                           height = 5) {
  # Load data
  data <- read.table(data_file, header = TRUE, sep = "\t")
  # Convert response to binary (0/1)
  data$Group <- ifelse(data[[response_col]] == response_levels[1], 0, 1)
  # Compute ROC curves
  rocs <- lapply(predictors, function(m) pROC::roc(data$Group, data[[m]], na.rm = TRUE))
  names(rocs) <- predictors
  # Decide which predictors can be smoothed (>=5 unique values)
  can_smooth_vec <- sapply(predictors, function(m) {
    x <- data[[m]]
    ux <- unique(x[!is.na(x)])
    length(ux) >= 5
  })
  # Pick a base index (first smoothable, else first predictor)
  base_idx <- if (any(can_smooth_vec)) which(can_smooth_vec)[1] else 1
  # Open PDF
  pdf(file = output_file, width = width, height = height)
  # Draw base ROC curve
  if (can_smooth_vec[base_idx]) {
    plot(smooth(rocs[[base_idx]]), col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  } else {
    plot(rocs[[base_idx]], col = colors[base_idx], lwd = 2, legacy.axes = TRUE)
  }
  # Overlay the rest
  for (i in seq_along(rocs)) {
    if (i == base_idx) next
    if (can_smooth_vec[i]) {
      ok <- try({
        plot(smooth(rocs[[i]]), col = colors[i], lwd = 2, add = TRUE)
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
      }
    } else {
      plot(rocs[[i]], col = colors[i], lwd = 2, add = TRUE)
    }
  }
  # Add legend with AUC values
  auc_labels <- sprintf("%s (AUC: %.2f)", predictors,
                        sapply(rocs, function(r) as.numeric(pROC::auc(r))))
  legend("bottomright", cex = 1.0, legend = auc_labels,
         col = colors, lty = 1, lwd = 2)
  
  dev.off()
  message("ROC curves saved to: ", output_file)
  invisible
measures <- c("PBIS.TP53","age","stage","tumor_purity_")
cols <- c("#0072B2","#008000","#C00000","#00B050")
plot_multi_roc(
  data_file   = "GSE123845_PBIS-TP53.txt",
  response_col = "Response",
  response_levels = c("CR/PR", "SD/PD"),
  predictors  = measures,
  colors      = cols,
  output_file = "Figure 4J.pdf"
)
