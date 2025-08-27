###Figure 3A
library(ggplot2)
plot_cindex <- function(data_file,
                        title = "C-index by Alpha across Groups",
                        output_file = NULL) {
  # Step 1: Load the data
  df <- read.table(data_file, header = TRUE, sep = "\t", row.names = 1)
  # Step 2: Build the plot
  p <- ggplot(df, aes(x = factor(alpha), y = C_index, fill = Group)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = C_index),
              position = position_dodge(0.8),
              vjust = -0.5, size = 3) +
    theme_minimal() +
    labs(title = title,
         x = "Alpha", 
         y = "C-index") +
    theme(legend.position = "top")
  
  # Step 3: Save the plot if output_file is provided
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 6, height = 4)
  }
  return(p)
}
# Show the plot in R
plot_cindex("Cindex_comparison.txt")
# Save the plot as a PDF
plot_cindex("Cindex_comparison.txt",
            title = "C-index Comparison Across Methods",
            output_file = "Figure 3A.pdf")

###Figure 3B
plot_grouped_correlation <- function(data_file,
                                     x_var,
                                     y_var,
                                     class_var = "class",
                                     group_colors = c("OGR" = "#E51718",
                                                      "RWG" = "#1D2D60",
                                                      "OSGR" = "#239B56"),
                                     output_file = "grouped_correlation_plot.pdf") {
  # Load data
  cor.data <- read.table(data_file, header = TRUE, sep = "\t", row.names = 1)
  # Extract variables
  xvar <- cor.data[[x_var]]
  yvar <- cor.data[[y_var]]
  # Size scale based on absolute differences
  cor.data$diff <- abs(xvar - yvar)
  cor.data$size <- scales::rescale(cor.data$diff, to = c(1.2, 3))
  # Group levels
  group_levels <- intersect(names(group_colors), unique(cor.data[[class_var]]))
  # Axis limits with slight expansion
  xlim <- range(xvar, na.rm = TRUE)
  ylim <- range(yvar, na.rm = TRUE)
  xlim <- xlim + c(-0.05, 0.05) * diff(xlim)
  ylim <- ylim + c(-0.05, 0.05) * diff(ylim)
  # Start plotting
  pdf(output_file, width = 7, height = 6)
  par(bty = "o", mgp = c(2, 0.5, 0), mar = c(4.1, 4.1, 2.1, 4.6), tcl = -0.25, font.main = 3)
  par(xpd = FALSE)
  plot(NULL, NULL, xlim = xlim, ylim = ylim,
     xlab = x_var, ylab = y_var, col = "white", main = "")
  # Background
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "#EAE9E9", border = FALSE)
  grid(col = "white", lty = 1, lwd = 1.5)
  # Loop through groups
  y_text <- max(ylim) - 0.02 * diff(ylim)
  for (grp in group_levels) {
    tmp <- subset(cor.data, cor.data[[class_var]] == grp)
    if (nrow(tmp) == 0) next
    x <- tmp[[x_var]]
    y <- tmp[[y_var]]
    # Choose correlation method
    method <- if (shapiro.test(x)$p.value > 0.05 && shapiro.test(y)$p.value > 0.05) "pearson" else "spearman"
    rho_test <- cor.test(x, y, method = method)
    rho_val <- round(rho_test$estimate, 2)
    pval_raw <- rho_test$p.value
    pval_fmt <- format.pval(pval_raw, digits = 2, eps = 1e-3)
    stars <- if (pval_raw <= 0.001) "***" else if (pval_raw <= 0.01) "**" else if (pval_raw <= 0.05) "*" else "ns"
    # Regression line
    reg <- lm(y ~ x)
    x_seq <- seq(xlim[1], xlim[2], length.out = 200)
    y_pred <- predict(reg, newdata = data.frame(x = x_seq))
    lines(x_seq, y_pred, lwd = 2, col = group_colors[grp])
    # Text annotation
    text(quantile(xlim, 0.03), y_text, adj = 0,
         bquote(.(grp)*": N = "*.(nrow(tmp))*"; "*rho*" = "*.(rho_val)*"; "*italic(P)*" = "*.(pval_fmt)*" ("*.(stars)*"); method = "*.(method)),
         col = group_colors[grp], cex = 0.85)
    y_text <- y_text - 0.06 * diff(ylim)
  }
  # Rugs
  rug(xvar[xvar >= xlim[1] & xvar <= xlim[2]], side = 3, col = "black", lwd = 0.6)
  rug(yvar[yvar >= ylim[1] & yvar <= ylim[2]], side = 4, col = "black", lwd = 0.6)
  # Legend for point sizes (shift indicator)
  unique_cex <- sort(unique(round(cor.data$size, 1)))
  legend_y_pos <- seq(ylim[1] + 0.85 * diff(ylim), ylim[1] + 0.65 * diff(ylim), length.out = length(unique_cex))
  par(xpd = TRUE)
  points(x = rep(par("usr")[2] + 0.02 * diff(xlim), length(unique_cex)), y = legend_y_pos,
         pch = 19, cex = unique_cex, col = "black")
  text(x = par("usr")[2] + 0.05 * diff(xlim),
       y = c(legend_y_pos[1] + 0.05 * diff(ylim), legend_y_pos),
       labels = c("Absolute\nVertical\nShift", round(scales::rescale(unique_cex, to = range(cor.data$diff)), 1)),
       adj = 0, cex = 0.8)
  # Legend for group colors
  legend_group_y <- seq(ylim[2], ylim[2] - 0.05 * length(group_levels) * diff(ylim), length.out = length(group_levels))
  points(x = rep(par("usr")[2] + 0.02 * diff(xlim), length(group_levels)), y = legend_group_y,
         pch = 19, cex = 1.6, col = group_colors[group_levels])
  text(x = rep(par("usr")[2] + 0.05 * diff(xlim), length(group_levels)),
       y = legend_group_y, labels = group_levels, adj = 0, cex = 0.8)
  # Border finishing
  par(new = TRUE, bty = "o")
  plot(-1, -1, col = "white", xlim = xlim, ylim = ylim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  invisible(dev.off())
  message("Plot saved to: ", output_file)
}
plot_grouped_correlation(
  data_file = "p53.Independent.DNA.Damage.Response.txt",
  x_var = "p53.Independent.DNA.Damage.Response",
  y_var = "PBIS.TP53",
  class_var = "class",
  output_file = "Figure 3B(1).pdf"
)
plot_grouped_correlation(
  data_file = "TP53.Regulates.Transcription.of.DNA.Repair.Genes.txt",
  x_var = "TP53.Regulates.Transcription.of.DNA.Repair.Genes",
  y_var = "PBIS.TP53",
  class_var = "class",
  output_file = "Figure 3B(2).pdf"
)
plot_grouped_correlation(
  data_file = "TP53.Regulates.Transcription.of.Cell.Cycle.Genes.txt",
  x_var = "TP53.Regulates.Transcription.of.Cell.Cycle.Genes",
  y_var = "PBIS.TP53",
  class_var = "class",
  output_file = "Figure 3B(3).pdf"
)

###Figure 3C
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
plot_group_comparison <- function(data_file,
                                  x_var = "Tumor",
                                  y_var = "Score",
                                  group_var = "Type",
                                  group_labels = c("Wild", "Mutation"),
                                  output_file = "group_comparison.pdf") {
  # --- Step 1: Load data
  df <- read.table(data_file, header = TRUE, sep = "\t")
  # --- Step 2: helper function to choose test automatically
  get_p_value <- function(sub_df) {
    x <- sub_df[[y_var]][sub_df[[group_var]] == group_labels[1]]
    y <- sub_df[[y_var]][sub_df[[group_var]] == group_labels[2]]
    if (length(x) < 3 || length(y) < 3) return(NA)
    p_norm_x <- shapiro.test(x)$p.value
    p_norm_y <- shapiro.test(y)$p.value
    p_var <- var.test(x, y)$p.value
    if (p_norm_x > 0.05 && p_norm_y > 0.05) {
      test <- if (p_var > 0.05) t.test(x, y, var.equal = TRUE) else t.test(x, y)
    } else {
      test <- wilcox.test(x, y)
    }
    return(test$p.value)
  }
  # --- Step 3: compute p-values per tumor type
  p_df <- df %>%
    group_by(.data[[x_var]]) %>%
    summarise(p = get_p_value(cur_data_all()), .groups = "drop") %>%
    mutate(p_label = case_when(
      is.na(p)       ~ "NA",
      p < 0.0001     ~ "****",
      p < 0.001      ~ "***",
      p < 0.01       ~ "**",
      p < 0.05       ~ "*",
      TRUE           ~ "ns"
    ))
  # --- Step 4: y-axis positions for annotations
  y_pos <- df %>%
    group_by(.data[[x_var]]) %>%
    summarise(y.position = max(.data[[y_var]], na.rm = TRUE) * 1.05, .groups = "drop")
  # Merge y positions with p-values
  p_df <- left_join(p_df, y_pos, by = x_var)
  # --- Step 5: make plot
  p <- ggplot(df, aes_string(x = x_var, y = y_var, color = group_var, fill = group_var)) +
    stat_boxplot(geom = "errorbar",
                 size = 0.6,
                 width = 0.4,
                 linetype = "solid",
                 position = position_dodge(1.0),
                 color = "black") +
    geom_boxplot(outlier.colour = NA,
                 size = 0.6,
                 width = 0.4,
                 linetype = "solid",
                 position = position_dodge(1.0),
                 color = "black") +
    scale_fill_jama(name = group_var) +
    coord_cartesian(ylim = c(min(df[[y_var]], na.rm = TRUE),
                             max(df[[y_var]], na.rm = TRUE) * 1.1)) +
    labs(x = x_var, y = y_var) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      legend.position = "top"
    ) +
    geom_text(data = p_df,
              aes_string(x = x_var, y = "y.position", label = "p_label"),
              inherit.aes = FALSE,
              size = 4,
              color = "black")
  # --- Step 6: save to file
  ggsave(output_file, plot = p, width = 15, height = 5)
  return(p)
}
plot_group_comparison(
  data_file = "PBIS-TP53_vs_TP53mutation.txt",
  x_var = "Tumor",
  y_var = "Score",
  group_var = "Type",
  group_labels = c("Wild", "Mutation"),
  output_file = "Figure 3C.pdf"
)


