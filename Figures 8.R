###Figure 8A (Figures 8B is interpreted in the same manner.)
library(Seurat)
library(Matrix)
library(GSVA)
library(CellChat) 
library(patchwork)
library(SeuratData)
cellchat <- readRDS("cellchat_Low.rds")
groupSize<-as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count,vertex.weight=groupSize,
                 weight.scale=T,label.edge=F,
                 title.name="Number of interactions")
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength")

###Figure 8C (Figures 8D are interpreted in the same manner.)
cellchat_filtered <- cellchat
cellchat_filtered@net$prob[cellchat_filtered@net$prob <= 0.05] <- 0
comm_all <- subsetCommunication(cellchat_filtered, slot.name = "net")
sources <- c("Cycling cell","Epithelial cell","Myeloid cell","Lymphoid cell")
target  <- "Stromal cell"
comm_keep <- subset(comm_all, source %in% sources & target %in% target)
pairLR.keep <- unique(comm_keep[, c("ligand","receptor","interaction_name")])
colnames(pairLR.keep) <- c("ligand","receptor","interaction_name")
netVisual_bubble(
  cellchat_filtered,
  sources.use = sources,
  targets.use = target,
  pairLR.use  = pairLR.keep,
  remove.isolate = TRUE,
  angle.x = 45
)

###Figure 8H (Figures 8E-G are interpreted in the same manner.)
library(Seurat)
library(ggpubr)
library(dplyr)
library(ggplot2)
get_p_value_expr <- function(df, value_col = "expr", group_col = "grp",
                             g1 = "Low", g2 = "High") {
  stopifnot(all(c(value_col, group_col) %in% names(df)))
  x <- df[[value_col]][df[[group_col]] == g1]
  y <- df[[value_col]][df[[group_col]] == g2]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3 ||
      sd(x) == 0 || sd(y) == 0) {
    test <- tryCatch(wilcox.test(x, y, exact = FALSE),
                     error = function(e) list(p.value = NA_real_))
    return(list(p.value = test$p.value, method = "wilcox (fallback)"))
  }
  p_norm_x <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
  p_norm_y <- tryCatch(shapiro.test(y)$p.value, error = function(e) NA_real_)
  p_var    <- tryCatch(var.test(x, y)$p.value,   error = function(e) NA_real_)
  if (!is.na(p_norm_x) && !is.na(p_norm_y) &&
      p_norm_x > 0.05 && p_norm_y > 0.05) {
    if (!is.na(p_var) && p_var > 0.05) {
      test <- t.test(x, y, var.equal = TRUE)
    } else {
      test <- t.test(x, y, var.equal = FALSE)  # Welch
    }
  } else {
    test <- wilcox.test(x, y, exact = FALSE)
  }
  list(p.value = test$p.value, method = test$method)
}
gene <- "FLT1"                           
celltype_of_interest <- "Stromal cell"    
group_var <- "group"                      
# -----------------------
obj_ct <- subset(scobj, subset = Cell_Type == celltype_of_interest)
obj_ct[[group_var]][,1] <- factor(obj_ct[[group_var]][,1], levels = c("Low","High"))
my_cols <- c("Low"="#4DA3FF", "High"="#E9892D")   # 自定义两组颜色（Low, High 顺序）
df <- FetchData(obj_ct, vars = c(gene, group_var)) %>%
  rename(expr = !!gene, grp = !!group_var) %>%
  filter(!is.na(expr), !is.na(grp))
df$grp <- factor(df$grp, levels = c("Low","High"))
stats_tbl <- df %>%
  group_by(grp) %>%
  summarise(n = n(),
            mean = mean(expr, na.rm=TRUE),
            median = median(expr, na.rm=TRUE),
            .groups="drop")
print(stats_tbl)   
res <- get_p_value_expr(df, value_col = "expr", group_col = "grp",
                        g1 = "Low", g2 = "High")
pval <- res$p.value
message(sprintf("P-value = %.4g | method = %s", pval, res$method))
p_label <- ifelse(is.na(pval), "NA",
           ifelse(pval < 1e-4, "****",
           ifelse(pval < 1e-3, "***",
           ifelse(pval < 1e-2, "**",
           ifelse(pval < 0.05, "*", "ns")))))
y_top   <- max(df$expr, na.rm = TRUE)
y_label <- y_top * 1.12
y_lim   <- y_top * 1.22
p <- VlnPlot(
  obj_ct, features = gene, group.by = group_var,
  pt.size = 0,                     # 关键：不画散点
  assay = "RNA", slot = "data"
) +
  theme_classic(base_size = 12) +
  labs(title = paste0("Violin | ", gene),
       x = NULL, y = "log-normalized expression") +
  coord_cartesian(ylim = c(NA, y_lim), clip = "off") +
  scale_fill_manual(values = my_cols, breaks = c("Low","High")) +
  scale_color_manual(values = my_cols, breaks = c("Low","High")) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = NA,
               color = "black", linewidth = 0.4)
anno <- data.frame(
  group1 = "Low", group2 = "High",
  y.position = y_label,
  label = p_label
)
p + ggpubr::stat_pvalue_manual(
      anno, label = "label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.02,
      inherit.aes = FALSE
    )

###Figure 8I
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
plot_km_survival <- function(df,
                             time_col = "OS.time",
                             event_col = "OS",
                             score_col = "Score",
                             cohort_name = "Cohort",
                             output_file = "KM_survival.pdf",
                             cut_method = c("mean", "median"),
                             palette = c("grey", "red"),
                             width = 5,
                             height = 5) {
  cut_method <- match.arg(cut_method)
  # --- Step 1: Dichotomize patients
  cut_value <- if (cut_method == "mean") {
    mean(df[[score_col]], na.rm = TRUE)
  } else {
    median(df[[score_col]], na.rm = TRUE)
  }
  df <- df %>%
    mutate(Group = ifelse(.data[[score_col]] >= cut_value, "High", "Low"))
  Group <- factor(df$Group, levels = c("Low", "High"))
  # --- Step 2: Survival object
  surv_obj <- Surv(time = df[[time_col]], event = df[[event_col]])

  # --- Step 3: Log-rank test
  log_rank_test <- survdiff(surv_obj ~ Group, data = df)
  p.val <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
  p_text <- ifelse(p.val < 0.001, "p < 0.001", paste0("p = ", round(p.val, 3)))
  # --- Step 4: Cox regression
  cox_model <- coxph(surv_obj ~ Group, data = df)
  summary_cox <- summary(cox_model)
  HR <- round(summary_cox$coefficients[2], 2)
  low95 <- round(summary_cox$conf.int[3], 2)
  up95 <- round(summary_cox$conf.int[4], 2)
  HR_text <- paste0("Hazard Ratio = ", HR)
  CI_text <- paste0("95% CI: ", low95, "-", up95)
  # --- Step 5: Group sizes
  n_low <- sum(Group == "Low")
  n_high <- sum(Group == "High")
  legend_labels <- c(paste0("Low (N=", n_low, ")"),
                     paste0("High (N=", n_high, ")"))
  # --- Step 6: KM plot
  p <- ggsurvplot(
    survfit(surv_obj ~ Group), 
    data = df,
    conf.int = FALSE,
    censor = FALSE,
    palette = palette,
    legend.title = "Group",
    legend.labs = legend_labels,
    font.legend = 11,
    pval = paste(p_text, HR_text, CI_text, sep = "\n"),
    ggtheme = theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        panel.grid = element_blank(),
        legend.position = c(0.25, 0.25),
        legend.justification = c("left", "bottom"),
        legend.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
      ),
    linetype = "solid",
    size = 1.5
  )
  # --- Step 7: Save PDF
  pdf(output_file, width = width, height = height)
  print(p$plot + ggtitle(cohort_name))
  dev.off()
  message("✅ Kaplan–Meier survival curve saved to: ", output_file)
  return(list(plot = p, cox = summary_cox, pval = p.val))
}
results <- plot_km_survival(
  df = res.cat,
  time_col = "OS.time",
  event_col = "OS",
  score_col = "Score",
  cohort_name = "KICH=KIRC",
  output_file = "Figure 8I.pdf",
  cut_method = "mean"   # or "median"
)

###Analysis in Figures 8J and 8K was performed via the HPA database







