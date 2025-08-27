###Figure 5B
library(Seurat)
library(devtools)
library(dplyr)
library(scRNAtoolVis)
library(SCP)
library(ggplot2)
library(forcats)
library(patchwork)
CellDimPlot(scobj,group.by="seurat_clusters",reduction="umap")

###Figure 5C
CellDimPlot(scobj,group.by="Cell_type",reduction="umap")

###Figure 5D
markerGene <- data.frame(
  gene = c("CD3D","CD2","CD3E","NKG7","CCL5","GNLY","PRF1","LYZ","RGS2","COL1A1","VWF","TOP2A",
"CDH4","SLC26A7","CLCNKB","NOX4","CLDN1","KRT8"),
  cluster = c("Lymphoid cell","Lymphoid cell","Lymphoid cell","Lymphoid cell","Lymphoid cell","Lymphoid cell","Lymphoid cell","Myeloid cell","Myeloid cell","Stromal cell","Stromal cell",
"Cycling cell","Epithelial cell","Epithelial cell","Epithelial cell","Epithelial cell","Epithelial cell","Epithelial cell")
)
jjDotPlot(object = scobj,
          markerGene = markerGene,
          anno = TRUE,
          plot.margin = c(3, 1, 1, 1))###Figure 5D

###Figure 5E
cal1_cols <- c(
  "Cycling cell"    = "#E0B0FF",
  "Epithelial cell" = "#A7C7E7",
  "Lymphoid cell"   = "#AFE1AF",
  "Myeloid cell"    = "#BDB5D5",
  "Stromal cell"    = "#FFB6C1"
)
tumor_type_map <- c(
  "T14_5" = "KIRC",
  "T52_6" = "KIRC",
  "T54_5" = "KIRC",
  "T7_8"  = "KICH",
  "T56_2" = "KICH",
  "T36_3" = "KICH"
)
scobj@meta.data$TumorType <- tumor_type_map[scobj@meta.data$orig.ident]
p2 <- scobj@meta.data %>%
  ggplot(aes(y = fct_rev(fct_infreq(Cell_Type)), fill = Cell_Type)) +
  geom_bar(stat = "count") +
  labs(x = "Cell count", y = NULL) +
  scale_fill_manual(name = "Cell Type", values = cal1_cols) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(size = 14, color = 'black'))
cell_counts <- as.data.frame(table(scobj@meta.data$TumorType, scobj@meta.data$Cell_Type))
colnames(cell_counts) <- c("TumorType", "CellType", "Count")
p3 <- ggplot(data = cell_counts, aes(x = TumorType, y = Count, fill = CellType)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = NULL, y = "Cell Type Frequency") +
  scale_fill_manual(name = "Cell Type", values = cal1_cols) +
  theme_bw(base_size = 16) +
  theme(axis.text = element_text(size = 14, color = 'black'))
p4 <- p2 + p3 + plot_layout(guides = "collect") &
  plot_annotation(
    title = "Atlas Composition by Tumor Type",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'))
  )

###Figure 5F
PBIS_TP53=read.table("PBIS-TP53.txt",header=T,sep="\t",check.names = F,row.names=1)
is_identical <- identical(colnames(scobj), rownames(PBIS_TP53))
print(is_identical)  # TRUE
scobj<- AddMetaData(scobj, metadata = PBIS_TP53)
my5color <- c(
  "Cycling cell"    = '#ee6666',
  "Stromal cell"    = '#fac858',
  "Epithelial cell" = '#5470c6',
  "Lymphoid cell"   = '#91cc75',
  "Myeloid cell"    = '#73c0de'
)
my_comparisons <- list(
  c("Cycling cell", "Stromal cell"),
  c("Cycling cell", "Epithelial cell"),
  c("Cycling cell", "Myeloid cell"),
  c("Cycling cell", "Lymphoid cell")
)
scobj@meta.data$Cell_Type <- factor(
  scobj@meta.data$Cell_Type,
  levels = c("Cycling cell", "Stromal cell", "Epithelial cell", "Myeloid cell", "Lymphoid cell")
)
library(ggpubr)
p.AddModuleScore <- ggviolin(
  data = scobj@meta.data,
  x = "Cell_Type",
  y = "PBIS.TP53",
  color = "Cell_Type",
  fill = "Cell_Type",
  add = "mean_sd",
  add.params = list(color = "black")
) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +  # 或 label = "p.format"
  scale_color_manual(values = my5color) +
  scale_fill_manual(values = my5color) +
  theme(
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = '', y = 'PBIS-TP53 Score')

###Figure 5G
CellDimPlot(scobj,group.by="CopyKAT",reduction="umap")

###Figure 5H
plot_copykat_stacked <- function(metadata_file,
                                 cell_type_col = "Cell_Type",
                                 group_col = "group",
                                 copykat_col = "copykat.pred",
                                 output_file = "copykat_vs_group_stacked.pdf",
                                 width = 5,
                                 height = 4) {
  # --- Step 1: load metadata
  df <- read.table(metadata_file, header = TRUE, sep = "\t")
  # --- Step 2: create "type" column (Cell_Type + group)
  df$type <- paste(df[[cell_type_col]], df[[group_col]])
  # --- Step 3: rename copykat column to "cell_groups"
  colnames(df)[colnames(df) == copykat_col] <- "cell_groups"
  # --- Step 4: filter out "not.defined"
  df <- df[df$cell_groups != "not.defined", ]
  # --- Step 5: summarize counts
  mydata <- df %>%
    group_by(type, cell_groups) %>%
    summarise(values = n(), .groups = "drop")
  # --- Step 6: calculate % within each type
  filtered_data <- mydata %>%
    group_by(type) %>%
    mutate(percent = values / sum(values) * 100)
  # --- Step 7: open PDF
  pdf(file = output_file, width = width, height = height, bg = "white")
  # --- Step 8: plot stacked barplot
  p <- ggplot(filtered_data,
              aes(x = type,
                  y = percent,
                  fill = factor(cell_groups,
                                levels = c("aneuploid", "diploid")))) +
    geom_col(width = 0.55, position = position_stack(reverse = TRUE), show.legend = TRUE) +
    geom_text(aes(label = paste0(round(percent, 1), "%")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 3) +
    scale_fill_manual(values = c("firebrick1", "dodgerblue3", "limegreen", "blueviolet", "orange1")) +
    theme_classic(base_size = 14) +
    theme(axis.line = element_line(size = 0.3),
          axis.ticks = element_line(size = 0.3),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 0.97, vjust = 0.95)) +
    labs(x = "",
         y = expression("Total of Tumor cells (%)"),
         fill = "Cluster")
  print(p)
  dev.off()
  message("Stacked barplot saved to: ", output_file)
  return(p)
}
plot_copykat_stacked(
  metadata_file = "metadata.txt",
  cell_type_col = "Cell_Type",
  group_col = "group",
  copykat_col = "copykat.pred",
  output_file = "Figure 5H.pdf"
)

###Figure 5I
CellDimPlot(scobj,group.by="SCEVAN",reduction="umap")

###Figure 5J
plot_copykat_stacked <- function(metadata_file,
                                 cell_type_col = "Cell_Type",
                                 group_col = "group",
                                 SCEVAN_col = "SCEVAN.pred",
                                 output_file = "SCEVAN_vs_group_stacked.pdf",
                                 width = 5,
                                 height = 4) {
  # --- Step 1: load metadata
  df <- read.table(metadata_file, header = TRUE, sep = "\t")
  # --- Step 2: create "type" column (Cell_Type + group)
  df$type <- paste(df[[cell_type_col]], df[[group_col]])
  # --- Step 3: rename SCEVAN column to "cell_groups"
  colnames(df)[colnames(df) == SCEVAN_col] <- "cell_groups"
  # --- Step 4: filter out "not.defined"
  df <- df[df$cell_groups != "Filtered", ]
  # --- Step 5: summarize counts
  mydata <- df %>%
    group_by(type, cell_groups) %>%
    summarise(values = n(), .groups = "drop")
  # --- Step 6: calculate % within each type
  filtered_data <- mydata %>%
    group_by(type) %>%
    mutate(percent = values / sum(values) * 100)
  # --- Step 7: open PDF
  pdf(file = output_file, width = width, height = height, bg = "white")
  # --- Step 8: plot stacked barplot
  p <- ggplot(filtered_data,
              aes(x = type,
                  y = percent,
                  fill = factor(cell_groups,
                                levels = c("tumor", "normal")))) +
    geom_col(width = 0.55, position = position_stack(reverse = TRUE), show.legend = TRUE) +
    geom_text(aes(label = paste0(round(percent, 1), "%")),
              position = position_stack(vjust = 0.5, reverse = TRUE),
              size = 3) +
    scale_fill_manual(values = c("firebrick1", "dodgerblue3", "limegreen", "blueviolet", "orange1")) +
    theme_classic(base_size = 14) +
    theme(axis.line = element_line(size = 0.3),
          axis.ticks = element_line(size = 0.3),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 0.97, vjust = 0.95)) +
    labs(x = "",
         y = expression("Total of Tumor cells (%)"),
         fill = "Cluster")
  print(p)
  dev.off()
  message("Stacked barplot saved to: ", output_file)
  return(p)
}
plot_copykat_stacked(
  metadata_file = "metadata.txt",
  cell_type_col = "Cell_Type",
  group_col = "group",
  copykat_col = "SCEVAN.pred",
  output_file = "Figure 5J.pdf"
)
