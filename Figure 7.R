###Figure 7A-C,E
library(CytoTRACE2)
library(Seurat)
library(ggplot2)
library(patchwork)
library(devtools)
result_sce <- cytotrace2(scobj,
                         is_seurat = TRUE,
                         slot_type = "counts",
                         species = "human",
                         seed = 1234)
result_sce<-readRDS("result_sce.rds")
annotation <- data.frame(phenotype = result_sce@meta.data$cell_groups) %>% set_rownames(., colnames(result_sce))
plots <- plotData(cytotrace2_result = result_sce, annotation = annotation, is_seurat = TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)

###Figure 7D
library(ggplot2)
plot_correlation_scatter <- function(df,
                                     x_col = "PBIS.TP53",
                                     y_col = "CytoTRACE2_Score",
                                     output_file = "correlation_scatter.pdf",
                                     sample_size = 5000,
                                     seed = 123,
                                     width = 5,
                                     height = 4) {
  set.seed(seed)
  # --- Step 1: sample for normality test
  sample_df <- df[sample(1:nrow(df), min(sample_size, nrow(df))), ]
  # --- Step 2: Shapiro-Wilk test
  shapiro_x <- shapiro.test(sample_df[[x_col]])
  shapiro_y <- shapiro.test(sample_df[[y_col]])
  # --- Step 3: choose correlation method
  method_used <- if (shapiro_x$p.value > 0.05 & shapiro_y$p.value > 0.05) {
    "pearson"
  } else {
    "spearman"
  }
  # --- Step 4: correlation test on full data
  cor_test <- cor.test(df[[x_col]], df[[y_col]], method = method_used)
  r_value <- round(cor_test$estimate, 3)
  p_value <- signif(cor_test$p.value, 3)
  cor_label <- paste0(toupper(method_used), " R = ", r_value, "\nP = ", p_value)
  # --- Step 5: plot
  pdf(file = output_file, width = width, height = height)
  p <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(color = "black", alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid", size = 1) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.2,
             label = cor_label, size = 5) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank()
    ) +
    labs(x = x_col, y = y_col)
  print(p)
  dev.off()
  message("✅ Scatter plot with correlation saved to: ", output_file)
  return(list(cor_test = cor_test, plot = p))
}
results <- plot_correlation_scatter(
  df,
  x_col = "PBIS.TP53",
  y_col = "CytoTRACE2_Score",
  output_file = "Figure 7D.pdf"
)

###Figure 7F
library(ggsci)
library(cowplot)
library(dplyr)
library(data.table)
library(tidyverse)
library(paletteer) 
library(vcd)
library(limma)
library(gplots)
library(Seurat)
library(clustree)
library(ggplot2)  
library(reshape2)  
library(dplyr)  
library(monocle)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(gridExtra)  
library(ggnewscale)
library(org.Mm.eg.db)  
library(ComplexHeatmap) 
library(ClusterGVis)
scRNAsub<-readRDS("result_sce.rds")
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')  
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,  
                          phenoData = pd,  
                          featureData = fd,  
                          expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 16, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
mycds <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
Idents(scRNAsub) <- scRNAsub$cell_groups
diff.wilcox <- FindAllMarkers(scRNAsub)
all.markers <- diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val < 0.05)
diff.genes <- subset(all.markers, p_val_adj < 0.05)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
mycds <- reduceDimension(mycds,
                           residualModelFormulaStr = "~orig.ident", #去除样本影响
                           reduction_method = 'DDRTree')
mycds <- orderCells(mycds)
save(mycds,file = "mycds.monocle2.Rdata")
load("mycds.monocle2.Rdata")
GM_state <- function(cds,x){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(mycds)$State, pData(mycds)$RNA_snn_res.0.4)[,x]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
mycds <- orderCells(mycds, root_state = GM_state(mycds,"10"))
write.table(pData(mycds),file="time.txt",sep="\t",row.names=F,quote=F)
 data_df <- t(reducedDimS(mycds)) %>%
  as.data.frame() %>%
    select_('Component 1' = 1, 'Component 2' = 2) %>%
    rownames_to_column("Cells") %>%
    mutate(pData(mycds)$State,
           pData(mycds)$Pseudotime,
           pData(mycds)$orig.ident,
           pData(mycds)$cell_groups)
colnames(data_df) <- c("cells","Component_1","Component_2","State",
                         "Pseudotime","orig.ident","cell_groups")
  reduced_dim_coords <- reducedDimK(mycds)
  ca_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  p1=ggplot() +
    geom_point_rast(data = data_df, aes(x = Component_1,
                                        y = Component_2,
                                        color =Pseudotime),size =0.5) +
    scale_color_viridis()+
    theme_bw()+
    #theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
    theme(
      axis.line = element_line(color = "black", size = 0.5),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.length = unit(0.1, "lines"),
      axis.ticks = element_blank(),
      axis.title = element_text(size=15),
    )

###Figure 7G
p2=ggplot() +
    geom_point_rast(data = data_df, aes(x = Component_1,
                                        y = Component_2, color = cell_groups),size =0.5) +
    scale_color_npg()+
    theme_bw()+
    #theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.length = unit(0.1, "lines"),
      axis.ticks = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.title = element_text(size=15),
    )+ facet_wrap(~cell_groups, nrow = 1)
  w=length(unique(data_df$cell_groups))*3

###Figure 7H
diff_test<-differentialGeneTest(mycds,cores=16,fullModelFormulaStr="~sm.ns(Pseudotime)")
gene_num=200
diff_genes<-diff_test%>%arrange(qval)%>%head(gene_num)%>%dplyr::select(gene_short_name)
num_clusters=4
df<-plot_pseudotime_heatmap2(mycds[diff_genes[,1],],
                               num_clusters=num_clusters,
                               cores=1,
                               use_gene_short_name=TRUE,
                               show_rownames=TRUE)
organism="hsa"
pvalue_cutoff=0.05
topn=5
seed=123456
enrich_db="org.Hs.eg.db"
enrich<-enrichCluster(object=df,
                        OrgDb=enrich_db,
                        type="BP",
                        organism=organism,
                        pvalueCutoff=pvalue_cutoff,
                        topn=topn,
                        seed=seed)
markGenes<-sample(unique(df$wide.res$gene),25,replace=FALSE)
visCluster(object=df,
           plot.type="both",
           column_names_rot=45,
           show_row_dend=FALSE,
           markGenes=sample(df$wide.res$gene,gene_num,replace=FALSE),
           markGenes.side="right",
           annoTerm.data=enrich,
           go.col=rep(jjAnno::useMyCol("calm",n=num_clusters),each=5),
           line.side="left")


