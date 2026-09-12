## ================================================================
## KDM3A expression across leukocyte subtypes
## ================================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# --- settings, edit before running ---
SEURAT_OBJECT_PATH <- "/path/to/scData/final.pop.call.integrated.full.seurat.Rds"

# --- Load data ---
final.pop.call.integrated.full.seurat <- readRDS(SEURAT_OBJECT_PATH)
DefaultAssay(final.pop.call.integrated.full.seurat) <- "RNA"

# --- 1. Map original identities to main leukocyte subtypes ---
current_idents <- as.character(Idents(final.pop.call.integrated.full.seurat))

subtype_map <- c(
  "CLEC9A+ cDC2"                                                                          = "Dendritic Cells",
  "CD1C+ cDC1"                                                                            = "Dendritic Cells",
  "CD79+ B Cells"                                                                         = "B Cells",
  "CD79+ Plasma B Cells"                                                                  = "B Cells",
  "CD79+ Class-switched Memory B Cells"                                                   = "B Cells",
  "KIT+ Mast Cells"                                                                       = "Mast Cells",
  "GYPA+ Erythroid Cells"                                                                 = "Erythroid Cells",
  "CD3+MKI67+ Proliferating T Cells"                                                      = "T Cells",
  "FOXP3+ T Cells"                                                                        = "T Cells",
  "CD8+ T Cells"                                                                          = "T Cells",
  "CD4+ T Cells"                                                                          = "T Cells",
  "ACTA2+ Smooth Muscle Cells"                                                            = "Stromal Cells",
  "CD34+ Endothelial Cells"                                                               = "Stromal Cells",
  "CD56+KLRC1+ NK Cells"                                                                  = "NK Cells",
  "CD56-CD16+ NK Cells"                                                                   = "NK Cells",
  "CD14+CD16+CD64+SLAN- Intermediate Monocytes"                                           = "Monocytes",
  "CD14+CD16+CD64-SLAN+ Non-Classical Monocytes"                                          = "Monocytes",
  "CD14+CD16-CD64+SLAN- Classical Monocytes"                                              = "Monocytes",
  "CD14+TREM2-OLR1+ABCA+ Foamy Macrophages"                                               = "Macrophages",
  "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages"         = "Macrophages",
  "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages"             = "Macrophages",
  "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages"    = "Macrophages",
  "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages"                                      = "Macrophages",
  "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages"                                 = "Macrophages",
  "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages"                      = "Macrophages",
  "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages"                                   = "Macrophages",
  "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages"            = "Macrophages"
)

missing_sub <- setdiff(unique(current_idents), names(subtype_map))
if (length(missing_sub) > 0) {
  message("WARNING - identities missing from subtype_map (will be NA): ")
  print(missing_sub)
}

subtype_order <- c(
  "Macrophages", "Monocytes", "Dendritic Cells", "Mast Cells",
  "T Cells", "NK Cells", "B Cells", "Erythroid Cells", "Stromal Cells"
)

final.pop.call.integrated.full.seurat$leukocyte_subtype <- factor(
  unname(subtype_map[current_idents]),
  levels = subtype_order
)

# Sanity check - confirm no empty/unmatched groups
print(table(final.pop.call.integrated.full.seurat$leukocyte_subtype, useNA = "ifany"))

# --- 2. Build the three panels ---

p0 <- DimPlot(
  final.pop.call.integrated.full.seurat,
  reduction = "umap",
  group.by = "leukocyte_subtype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle("Leukocyte Subtypes") +
  scale_color_manual(
    values = c(
      "Macrophages"       = "dodgerblue4",
      "Monocytes"         = "firebrick4",
      "Dendritic Cells"   = "#F781BF",
      "Mast Cells"        = "#984EA3",
      "T Cells"           = "#1B9E77",
      "NK Cells"          = "#00A6CA",
      "B Cells"           = "#FF7F00",
      "Erythroid Cells"   = "#E41A1C",
      "Stromal Cells"     = "#4DAF4A"
    ))

p1 <- FeaturePlot(
  final.pop.call.integrated.full.seurat,
  features = "KDM3A",
  reduction = "umap",
  cols = c("grey90", "dodgerblue4"),
  pt.size = 0.4,
  order = TRUE
) + ggtitle("KDM3A Expression")

p2_subtype <- DotPlot(
  final.pop.call.integrated.full.seurat,
  features = "KDM3A",
  group.by = "leukocyte_subtype",
  scale = FALSE,        # plot raw average expression, not z-scored
  dot.scale = 12,
  scale.min = 0,
  scale.max = 20
) +
  coord_flip() +
  ggtitle("KDM3A across main leukocyte subtypes") +
  scale_color_gradient(low = "grey90", high = "dodgerblue4", limits = c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- 3. Combined layout: p0 + p1 on top, p2_subtype spanning below ---

(p0 + p1) / p2_subtype +
  plot_layout(heights = c(1.2, 1))

# --- 4. Summary stats table (raw, unscaled - matches dot plot colors) ---

expr_subtype <- FetchData(final.pop.call.integrated.full.seurat, vars = "KDM3A")
expr_subtype$group <- final.pop.call.integrated.full.seurat$leukocyte_subtype

summary_stats_subtype <- expr_subtype %>%
  group_by(group) %>%
  summarise(
    mean_expr = mean(KDM3A),
    pct_expr  = mean(KDM3A > 0) * 100,
    n_cells   = n()
  ) %>%
  arrange(desc(mean_expr))

print(summary_stats_subtype)
