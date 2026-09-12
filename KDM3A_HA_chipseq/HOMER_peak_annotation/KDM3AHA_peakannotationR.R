library(dplyr)
library(tidyr)
library(ggplot2)

# annotation distribution across the three KDM3A peak sets
# R. Siebeler

ANNOT_DIR <- "/path/to/bedtools/Peak_annotation"

# tidy up a HOMER annotation column and count occurrences
clean_and_count <- function(df) {
  df$Annotation <- gsub("\\(.*?\\)", "", df$Annotation)      # drop stuff in parens
  df$Annotation <- gsub("\\s*\\.\\d+$", "", df$Annotation)   # drop .1/.2 suffixes
  df$Annotation <- gsub("\\s+", " ", df$Annotation)
  df$Annotation <- trimws(df$Annotation)
  df <- df[!is.na(df$Annotation), ]
  table(df$Annotation)
}

df1 <- read.delim(file.path(ANNOT_DIR, "KDM3A_CTRL_specificAnnotated_peaks.txt"), header = TRUE, stringsAsFactors = FALSE)
df2 <- read.delim(file.path(ANNOT_DIR, "KDM3A_LPS_specificAnnotated_peaks.txt"), header = TRUE, stringsAsFactors = FALSE)
df3 <- read.delim(file.path(ANNOT_DIR, "intersect_KDM3A_CTRL_LPSAnnotated_peaks.txt"), header = TRUE, stringsAsFactors = FALSE)

df1_counts <- as.data.frame(clean_and_count(df1))
df2_counts <- as.data.frame(clean_and_count(df2))
df3_counts <- as.data.frame(clean_and_count(df3))

colnames(df1_counts) <- c("Annotation", "CTRL")
colnames(df2_counts) <- c("Annotation", "LPS")
colnames(df3_counts) <- c("Annotation", "Shared")

# merge the three sets by annotation type, 0 where a category is missing
annotation_counts <- full_join(df1_counts, df2_counts, by = "Annotation") %>%
  full_join(df3_counts, by = "Annotation") %>%
  replace(is.na(.), 0)

total_df1 <- sum(df1_counts$CTRL)
total_df2 <- sum(df2_counts$LPS)
total_df3 <- sum(df3_counts$Shared)

annotation_counts$CTRL_percent   <- (annotation_counts$CTRL / total_df1) * 100
annotation_counts$LPS_percent    <- (annotation_counts$LPS / total_df2) * 100
annotation_counts$Shared_percent <- (annotation_counts$Shared / total_df3) * 100

print(annotation_counts)

annotation_counts_long <- annotation_counts %>%
  pivot_longer(cols = ends_with("_percent"),
               names_to = "Dataframe",
               values_to = "Percentage") %>%
  mutate(Annotation = factor(Annotation, levels = unique(annotation_counts$Annotation)))

annotation_colors <- c(
  "Intergenic"   = "#1f77b4",
  "exon"         = "#ff7f0e",
  "intron"       = "#2ca02c",
  "promoter-TSS" = "#4f6d7a",
  "3' UTR"       = "#7f8c8d",
  "5' UTR"       = "#34495e",
  "TTS"          = "#f39c12",
  "non-coding"   = "#2d3e50"
)

ggplot(annotation_counts_long, aes(x = Dataframe, y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Annotation Percentages by Dataframe",
    x = "Dataframe",
    y = "Percentage (%)"
  ) +
  scale_fill_manual(values = annotation_colors) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white")
  )
