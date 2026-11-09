# ============================================================
# 02_metabolomics_QC_evaluation.R
# QC evaluation for untargeted metabolomics
#
# This script evaluates analytical stability only.
# No QC-based intensity-drift correction is performed.
# Samples were analyzed in the recorded preassigned sequence;
# no randomization is assumed or introduced here.
# ============================================================

rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# 1. Paths
# ============================================================

project_dir <- "F:/科研/有机生态低碳茶叶分析"

base_dir <- file.path(
  project_dir,
  "代谢组结果/1.Preprocess/org_mix"
)

abund_file <- file.path(base_dir, "metab_abund.txt")
rsd_file   <- file.path(base_dir, "metab_abund.txt_rsd.txt")
order_file <- file.path(base_dir, "order.txt")

out_dir <- file.path(
  project_dir,
  "Paper/Revised2/02_QC"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read files
# ============================================================

abund_raw <- fread(
  abund_file,
  data.table = FALSE,
  check.names = FALSE
)

rsd_raw <- fread(
  rsd_file,
  data.table = FALSE,
  check.names = FALSE
)

order_raw <- fread(
  order_file,
  data.table = FALSE,
  check.names = FALSE
)

# ============================================================
# 3. Injection sequence
# ============================================================

colnames(order_raw)[1:2] <- c("Sample", "InjectionOrder")

order_df <- order_raw %>%
  transmute(
    Sample = as.character(Sample),
    InjectionOrder = as.numeric(InjectionOrder)
  ) %>%
  arrange(InjectionOrder)

sample_cols <- intersect(order_df$Sample, colnames(abund_raw))

if (length(sample_cols) == 0) {
  stop("No sample columns matched between metab_abund.txt and order.txt.")
}

order_df <- order_df %>%
  filter(Sample %in% sample_cols) %>%
  arrange(InjectionOrder)

sample_cols <- order_df$Sample
qc_cols <- sample_cols[grepl("^QC", sample_cols)]
study_cols <- setdiff(sample_cols, qc_cols)

cat("\nMatched samples:", length(sample_cols), "\n")
cat("QC samples:", paste(qc_cols, collapse = ", "), "\n")
cat("Study samples:", length(study_cols), "\n")

fwrite(
  order_df,
  file.path(out_dir, "analytical_injection_sequence.csv")
)

# ============================================================
# 4. Abundance matrix
# ============================================================

abund_df <- abund_raw[, sample_cols, drop = FALSE]

abund_df[] <- lapply(
  abund_df,
  function(x) suppressWarnings(as.numeric(x))
)

abund_mat <- as.matrix(abund_df)

if ("metab_id" %in% colnames(abund_raw)) {
  rownames(abund_mat) <- abund_raw$metab_id
}

# ============================================================
# 5. Panel A: QC feature RSD cumulative distribution
# ============================================================

rsd_df <- data.frame(
  Feature = rsd_raw$V1,
  RSD = as.numeric(rsd_raw$V2) * 100,
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(RSD),
    is.finite(RSD),
    RSD >= 0
  )

rsd_curve <- rsd_df %>%
  arrange(RSD) %>%
  mutate(
    CumPercent = row_number() / n() * 100
  )

pct_rsd_15 <- mean(rsd_df$RSD <= 15, na.rm = TRUE) * 100
pct_rsd_20 <- mean(rsd_df$RSD <= 20, na.rm = TRUE) * 100
pct_rsd_30 <- mean(rsd_df$RSD <= 30, na.rm = TRUE) * 100

pA <- ggplot(
  rsd_curve,
  aes(x = RSD, y = CumPercent)
) +
  geom_line(linewidth = 0.8) +
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 32,
    y = pct_rsd_30,
    label = paste0("RSD <=30%: ", round(pct_rsd_30, 1), "%"),
    hjust = 0,
    vjust = -0.4,
    size = 3.5
  ) +
  coord_cartesian(
    xlim = c(0, 100),
    ylim = c(0, 100)
  ) +
  labs(
    x = "QC feature RSD (%)",
    y = "Cumulative percentage of features (%)"
  ) +
  theme_classic(base_size = 12)

# ============================================================
# 6. Panel B: relative total QC signal across analytical sequence
#    Descriptive diagnostic only; no drift correction is applied.
# ============================================================

signal_df <- data.frame(
  Sample = sample_cols,
  TotalSignal = colSums(
    abund_mat[, sample_cols, drop = FALSE],
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  left_join(order_df, by = "Sample") %>%
  mutate(
    Type = ifelse(grepl("^QC", Sample), "QC", "Study")
  ) %>%
  arrange(InjectionOrder)

signal_qc_df <- signal_df %>%
  filter(Type == "QC") %>%
  mutate(
    RelativeSignal =
      TotalSignal / mean(TotalSignal, na.rm = TRUE) * 100
  )

pB <- ggplot(
  signal_qc_df,
  aes(x = InjectionOrder, y = RelativeSignal)
) +
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    aes(label = Sample),
    size = 3.2,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    breaks = signal_qc_df$InjectionOrder
  ) +
  labs(
    x = "Injection order",
    y = "Relative total signal (%)"
  ) +
  theme_classic(base_size = 12)

# ============================================================
# 7. Panel C: PCA with QC samples highlighted
#
# This PCA is a QC diagnostic generated from metab_abund.txt.
# The transformation here is used only for QC visualization
# and is independent of the log2 analytical matrix used in
# downstream differential analysis.
# ============================================================

pca_mat <- abund_mat[, sample_cols, drop = FALSE]
pca_mat[!is.finite(pca_mat)] <- 0

sample_sums <- colSums(pca_mat, na.rm = TRUE)
sample_sums[sample_sums == 0] <- NA

pca_norm <- sweep(
  pca_mat,
  2,
  sample_sums,
  FUN = "/"
) * 1e6

pca_norm[!is.finite(pca_norm)] <- 0

pca_log <- log10(pca_norm + 1)

keep_feature <- apply(
  pca_log,
  1,
  function(x) sd(x, na.rm = TRUE) > 0
)

pca_log <- pca_log[keep_feature, , drop = FALSE]

pca_res <- prcomp(
  t(pca_log),
  center = TRUE,
  scale. = TRUE
)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(order_df, by = "Sample") %>%
  mutate(
    Type = ifelse(grepl("^QC", Sample), "QC", "Study"),
    QC_label = ifelse(
      Type == "QC",
      paste0(Sample, " (", InjectionOrder, ")"),
      NA_character_
    )
  )

pC <- ggplot() +
  geom_point(
    data = pca_df %>% filter(Type == "Study"),
    aes(x = PC1, y = PC2),
    color = "grey75",
    alpha = 0.30,
    size = 2.2
  ) +
  geom_point(
    data = pca_df %>% filter(Type == "QC"),
    aes(x = PC1, y = PC2),
    shape = 17,
    size = 3.5
  ) +
  ggrepel::geom_text_repel(
    data = pca_df %>% filter(Type == "QC"),
    aes(x = PC1, y = PC2, label = QC_label),
    size = 3.2,
    max.overlaps = Inf
  ) +
  labs(
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_classic(base_size = 12)

# ============================================================
# 8. Combine and export
# ============================================================

combined_plot <- pA / (pB + pC) +
  patchwork::plot_layout(heights = c(1, 1.05)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(out_dir, "Supplementary_Figure_QC_evaluation.pdf"),
  combined_plot,
  width = 15,
  height = 4.8
)

ggsave(
  file.path(out_dir, "Supplementary_Figure_QC_evaluation.png"),
  combined_plot,
  width = 15,
  height = 4.8,
  dpi = 600
)

fwrite(rsd_df, file.path(out_dir, "PanelA_RSD_source_data.csv"))
fwrite(signal_qc_df, file.path(out_dir, "PanelB_QC_relative_signal.csv"))
fwrite(pca_df, file.path(out_dir, "PanelC_PCA_scores.csv"))

qc_summary <- data.frame(
  Metric = c(
    "No. of QC samples",
    "No. of QC-evaluated features",
    "Features with RSD <=15% (%)",
    "Features with RSD <=20% (%)",
    "Features with RSD <=30% (%)",
    "QC-based intensity drift correction"
  ),
  Value = c(
    as.character(length(qc_cols)),
    as.character(nrow(rsd_df)),
    as.character(round(pct_rsd_15, 1)),
    as.character(round(pct_rsd_20, 1)),
    as.character(round(pct_rsd_30, 1)),
    "No"
  )
)

fwrite(
  qc_summary,
  file.path(out_dir, "QC_summary_statistics.csv")
)

cat("\nDone.\nOutput directory:\n", out_dir, "\n")
