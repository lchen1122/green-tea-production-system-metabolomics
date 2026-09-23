# ============================================================
# QC evaluation for untargeted metabolomics
#
# A: cumulative distribution of QC feature RSD
# B: relative total QC signal across analytical sequence
# C: PCA of analytical sequence with QC samples highlighted
# ============================================================

rm(list = ls())
gc()

# ============================================================
# 0. Packages
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# 1. File paths
# ============================================================

base_dir <- "F:/科研/有机生态低碳茶叶分析/代谢组结果/1.Preprocess/org_mix"

abund_file <- file.path(
  base_dir,
  "metab_abund.txt"
)

rsd_file <- file.path(
  base_dir,
  "metab_abund.txt_rsd.txt"
)

order_file <- file.path(
  base_dir,
  "order.txt"
)

out_dir <- file.path(
  base_dir,
  "QC_review_output"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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
# 3. Prepare injection order
# ============================================================

colnames(order_raw)[1:2] <- c(
  "Sample",
  "InjectionOrder"
)

order_df <- order_raw %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    InjectionOrder = as.numeric(InjectionOrder)
  ) %>%
  dplyr::arrange(
    InjectionOrder
  )

print(order_df)

# ============================================================
# 4. Identify sample columns
# ============================================================

sample_cols <- intersect(
  order_df$Sample,
  colnames(abund_raw)
)

if (length(sample_cols) == 0) {
  stop("No sample columns matched between metab_abund.txt and order.txt.")
}

order_df <- order_df %>%
  dplyr::filter(
    Sample %in% sample_cols
  ) %>%
  dplyr::arrange(
    InjectionOrder
  )

sample_cols <- order_df$Sample

qc_cols <- sample_cols[
  grepl("^QC", sample_cols)
]

study_cols <- sample_cols[
  !grepl("^QC", sample_cols)
]

cat(
  "\nMatched samples:",
  length(sample_cols),
  "\n"
)

cat(
  "QC samples:",
  paste(qc_cols, collapse = ", "),
  "\n"
)

# ============================================================
# 5. Abundance matrix
# ============================================================

abund_df <- abund_raw[
  ,
  sample_cols,
  drop = FALSE
]

abund_df[] <- lapply(
  abund_df,
  function(x) {
    suppressWarnings(
      as.numeric(x)
    )
  }
)

abund_mat <- as.matrix(
  abund_df
)

if ("metab_id" %in% colnames(abund_raw)) {
  rownames(abund_mat) <- abund_raw$metab_id
}

# ============================================================
# 6. PANEL A
# QC feature RSD cumulative distribution
#
# metab_abund.txt_rsd.txt:
# V1 = feature ID
# V2 = RSD expressed as a proportion
# e.g. 0.30 = 30%
# ============================================================

rsd_df <- data.frame(
  Feature = rsd_raw$V1,
  RSD = as.numeric(rsd_raw$V2) * 100,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    !is.na(RSD),
    is.finite(RSD),
    RSD >= 0
  )

rsd_curve <- rsd_df %>%
  dplyr::arrange(
    RSD
  ) %>%
  dplyr::mutate(
    CumPercent =
      dplyr::row_number() /
      dplyr::n() * 100
  )

pct_rsd_15 <- mean(
  rsd_df$RSD <= 15,
  na.rm = TRUE
) * 100

pct_rsd_20 <- mean(
  rsd_df$RSD <= 20,
  na.rm = TRUE
) * 100

pct_rsd_30 <- mean(
  rsd_df$RSD <= 30,
  na.rm = TRUE
) * 100

cat(
  "\nRSD <=15%:",
  round(pct_rsd_15, 1),
  "%\n"
)

cat(
  "RSD <=20%:",
  round(pct_rsd_20, 1),
  "%\n"
)

cat(
  "RSD <=30%:",
  round(pct_rsd_30, 1),
  "%\n"
)

pA <- ggplot(
  rsd_curve,
  aes(
    x = RSD,
    y = CumPercent
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 32,
    y = pct_rsd_30,
    label = paste0(
      "RSD ≤30%: ",
      round(pct_rsd_30, 1),
      "%"
    ),
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
  theme_classic(
    base_size = 12
  )

pA

# ============================================================
# 7. PANEL B
# Relative total QC signal across analytical sequence
#
# Each QC total signal is expressed relative to the
# mean of the three QC injections.
# ============================================================

signal_df <- data.frame(
  Sample = sample_cols,
  TotalSignal = colSums(
    abund_mat[
      ,
      sample_cols,
      drop = FALSE
    ],
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    order_df,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Type = ifelse(
      grepl("^QC", Sample),
      "QC",
      "Study"
    )
  ) %>%
  dplyr::arrange(
    InjectionOrder
  )

signal_qc_df <- signal_df %>%
  dplyr::filter(
    Type == "QC"
  ) %>%
  dplyr::mutate(
    RelativeSignal =
      TotalSignal /
      mean(
        TotalSignal,
        na.rm = TRUE
      ) * 100
  )

print(signal_qc_df)

pB <- ggplot(
  signal_qc_df,
  aes(
    x = InjectionOrder,
    y = RelativeSignal
  )
) +
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  ggrepel::geom_text_repel(
    aes(
      label = Sample
    ),
    size = 3.2,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    breaks = signal_qc_df$InjectionOrder
  ) +
  scale_y_continuous(
    limits = c(75,125),
    breaks = seq(
      75,
      125,
      10
    )
  ) +
  labs(
    x = "Injection order",
    y = "Relative total signal (%)"
  ) +
  theme_classic(
    base_size = 12
  )

pB

# ============================================================
# 8. PANEL C
# PCA of analytical sequence
#
# Starting from the untransformed abundance matrix:
# 1. sum normalization
# 2. log10 transformation
# 3. PCA
# ============================================================

pca_mat <- abund_mat[
  ,
  sample_cols,
  drop = FALSE
]

pca_mat[
  !is.finite(pca_mat)
] <- 0

# ------------------------------------------------------------
# 8.1 Sum normalization
# ------------------------------------------------------------

sample_sums <- colSums(
  pca_mat,
  na.rm = TRUE
)

sample_sums[
  sample_sums == 0
] <- NA

pca_norm <- sweep(
  pca_mat,
  2,
  sample_sums,
  FUN = "/"
)

# multiply by a constant only for numerical convenience
pca_norm <- pca_norm * 1e6

pca_norm[
  !is.finite(pca_norm)
] <- 0

# ------------------------------------------------------------
# 8.2 log10 transformation
# ------------------------------------------------------------

pca_log <- log10(
  pca_norm + 1
)

# ------------------------------------------------------------
# 8.3 Remove zero-variance features
# ------------------------------------------------------------

keep_feature <- apply(
  pca_log,
  1,
  function(x) {
    sd(
      x,
      na.rm = TRUE
    ) > 0
  }
)

pca_log <- pca_log[
  keep_feature,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 8.4 PCA
# ------------------------------------------------------------

pca_res <- prcomp(
  t(pca_log),
  center = TRUE,
  scale. = TRUE
)

pca_var <- (
  pca_res$sdev^2 /
    sum(pca_res$sdev^2)
) * 100

pca_df <- data.frame(
  Sample = rownames(
    pca_res$x
  ),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    order_df,
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Type = ifelse(
      grepl("^QC", Sample),
      "QC",
      "Study"
    ),
    QC_label = ifelse(
      Type == "QC",
      paste0(
        Sample,
        " (",
        InjectionOrder,
        ")"
      ),
      NA
    )
  )

# ============================================================
# 9. Plot PCA
# ============================================================

pC <- ggplot() +
  
  # study samples faded
  geom_point(
    data = pca_df %>%
      dplyr::filter(
        Type == "Study"
      ),
    aes(
      x = PC1,
      y = PC2
    ),
    color = "grey75",
    alpha = 0.30,
    size = 2.2
  ) +
  
  # QC highlighted
  geom_point(
    data = pca_df %>%
      dplyr::filter(
        Type == "QC"
      ),
    aes(
      x = PC1,
      y = PC2
    ),
    shape = 17,
    size = 3.5
  ) +
  
  # QC labels
  ggrepel::geom_text_repel(
    data = pca_df %>%
      dplyr::filter(
        Type == "QC"
      ),
    aes(
      x = PC1,
      y = PC2,
      label = QC_label
    ),
    size = 3.2,
    max.overlaps = Inf
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(
        pca_var[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        pca_var[2],
        1
      ),
      "%)"
    )
  ) +
  
  theme_classic(
    base_size = 12
  )

pC

# ============================================================
# 10. Combine A + B + C
# ============================================================

combined_plot <- pA / (pB + pC) +
  patchwork::plot_layout(
    heights = c(1, 1.05)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
  )
combined_plot

# ============================================================
# 11. Save figures
# ============================================================

ggsave(
  file.path(
    out_dir,
    "FigureA_QC_RSD_cumulative.pdf"
  ),
  pA,
  width = 5.5,
  height = 4.8
)

ggsave(
  file.path(
    out_dir,
    "FigureB_QC_relative_signal.pdf"
  ),
  pB,
  width = 5.5,
  height = 4.8
)

ggsave(
  file.path(
    out_dir,
    "FigureC_QC_PCA.pdf"
  ),
  pC,
  width = 5.5,
  height = 4.8
)

ggsave(
  file.path(
    out_dir,
    "Figure_ABC_QC_evaluation_combined.pdf"
  ),
  combined_plot,
  width = 15,
  height = 4.8
)

ggsave(
  file.path(
    out_dir,
    "Figure_ABC_QC_evaluation_combined.png"
  ),
  combined_plot,
  width = 15,
  height = 4.8,
  dpi = 600
)

# ============================================================
# 12. Export source data
# ============================================================

fwrite(
  rsd_df,
  file.path(
    out_dir,
    "PanelA_RSD_source_data.csv"
  )
)

fwrite(
  signal_qc_df,
  file.path(
    out_dir,
    "PanelB_QC_relative_signal.csv"
  )
)

fwrite(
  pca_df,
  file.path(
    out_dir,
    "PanelC_PCA_scores.csv"
  )
)

# ============================================================
# 13. Summary
# ============================================================

qc_summary <- data.frame(
  Metric = c(
    "No. of QC samples",
    "No. of QC-evaluated features",
    "Features with RSD <=15%",
    "Features with RSD <=20%",
    "Features with RSD <=30%"
  ),
  Value = c(
    length(qc_cols),
    nrow(rsd_df),
    round(pct_rsd_15, 1),
    round(pct_rsd_20, 1),
    round(pct_rsd_30, 1)
  )
)

print(qc_summary)

fwrite(
  qc_summary,
  file.path(
    out_dir,
    "QC_summary_statistics.csv"
  )
)

cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")