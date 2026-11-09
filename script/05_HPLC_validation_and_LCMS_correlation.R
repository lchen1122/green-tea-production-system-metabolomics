# ============================================================
# 05_HPLC_validation_and_LCMS_correlation.R
# Targeted HPLC validation of EGCG, ECG, EGC, EC and caffeine
# and per-sample correlation with untargeted LC-MS signals
#
# LC-MS representative-feature selection follows the platform
# priority rule:
# 1) experimental Fragmentation Score > theoretical score
# 2) within the same evidence level, higher score
# 3) if tied, lower QC RSD
# 4) if still tied, higher overall abundance
# ============================================================

rm(list = ls())
gc()

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. Paths
# ============================================================

project_dir <- "F:/科研/有机生态低碳茶叶分析"

metab_file <- file.path(
  project_dir,
  "代谢组结果/1.Preprocess/mix/all_metab_data.txt"
)

hplc_file <- file.path(
  project_dir,
  "Paper/HPLC.xlsx"
)

out_dir <- file.path(
  project_dir,
  "Paper/Revised2/05_HPLC_validation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read untargeted metabolomics data
# ============================================================

df <- fread(
  metab_file,
  data.table = FALSE,
  check.names = FALSE
) %>%
  mutate(
    `Fragmentation Score` =
      suppressWarnings(as.numeric(`Fragmentation Score`)),
    `Theoretical Fragmentation Score` =
      suppressWarnings(as.numeric(`Theoretical Fragmentation Score`)),
    `Mass Error (ppm)` =
      suppressWarnings(as.numeric(`Mass Error (ppm)`))
  )

bio_cols <- grep(
  "^(CON|ELC|ORG)_",
  colnames(df),
  value = TRUE
)

qc_cols <- grep(
  "^QC",
  colnames(df),
  value = TRUE
)

# ============================================================
# 3. Define HPLC targets and annotation synonyms
# ============================================================

target_map <- data.frame(
  Metabolite = c(
    "Epigallocatechin-3-Gallate",
    "Epigallocatechin Gallate",
    "Epicatechin-3-Gallate",
    "(-)-Epicatechin 3-O-Gallate",
    "Epigallocatechin",
    "Epicatechin",
    "Caffeine"
  ),
  HPLC_target = c(
    "EGCG",
    "EGCG",
    "ECG",
    "ECG",
    "EGC",
    "EC",
    "Caffeine"
  ),
  stringsAsFactors = FALSE
)

target_df <- df %>%
  inner_join(target_map, by = "Metabolite")

if (nrow(target_df) == 0) {
  stop("No LC-MS features matched the HPLC target annotations.")
}

# ============================================================
# 4. Rank duplicate LC-MS annotations using platform rule
# ============================================================

calc_qc_rsd <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)

  if (!is.finite(m) || m == 0 || !is.finite(s)) {
    return(Inf)
  }

  s / m
}

if (length(qc_cols) > 1) {
  target_df$QC_RSD <- apply(
    target_df[, qc_cols, drop = FALSE],
    1,
    calc_qc_rsd
  )
} else {
  target_df$QC_RSD <- Inf
}

target_df$OverallMean <- rowMeans(
  as.data.frame(
    lapply(
      target_df[, bio_cols, drop = FALSE],
      as.numeric
    )
  ),
  na.rm = TRUE
)

target_ranked <- target_df %>%
  mutate(
    Evidence_rank = case_when(
      level == "B(i)" &
        !is.na(`Fragmentation Score`) &
        `Fragmentation Score` > 35 &
        abs(`Mass Error (ppm)`) < 10 ~ 1L,

      level == "B(ii)" &
        !is.na(`Theoretical Fragmentation Score`) &
        `Theoretical Fragmentation Score` > 40 &
        abs(`Mass Error (ppm)`) < 10 ~ 2L,

      TRUE ~ 3L
    ),

    Match_score = case_when(
      Evidence_rank == 1L ~ `Fragmentation Score`,
      Evidence_rank == 2L ~ `Theoretical Fragmentation Score`,
      TRUE ~ -Inf
    )
  ) %>%
  arrange(
    HPLC_target,
    Evidence_rank,
    desc(Match_score),
    QC_RSD,
    desc(OverallMean)
  )

representative_features <- target_ranked %>%
  group_by(HPLC_target) %>%
  slice_head(n = 1) %>%
  ungroup()

print(
  representative_features %>%
    select(
      HPLC_target,
      metab_id,
      Metabolite,
      Mode,
      level,
      `Mass Error (ppm)`,
      `Fragmentation Score`,
      `Theoretical Fragmentation Score`,
      QC_RSD,
      OverallMean
    )
)

write.csv(
  representative_features %>%
    select(
      HPLC_target,
      metab_id,
      Metabolite,
      `m/z`,
      `Retention time`,
      Mode,
      Adducts,
      Formula,
      level,
      `Mass Error (ppm)`,
      `Fragmentation Score`,
      `Theoretical Fragmentation Score`,
      QC_RSD,
      OverallMean
    ),
  file.path(out_dir, "LCMS_representative_features_for_HPLC_targets.csv"),
  row.names = FALSE
)

# ============================================================
# 5. LC-MS per-sample data for the five targets
# ============================================================

lcms_long <- representative_features %>%
  select(
    HPLC_target,
    metab_id,
    all_of(bio_cols)
  ) %>%
  pivot_longer(
    cols = all_of(bio_cols),
    names_to = "Sample",
    values_to = "LCMS_abundance"
  ) %>%
  mutate(
    LCMS_abundance = as.numeric(LCMS_abundance),
    Group = factor(
      sub("_\\d+$", "", Sample),
      levels = c("CON", "ELC", "ORG")
    )
  )

# ============================================================
# 6. Read HPLC data and average technical duplicates per sample
# ============================================================

hplc_raw <- read_excel(
  hplc_file
)

hplc_sample <- hplc_raw %>%
  transmute(
    Sample = as.character(Group),
    HPLC_target = as.character(Class),
    HPLC_mg_g = as.numeric(Value) * 10,
    Group = factor(
      sub(
        "_\\d+$",
        "",
        Sample
      ),
      levels = c(
        "CON",
        "ELC",
        "ORG"
      )
    )
  ) %>%
  filter(
    HPLC_target %in% c(
      "EGCG",
      "ECG",
      "EGC",
      "EC",
      "Caffeine"
    )
  )

nrow(hplc_sample)

table(
  hplc_sample$HPLC_target
)

hplc_sample %>%
  count(
    Sample,
    HPLC_target
  ) %>%
  filter(
    n != 1
  )



# ============================================================
# 7. Match HPLC and LC-MS per sample
# ============================================================

paired <- inner_join(
  hplc_sample,
  lcms_long,
  by = c("Sample", "HPLC_target", "Group")
)

write.csv(
  paired,
  file.path(out_dir, "HPLC_LCMS_per_sample_paired_data.csv"),
  row.names = FALSE
)

# ============================================================
# 8. Pearson and Spearman correlations
# ============================================================

safe_cor_test <- function(x, y, method) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(
      data.frame(
        n = length(x),
        estimate = NA_real_,
        p_value = NA_real_
      )
    )
  }

  ct <- cor.test(
    x,
    y,
    method = method,
    exact = FALSE
  )

  data.frame(
    n = length(x),
    estimate = unname(ct$estimate),
    p_value = ct$p.value
  )
}

correlation_results <- lapply(
  c("EGCG", "ECG", "EGC", "EC", "Caffeine"),
  function(target) {
    dat <- paired %>%
      filter(HPLC_target == target)

    pearson <- safe_cor_test(
      dat$HPLC_mg_g,
      dat$LCMS_abundance,
      "pearson"
    )

    spearman <- safe_cor_test(
      dat$HPLC_mg_g,
      dat$LCMS_abundance,
      "spearman"
    )

    bind_rows(
      cbind(
        HPLC_target = target,
        Method = "Pearson",
        pearson
      ),
      cbind(
        HPLC_target = target,
        Method = "Spearman",
        spearman
      )
    )
  }
) %>%
  bind_rows()

print(correlation_results)

write.csv(
  correlation_results,
  file.path(out_dir, "HPLC_LCMS_per_sample_correlations.csv"),
  row.names = FALSE
)

# ============================================================
# 9. Groupwise trend plots
# ============================================================

lcms_scaled <- lcms_long %>%
  group_by(HPLC_target) %>%
  mutate(
    Scaled_abundance = {
      lo <- min(LCMS_abundance, na.rm = TRUE)
      hi <- max(LCMS_abundance, na.rm = TRUE)

      if (hi == lo) {
        rep(0.5, n())
      } else {
        (LCMS_abundance - lo) / (hi - lo)
      }
    }
  ) %>%
  ungroup()

lcms_summary <- lcms_scaled %>%
  group_by(HPLC_target, Group) %>%
  summarise(
    Mean = mean(Scaled_abundance, na.rm = TRUE),
    SD = sd(Scaled_abundance, na.rm = TRUE),
    n = sum(!is.na(Scaled_abundance)),
    .groups = "drop"
  )

hplc_summary <- hplc_sample %>%
  group_by(HPLC_target, Group) %>%
  summarise(
    Mean = mean(HPLC_mg_g, na.rm = TRUE),
    SD = sd(HPLC_mg_g, na.rm = TRUE),
    n = sum(!is.na(HPLC_mg_g)),
    .groups = "drop"
  )

target_levels <- c("EGCG", "ECG", "EGC", "EC", "Caffeine")

lcms_summary$HPLC_target <- factor(
  lcms_summary$HPLC_target,
  levels = target_levels
)

hplc_summary$HPLC_target <- factor(
  hplc_summary$HPLC_target,
  levels = target_levels
)

p_lcms <- ggplot(
  lcms_summary,
  aes(x = Group, y = Mean)
) +
  geom_col(width = 0.68, color = "black") +
  geom_errorbar(
    aes(
      ymin = pmax(Mean - SD, 0),
      ymax = pmin(Mean + SD, 1)
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  facet_wrap(~ HPLC_target, nrow = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Scaled LC-MS abundance"
  ) +
  theme_classic(base_size = 12)

p_hplc <- ggplot(
  hplc_summary,
  aes(x = Group, y = Mean)
) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.12,
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ HPLC_target,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "HPLC content (mg/g)"
  ) +
  theme_classic(base_size = 12)

p_combined <- p_lcms / p_hplc +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(out_dir, "Supplementary_Figure_HPLC_LCMS_groupwise.pdf"),
  p_combined,
  width = 10,
  height = 7.2
)

# ============================================================
# 10. Per-sample correlation scatter plots
# ============================================================

cor_label <- correlation_results %>%
  filter(Method == "Pearson") %>%
  mutate(
    label = ifelse(
      is.na(estimate),
      "Pearson r = NA",
      paste0(
        "Pearson r = ",
        sprintf("%.2f", estimate),
        "\nP = ",
        format.pval(p_value, digits = 2, eps = 0.001)
      )
    )
  )

plot_dat <- paired %>%
  left_join(
    cor_label %>%
      select(HPLC_target, label),
    by = "HPLC_target"
  ) %>%
  mutate(
    HPLC_target = factor(
      HPLC_target,
      levels = target_levels
    )
  )

p_cor <- ggplot(
  plot_dat,
  aes(
    x = HPLC_mg_g,
    y = LCMS_abundance
  )
) +
  geom_point(size = 2.2) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ HPLC_target,
    scales = "free",
    nrow = 1
  ) +
  geom_text(
    data = cor_label %>%
      mutate(
        HPLC_target = factor(
          HPLC_target,
          levels = target_levels
        )
      ),
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.2
  ) +
  labs(
    x = "HPLC content (mg/g)",
    y = "LC-MS abundance"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(out_dir, "Supplementary_Figure_HPLC_LCMS_correlations.pdf"),
  p_cor,
  width = 12,
  height = 3.8
)

write.csv(
  hplc_summary,
  file.path(out_dir, "HPLC_group_summary_mg_g.csv"),
  row.names = FALSE
)

write.csv(
  lcms_summary,
  file.path(out_dir, "LCMS_group_summary_scaled.csv"),
  row.names = FALSE
)

cat("\nDone.\nOutput directory:\n", out_dir, "\n")
