# ============================================================
# 04_metabolomics_multivariate_and_differential_analysis.R
# Main untargeted-metabolomics analysis for CON / ELC / ORG
#
# IMPORTANT:
# all_metab_data.txt is treated here as the final log2-scale
# analytical matrix supplied after platform preprocessing.
# Therefore no additional log2(x + 1) transformation is applied.
# ============================================================

rm(list = ls())
gc()

library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(permute)
library(limma)
library(ComplexUpset)
library(pheatmap)
library(emmeans)
library(patchwork)
library(grid)

set.seed(20261106)

# ============================================================
# 1. Paths
# ============================================================

project_dir <- "F:/科研/有机生态低碳茶叶分析"

data_file <- file.path(
  project_dir,
  "代谢组结果/1.Preprocess/mix/all_metab_data.txt"
)

farm_file <- file.path(
  project_dir,
  "代谢组结果/order.xlsx"
)

annotation_file <- file.path(
  project_dir,
  "Paper/Revised2/01_annotation/Supplementary_Table_S2_revised.csv"
)

out_dir <- file.path(
  project_dir,
  "Paper/Revised2/04_metabolomics"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read metabolomics and annotation-confidence data
# ============================================================

metab_df <- fread(
  data_file,
  data.table = FALSE,
  check.names = FALSE
)

if (file.exists(annotation_file)) {
  anno_review <- fread(
    annotation_file,
    data.table = FALSE,
    check.names = FALSE
  ) %>%
    select(
      metab_id,
      MSI_level,
      Annotation_status,
      Annotation_evidence_type
    )

  metab_df <- metab_df %>%
    left_join(anno_review, by = "metab_id")
} else {
  warning("Revised annotation file not found. Run 01_metabolomics_annotation_confidence.R first.")
  metab_df$MSI_level <- NA_character_
  metab_df$Annotation_status <- NA_character_
  metab_df$Annotation_evidence_type <- NA_character_
}

# ============================================================
# 3. Farm metadata
# ============================================================

order_df <- read_excel(
  farm_file,
  col_names = FALSE
)

colnames(order_df)[1:2] <- c("Sample", "Farm_code")

order_df <- order_df %>%
  transmute(
    Sample = as.character(Sample),
    Farm = paste0("Farm ", str_extract(Farm_code, "\\d+$")),
    FarmBlock = paste0("F", str_extract(Farm_code, "\\d+$"))
  )

bio_cols <- grep(
  "^(CON|ELC|ORG)_",
  colnames(metab_df),
  value = TRUE
)

bio_info <- data.frame(
  Sample = bio_cols,
  Group = case_when(
    grepl("^CON_", bio_cols) ~ "CON",
    grepl("^ELC_", bio_cols) ~ "ELC",
    grepl("^ORG_", bio_cols) ~ "ORG"
  ),
  stringsAsFactors = FALSE
) %>%
  left_join(order_df, by = "Sample") %>%
  mutate(
    Group = factor(Group, levels = c("CON", "ELC", "ORG")),
    Farm = factor(Farm, levels = c("Farm 1", "Farm 2", "Farm 3")),
    FarmBlock = factor(FarmBlock, levels = c("F1", "F2", "F3"))
  )

stopifnot(!any(is.na(bio_info$FarmBlock)))

write.csv(
  bio_info,
  file.path(out_dir, "sample_metadata.csv"),
  row.names = FALSE
)

print(table(bio_info$Group, bio_info$Farm))

# ============================================================
# 4. Expression matrix
# ============================================================

expr <- as.matrix(
  metab_df[, bio_info$Sample, drop = FALSE]
)

mode(expr) <- "numeric"
rownames(expr) <- metab_df$metab_id

if (any(!is.finite(expr))) {
  stop("Non-finite values detected in final analytical matrix.")
}

# No additional logarithmic transformation here.
expr_log2 <- expr

# ============================================================
# 5. PCA
# ============================================================

pca_res <- prcomp(
  t(expr_log2),
  center = TRUE,
  scale. = TRUE
)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
) %>%
  left_join(bio_info, by = "Sample")

group_colors <- c(
  "CON" = "#E64B35",
  "ELC" = "#00A087",
  "ORG" = "#4DBBD5"
)

farm_shapes <- c(
  "Farm 1" = 16,
  "Farm 2" = 17,
  "Farm 3" = 15
)

p_fig2a <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Farm
  )
) +
  geom_point(size = 3.5) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = farm_shapes) +
  labs(
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
    color = NULL,
    shape = "Farm"
  ) +
  theme_bw(base_size = 12)
p_fig2a
ggsave(
  file.path(out_dir, "Figure2A_PCA_farm.pdf"),
  p_fig2a,
  width = 6.5,
  height = 5.2
)

write.csv(
  pca_df,
  file.path(out_dir, "Figure2A_PCA_scores.csv"),
  row.names = FALSE
)

# ============================================================
# 6. CON_7-CON_10 sensitivity dataset and PCA
# ============================================================

drop_samples <- intersect(
  paste0("CON_", 7:10),
  colnames(expr_log2)
)

keep_samples <- setdiff(
  colnames(expr_log2),
  drop_samples
)

expr_keep <- expr_log2[
  ,
  keep_samples,
  drop = FALSE
]

bio_info_keep <- bio_info %>%
  filter(
    Sample %in% keep_samples
  )

bio_info_keep <- bio_info_keep[
  match(
    keep_samples,
    bio_info_keep$Sample
  ),
]

stopifnot(
  identical(
    colnames(expr_keep),
    bio_info_keep$Sample
  )
)

pca_keep <- prcomp(
  t(expr_keep),
  center = TRUE,
  scale. = TRUE
)

pca_keep_var <- pca_keep$sdev^2 / sum(pca_keep$sdev^2) * 100

pca_keep_df <- data.frame(
  Sample = rownames(pca_keep$x),
  PC1 = pca_keep$x[, 1],
  PC2 = pca_keep$x[, 2]
) %>%
  left_join(
    bio_info_keep,
    by = "Sample"
  )

p_pca_sensitivity <- ggplot(
  pca_keep_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Farm,
    label = Sample
  )
) +
  geom_point(size = 3.2) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = farm_shapes) +
  labs(
    title = "Without CON_7-CON_10",
    x = paste0("PC1 (", round(pca_keep_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_keep_var[2], 1), "%)"),
    color = NULL,
    shape = "Farm"
  ) +
  theme_bw(base_size = 12)

p_pca_sensitivity

ggsave(
  file.path(out_dir, "Sensitivity_PCA_without_CON_7_CON_10.pdf"),
  p_pca_sensitivity,
  width = 8,
  height = 6
)

# ============================================================
# 7. Distance matrices
# ============================================================

dist_all <- dist(
  t(expr_log2),
  method = "euclidean"
)

dist_keep <- dist(
  t(expr_keep),
  method = "euclidean"
)

# ============================================================
# 8. PERMANOVA robustness analyses
# ============================================================

# Group only, unrestricted permutation
adonis_free <- adonis2(
  dist_all ~ Group,
  data = bio_info,
  permutations = 9999
)

# Without CON_7-CON_10
adonis_keep <- adonis2(
  dist_keep ~ Group,
  data = bio_info_keep,
  permutations = 9999
)

# Group only, permutation restricted within farm
perm_farm <- how(nperm = 9999)
setBlocks(perm_farm) <- bio_info$FarmBlock

adonis_blocked <- adonis2(
  dist_all ~ Group,
  data = bio_info,
  permutations = perm_farm
)

# Farm-adjusted model, permutation restricted within farm
adonis_farm <- adonis2(
  dist_all ~ FarmBlock + Group,
  data = bio_info,
  permutations = perm_farm,
  by = "margin"
)

permanova_summary <- data.frame(
  Model = c(
    "Group only, unrestricted permutation",
    "Group only, within-farm permutation",
    "Farm-adjusted, within-farm permutation"
  ),
  R2 = c(
    adonis_free["Model", "R2"],
    adonis_blocked["Model", "R2"],
    adonis_farm["Group", "R2"]
  ),
  F = c(
    adonis_free["Model", "F"],
    adonis_blocked["Model", "F"],
    adonis_farm["Group", "F"]
  ),
  P = c(
    adonis_free["Model", "Pr(>F)"],
    adonis_blocked["Model", "Pr(>F)"],
    adonis_farm["Group", "Pr(>F)"]
  )
)

print(permanova_summary)

write.csv(
  permanova_summary,
  file.path(out_dir, "PERMANOVA_summary.csv"),
  row.names = FALSE
)

permanova_sensitivity <- data.frame(
  Analysis = c(
    "All samples",
    "Without CON_7-CON_10"
  ),
  R2 = c(
    adonis_free["Model", "R2"],
    adonis_keep["Model", "R2"]
  ),
  F = c(
    adonis_free["Model", "F"],
    adonis_keep["Model", "F"]
  ),
  P = c(
    adonis_free["Model", "Pr(>F)"],
    adonis_keep["Model", "Pr(>F)"]
  )
)

print(permanova_sensitivity)

write.csv(
  permanova_sensitivity,
  file.path(out_dir, "PERMANOVA_sensitivity_CON7_CON10.csv"),
  row.names = FALSE
)

# ============================================================
# 9. Multivariate dispersion robustness analyses
# ============================================================

# All samples
bd_all <- betadisper(
  dist_all,
  group = bio_info$Group
)

bd_all_anova <- anova(bd_all)

bd_all_perm <- permutest(
  bd_all,
  permutations = 9999
)

# Without CON_7-CON_10
bd_keep <- betadisper(
  dist_keep,
  group = bio_info_keep$Group
)

bd_keep_anova <- anova(bd_keep)

bd_keep_perm <- permutest(
  bd_keep,
  permutations = 9999
)

# Permutations restricted within farm
perm_farm_disp <- how(nperm = 9999)
setBlocks(perm_farm_disp) <- bio_info$FarmBlock

bd_farm_perm <- permutest(
  bd_all,
  permutations = perm_farm_disp
)

print(bd_all_anova)
print(bd_all_perm)
print(bd_keep_anova)
print(bd_keep_perm)
print(bd_farm_perm)

# Distance-to-centroid data
disp_df <- data.frame(
  Sample = names(bd_all$distances),
  DistanceToCentroid = bd_all$distances,
  stringsAsFactors = FALSE
) %>%
  left_join(
    bio_info %>%
      select(
        Sample,
        Group,
        Farm,
        FarmBlock
      ),
    by = "Sample"
  )

dispersion_summary <- disp_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean = mean(DistanceToCentroid),
    median = median(DistanceToCentroid),
    sd = sd(DistanceToCentroid),
    .groups = "drop"
  )

print(dispersion_summary)

write.csv(
  dispersion_summary,
  file.path(out_dir, "Dispersion_summary.csv"),
  row.names = FALSE
)

# Figure 2B
p_fig2b <- ggplot(
  disp_df,
  aes(
    x = Group,
    y = DistanceToCentroid,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2
  ) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = NULL,
    y = "Distance to centroid"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none"
  )

p_fig2b

ggsave(
  file.path(out_dir, "Figure2B_distance_to_centroid.pdf"),
  p_fig2b,
  width = 5,
  height = 5
)

ggsave(
  file.path(out_dir, "Figure2B_distance_to_centroid.png"),
  p_fig2b,
  width = 5,
  height = 5,
  dpi = 600
)

# Farm-adjusted dispersion model
disp_model_original <- lm(
  DistanceToCentroid ~ Group,
  data = disp_df
)

disp_model_farm <- lm(
  DistanceToCentroid ~ FarmBlock + Group,
  data = disp_df
)

disp_original_anova <- anova(disp_model_original)
disp_farm_anova <- anova(disp_model_farm)

disp_pairwise_farm <- emmeans(
  disp_model_farm,
  pairwise ~ Group,
  adjust = "tukey"
)

print(disp_original_anova)
print(disp_farm_anova)
print(disp_pairwise_farm)

write.csv(
  as.data.frame(disp_original_anova),
  file.path(out_dir, "Dispersion_original_ANOVA.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(disp_farm_anova),
  file.path(out_dir, "Dispersion_farm_adjusted_ANOVA.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(disp_pairwise_farm$emmeans),
  file.path(out_dir, "Dispersion_farm_adjusted_emmeans.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(disp_pairwise_farm$contrasts),
  file.path(out_dir, "Dispersion_farm_adjusted_pairwise_Tukey.csv"),
  row.names = FALSE
)

# ============================================================
# 10. limma differential analysis: group-only model
# ============================================================

design <- model.matrix(
  ~ 0 + Group,
  data = bio_info
)

colnames(design) <- levels(bio_info$Group)

contrast_matrix <- makeContrasts(
  CON_vs_ELC = CON - ELC,
  CON_vs_ORG = CON - ORG,
  ELC_vs_ORG = ELC - ORG,
  levels = design
)

fit <- lmFit(
  expr_log2,
  design
)

fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- eBayes(fit2)

get_limma <- function(coef_name) {
  x <- topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  x$metab_id <- rownames(x)
  x
}

res_CON_ELC <- get_limma("CON_vs_ELC")
res_CON_ORG <- get_limma("CON_vs_ORG")
res_ELC_ORG <- get_limma("ELC_vs_ORG")

# ============================================================
# 11. Farm-adjusted limma model
# ============================================================

design_farm <- model.matrix(
  ~ 0 + Group + FarmBlock,
  data = bio_info
)

contrast_matrix_farm <- makeContrasts(
  CON_vs_ELC = GroupCON - GroupELC,
  CON_vs_ORG = GroupCON - GroupORG,
  ELC_vs_ORG = GroupELC - GroupORG,
  levels = design_farm
)

fit_farm <- lmFit(
  expr_log2,
  design_farm
)

fit2_farm <- contrasts.fit(
  fit_farm,
  contrast_matrix_farm
)

fit2_farm <- eBayes(fit2_farm)

get_limma_farm <- function(coef_name) {
  x <- topTable(
    fit2_farm,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  x$metab_id <- rownames(x)
  x
}

res_CON_ELC_farm <- get_limma_farm("CON_vs_ELC")
res_CON_ORG_farm <- get_limma_farm("CON_vs_ORG")
res_ELC_ORG_farm <- get_limma_farm("ELC_vs_ORG")

# ============================================================
# 12. Differential-feature thresholds and robustness counts
# ============================================================

p_cutoff <- 0.05
fdr05_cutoff <- 0.05
fdr10_cutoff <- 0.10
logfc_cutoff <- 0.58  # approximately 1.5-fold on log2 scale

count_features <- function(res, p_col, p_thr) {
  sum(
    res[[p_col]] < p_thr &
      abs(res$logFC) > logfc_cutoff,
    na.rm = TRUE
  )
}

diff_count_summary <- data.frame(
  Comparison = c(
    "CON vs ELC",
    "CON vs ORG",
    "ELC vs ORG"
  ),
  Nominal_P005 = c(
    count_features(res_CON_ELC, "P.Value", 0.05),
    count_features(res_CON_ORG, "P.Value", 0.05),
    count_features(res_ELC_ORG, "P.Value", 0.05)
  ),
  Farm_adjusted_Nominal_P005 = c(
    count_features(res_CON_ELC_farm, "P.Value", 0.05),
    count_features(res_CON_ORG_farm, "P.Value", 0.05),
    count_features(res_ELC_ORG_farm, "P.Value", 0.05)
  ),
  FDR005 = c(
    count_features(res_CON_ELC, "adj.P.Val", 0.05),
    count_features(res_CON_ORG, "adj.P.Val", 0.05),
    count_features(res_ELC_ORG, "adj.P.Val", 0.05)
  ),
  Farm_adjusted_FDR005 = c(
    count_features(res_CON_ELC_farm, "adj.P.Val", 0.05),
    count_features(res_CON_ORG_farm, "adj.P.Val", 0.05),
    count_features(res_ELC_ORG_farm, "adj.P.Val", 0.05)
  ),
  FDR010 = c(
    count_features(res_CON_ELC, "adj.P.Val", 0.10),
    count_features(res_CON_ORG, "adj.P.Val", 0.10),
    count_features(res_ELC_ORG, "adj.P.Val", 0.10)
  ),
  Farm_adjusted_FDR010 = c(
    count_features(res_CON_ELC_farm, "adj.P.Val", 0.10),
    count_features(res_CON_ORG_farm, "adj.P.Val", 0.10),
    count_features(res_ELC_ORG_farm, "adj.P.Val", 0.10)
  )
)

print(diff_count_summary)

write.csv(
  diff_count_summary,
  file.path(out_dir, "differential_feature_count_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 13. Original vs farm-adjusted logFC agreement
# ============================================================

compare_logFC <- function(old_res, farm_res) {

  tmp <- inner_join(
    old_res %>%
      select(
        metab_id,
        logFC_original = logFC
      ),
    farm_res %>%
      select(
        metab_id,
        logFC_farm_adjusted = logFC
      ),
    by = "metab_id"
  )

  data.frame(
    Pearson = cor(
      tmp$logFC_original,
      tmp$logFC_farm_adjusted,
      use = "complete.obs",
      method = "pearson"
    ),
    Spearman = cor(
      tmp$logFC_original,
      tmp$logFC_farm_adjusted,
      use = "complete.obs",
      method = "spearman"
    )
  )
}

logFC_compare <- bind_rows(
  "CON vs ELC" = compare_logFC(
    res_CON_ELC,
    res_CON_ELC_farm
  ),
  "CON vs ORG" = compare_logFC(
    res_CON_ORG,
    res_CON_ORG_farm
  ),
  "ELC vs ORG" = compare_logFC(
    res_ELC_ORG,
    res_ELC_ORG_farm
  ),
  .id = "Comparison"
)

print(logFC_compare)

write.csv(
  logFC_compare,
  file.path(out_dir, "limma_logFC_agreement.csv"),
  row.names = FALSE
)

# ============================================================
# 14. Farm-stratified label permutation for differential counts
# ============================================================

n_perm <- 9999

count_from_fit <- function(expr_mat, group_vector) {
  group_factor <- factor(
    group_vector,
    levels = c("CON", "ELC", "ORG")
  )

  design_p <- model.matrix(~ 0 + group_factor)
  colnames(design_p) <- c("CON", "ELC", "ORG")

  contrast_p <- makeContrasts(
    CON_vs_ELC = CON - ELC,
    CON_vs_ORG = CON - ORG,
    ELC_vs_ORG = ELC - ORG,
    levels = design_p
  )

  fit_p <- lmFit(expr_mat, design_p)
  fit_p <- contrasts.fit(fit_p, contrast_p)
  fit_p <- eBayes(fit_p)

  sapply(
    1:3,
    function(k) {
      tt <- topTable(
        fit_p,
        coef = k,
        number = Inf,
        adjust.method = "BH",
        sort.by = "none"
      )
      sum(
        tt$P.Value < p_cutoff &
          abs(tt$logFC) > logfc_cutoff,
        na.rm = TRUE
      )
    }
  )
}

observed_counts <- c(
  count_features(res_CON_ELC, "P.Value", 0.05),
  count_features(res_CON_ORG, "P.Value", 0.05),
  count_features(res_ELC_ORG, "P.Value", 0.05)
)

null_counts <- matrix(
  NA_integer_,
  nrow = n_perm,
  ncol = 3,
  dimnames = list(
    NULL,
    c("CON_vs_ELC", "CON_vs_ORG", "ELC_vs_ORG")
  )
)

for (b in seq_len(n_perm)) {
  group_perm <- as.character(bio_info$Group)

  for (farm in levels(bio_info$FarmBlock)) {
    idx <- which(bio_info$FarmBlock == farm)
    group_perm[idx] <- sample(group_perm[idx], length(idx), replace = FALSE)
  }

  null_counts[b, ] <- count_from_fit(expr_log2, group_perm)

  if (b %% 500 == 0) {
    cat("Differential-feature permutation:", b, "/", n_perm, "\n")
  }
}

perm_summary <- data.frame(
  Comparison = c("CON vs ELC", "CON vs ORG", "ELC vs ORG"),
  Observed = observed_counts,
  Null_mean = colMeans(null_counts),
  Null_SD = apply(null_counts, 2, sd),
  Null_median = apply(null_counts, 2, median),
  Null_95th = apply(null_counts, 2, quantile, probs = 0.95),
  Null_99th = apply(null_counts, 2, quantile, probs = 0.99),
  Empirical_P = sapply(
    1:3,
    function(j) {
      (1 + sum(null_counts[, j] >= observed_counts[j])) /
        (n_perm + 1)
    }
  )
)

print(perm_summary)

write.csv(
  as.data.frame(null_counts),
  file.path(out_dir, "farm_stratified_permutation_null_counts.csv"),
  row.names = FALSE
)

write.csv(
  perm_summary,
  file.path(out_dir, "farm_stratified_permutation_summary.csv"),
  row.names = FALSE
)

perm_long <- as.data.frame(null_counts) %>%
  pivot_longer(
    everything(),
    names_to = "Comparison",
    values_to = "Count"
  ) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c("CON_vs_ELC", "CON_vs_ORG", "ELC_vs_ORG"),
      labels = c("CON vs ELC", "CON vs ORG", "ELC vs ORG")
    )
  )

obs_df <- data.frame(
  Comparison = factor(
    c("CON vs ELC", "CON vs ORG", "ELC vs ORG"),
    levels = c("CON vs ELC", "CON vs ORG", "ELC vs ORG")
  ),
  Observed = observed_counts
)

p_perm <- ggplot(
  perm_long,
  aes(x = Count)
) +
  geom_histogram(
    bins = 40,
    color = "black",
    fill = "grey80"
  ) +
  geom_vline(
    data = obs_df,
    aes(xintercept = Observed),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(~ Comparison, scales = "free_x", nrow = 1) +
  labs(
    x = "Differential-feature count under permuted labels",
    y = "Frequency"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(out_dir, "Supplementary_Figure_farm_stratified_permutation.pdf"),
  p_perm,
  width = 11,
  height = 4
)

# ============================================================
# 15. Figure 2C: UpSet plot using exploratory threshold
# ============================================================

sig_CON_ELC <- res_CON_ELC %>%
  filter(
    P.Value < p_cutoff,
    abs(logFC) > logfc_cutoff
  ) %>%
  pull(metab_id)

sig_CON_ORG <- res_CON_ORG %>%
  filter(
    P.Value < p_cutoff,
    abs(logFC) > logfc_cutoff
  ) %>%
  pull(metab_id)

sig_ELC_ORG <- res_ELC_ORG %>%
  filter(
    P.Value < p_cutoff,
    abs(logFC) > logfc_cutoff
  ) %>%
  pull(metab_id)

all_sig <- unique(
  c(sig_CON_ELC, sig_CON_ORG, sig_ELC_ORG)
)

upset_df <- data.frame(
  metab_id = all_sig,
  CON_vs_ELC = all_sig %in% sig_CON_ELC,
  CON_vs_ORG = all_sig %in% sig_CON_ORG,
  ELC_vs_ORG = all_sig %in% sig_ELC_ORG
)

p_fig2c <- upset(
  upset_df,
  intersect = c(
    "CON_vs_ELC",
    "CON_vs_ORG",
    "ELC_vs_ORG"
  ),
  name = "Differential features",
  sort_intersections_by = "cardinality"
)

p_fig2c

ggsave(
  file.path(out_dir, "Figure2C_UpSet.pdf"),
  p_fig2c,
  width = 7,
  height = 5.5
)

# ============================================================
# 16. Figure 2D
# High-confidence annotated differential metabolites
# ============================================================


# ============================================================
# 16.1 Retain named MSI Level 2 annotations
# ============================================================

annotated_ids <- metab_df %>%
  filter(
    !is.na(Metabolite),
    Metabolite != "",
    !grepl(
      "^meta",
      Metabolite,
      ignore.case = TRUE
    ),
    MSI_level == "Level 2"
  ) %>%
  pull(metab_id)

cat(
  "\nNamed MSI Level 2 features:",
  length(annotated_ids),
  "\n"
)


# ============================================================
# 16.2 Calculate mean abundance by production system
# ============================================================

mean_df <- data.frame(
  metab_id = rownames(expr_log2),
  
  mean_CON = rowMeans(
    expr_log2[
      ,
      grepl(
        "^CON_",
        colnames(expr_log2)
      ),
      drop = FALSE
    ],
    na.rm = TRUE
  ),
  
  mean_ELC = rowMeans(
    expr_log2[
      ,
      grepl(
        "^ELC_",
        colnames(expr_log2)
      ),
      drop = FALSE
    ],
    na.rm = TRUE
  ),
  
  mean_ORG = rowMeans(
    expr_log2[
      ,
      grepl(
        "^ORG_",
        colnames(expr_log2)
      ),
      drop = FALSE
    ],
    na.rm = TRUE
  )
)


# ============================================================
# 16.3 Extract farm-adjusted limma statistics
# ============================================================

r_CON_ELC <- res_CON_ELC_farm %>%
  filter(
    metab_id %in% annotated_ids
  ) %>%
  select(
    metab_id,
    FDR_CON_ELC = adj.P.Val,
    FC_CON_ELC = logFC
  )


r_CON_ORG <- res_CON_ORG_farm %>%
  filter(
    metab_id %in% annotated_ids
  ) %>%
  select(
    metab_id,
    FDR_CON_ORG = adj.P.Val,
    FC_CON_ORG = logFC
  )


r_ELC_ORG <- res_ELC_ORG_farm %>%
  filter(
    metab_id %in% annotated_ids
  ) %>%
  select(
    metab_id,
    FDR_ELC_ORG = adj.P.Val,
    FC_ELC_ORG = logFC
  )


heatmap_candidates <- r_CON_ELC %>%
  full_join(
    r_CON_ORG,
    by = "metab_id"
  ) %>%
  full_join(
    r_ELC_ORG,
    by = "metab_id"
  ) %>%
  left_join(
    mean_df,
    by = "metab_id"
  ) %>%
  left_join(
    metab_df %>%
      select(
        metab_id,
        Metabolite,
        MSI_level,
        Annotation_status,
        Annotation_evidence_type
      ),
    by = "metab_id"
  )


# ============================================================
# 16.4 Define Figure 2D selection thresholds
# ============================================================

heatmap_fdr_cutoff <- 0.10
heatmap_logfc_cutoff <- 0.58


# ============================================================
# 16.5 Build production-system-specific candidate pools
#
# support:
# number of relevant pairwise comparisons satisfying
# BH-FDR < 0.10 and |log2FC| > 0.58
#
# dominance:
# abundance in the target system minus the larger mean
# abundance of the other two systems
# ============================================================

make_pool <- function(
    target_group
) {
  
  if (
    target_group == "CON"
  ) {
    
    out <- heatmap_candidates %>%
      mutate(
        
        sig_vs_ELC =
          FDR_CON_ELC < heatmap_fdr_cutoff &
          FC_CON_ELC > heatmap_logfc_cutoff,
        
        sig_vs_ORG =
          FDR_CON_ORG < heatmap_fdr_cutoff &
          FC_CON_ORG > heatmap_logfc_cutoff,
        
        support =
          as.integer(sig_vs_ELC) +
          as.integer(sig_vs_ORG),
        
        best_FDR = pmin(
          FDR_CON_ELC,
          FDR_CON_ORG,
          na.rm = TRUE
        ),
        
        effect_score =
          pmax(
            FC_CON_ELC,
            0
          ) +
          pmax(
            FC_CON_ORG,
            0
          ),
        
        dominance =
          mean_CON -
          pmax(
            mean_ELC,
            mean_ORG
          ),
        
        GroupRep = "CON"
      )
    
  } else if (
    target_group == "ELC"
  ) {
    
    out <- heatmap_candidates %>%
      mutate(
        
        sig_vs_CON =
          FDR_CON_ELC < heatmap_fdr_cutoff &
          FC_CON_ELC < -heatmap_logfc_cutoff,
        
        sig_vs_ORG =
          FDR_ELC_ORG < heatmap_fdr_cutoff &
          FC_ELC_ORG > heatmap_logfc_cutoff,
        
        support =
          as.integer(sig_vs_CON) +
          as.integer(sig_vs_ORG),
        
        best_FDR = pmin(
          FDR_CON_ELC,
          FDR_ELC_ORG,
          na.rm = TRUE
        ),
        
        effect_score =
          pmax(
            -FC_CON_ELC,
            0
          ) +
          pmax(
            FC_ELC_ORG,
            0
          ),
        
        dominance =
          mean_ELC -
          pmax(
            mean_CON,
            mean_ORG
          ),
        
        GroupRep = "ELC"
      )
    
  } else if (
    target_group == "ORG"
  ) {
    
    out <- heatmap_candidates %>%
      mutate(
        
        sig_vs_CON =
          FDR_CON_ORG < heatmap_fdr_cutoff &
          FC_CON_ORG < -heatmap_logfc_cutoff,
        
        sig_vs_ELC =
          FDR_ELC_ORG < heatmap_fdr_cutoff &
          FC_ELC_ORG < -heatmap_logfc_cutoff,
        
        support =
          as.integer(sig_vs_CON) +
          as.integer(sig_vs_ELC),
        
        best_FDR = pmin(
          FDR_CON_ORG,
          FDR_ELC_ORG,
          na.rm = TRUE
        ),
        
        effect_score =
          pmax(
            -FC_CON_ORG,
            0
          ) +
          pmax(
            -FC_ELC_ORG,
            0
          ),
        
        dominance =
          mean_ORG -
          pmax(
            mean_CON,
            mean_ELC
          ),
        
        GroupRep = "ORG"
      )
    
  } else {
    
    stop(
      "target_group must be CON, ELC, or ORG."
    )
  }
  
  out
}


# ============================================================
# 16.6 Retain enriched metabolites
# ============================================================

con_pool <- make_pool(
  "CON"
) %>%
  filter(
    support >= 1,
    dominance > 0
  ) %>%
  arrange(
    desc(support),
    best_FDR,
    desc(effect_score),
    desc(dominance)
  )


elc_pool <- make_pool(
  "ELC"
) %>%
  filter(
    support >= 1,
    dominance > 0
  ) %>%
  arrange(
    desc(support),
    best_FDR,
    desc(effect_score),
    desc(dominance)
  )


org_pool <- make_pool(
  "ORG"
) %>%
  filter(
    support >= 1,
    dominance > 0
  ) %>%
  arrange(
    desc(support),
    best_FDR,
    desc(effect_score),
    desc(dominance)
  )


cat(
  "\nFigure 2D candidate numbers:\n"
)

cat(
  "CON:",
  nrow(con_pool),
  "\n"
)

cat(
  "ELC:",
  nrow(elc_pool),
  "\n"
)

cat(
  "ORG:",
  nrow(org_pool),
  "\n"
)


# ============================================================
# 16.7 Combine all metabolites passing the criteria
# ============================================================

selected_metabolites <- bind_rows(
  con_pool,
  elc_pool,
  org_pool
)


stopifnot(
  !anyDuplicated(
    selected_metabolites$metab_id
  )
)


cat(
  "Total selected metabolites:",
  nrow(selected_metabolites),
  "\n"
)


print(
  table(
    selected_metabolites$GroupRep
  )
)


write.csv(
  selected_metabolites,
  file.path(
    out_dir,
    "Figure2D_selected_features.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16.8 Define sample order: CON -> ELC -> ORG
# ============================================================

order_group_samples <- function(
    sample_names,
    group_name
) {
  
  x <- sample_names[
    grepl(
      paste0("^", group_name, "_"),
      sample_names
    )
  ]
  
  x <- as.character(x)
  
  x[
    order(
      as.integer(
        sub(
          paste0("^", group_name, "_"),
          "",
          x
        )
      )
    )
  ]
}


sample_order <- c(
  order_group_samples(
    colnames(expr_log2),
    "CON"
  ),
  order_group_samples(
    colnames(expr_log2),
    "ELC"
  ),
  order_group_samples(
    colnames(expr_log2),
    "ORG"
  )
)

sample_order <- as.character(
  sample_order
)


# Check order explicitly
cat(
  "\nFigure 2D sample order:\n"
)

print(
  sample_order
)


stopifnot(
  all(
    grepl(
      "^CON_",
      sample_order[
        1:sum(grepl("^CON_", sample_order))
      ]
    )
  )
)

stopifnot(
  identical(
    unique(
      sub(
        "_\\d+$",
        "",
        sample_order
      )
    ),
    c(
      "CON",
      "ELC",
      "ORG"
    )
  )
)


# 16.9 Extract abundance matrix
mat_fig2d <- expr_log2[
  selected_metabolites$metab_id,
  sample_order,
  drop = FALSE
]

# 16.10 Row-wise z-score
mat_fig2d_z <- t(
  scale(
    t(mat_fig2d)
  )
)

mat_fig2d_z[
  !is.finite(mat_fig2d_z)
] <- 0

# 16.11 Replace metab_id with metabolite names
metabolite_labels <- selected_metabolites$Metabolite[
  match(
    rownames(mat_fig2d_z),
    selected_metabolites$metab_id
  )
]

stopifnot(
  !any(is.na(metabolite_labels))
)

rownames(mat_fig2d_z) <- make.unique(
  metabolite_labels
)

# ============================================================
# 16.12 Column annotation
# ============================================================

annotation_col_fig2d <- data.frame(
  Group = factor(
    sub(
      "_\\d+$",
      "",
      sample_order
    ),
    levels = c(
      "CON",
      "ELC",
      "ORG"
    )
  )
)

rownames(
  annotation_col_fig2d
) <- sample_order





# ============================================================
# 16.13 Row annotation
# ============================================================

annotation_row_fig2d <- data.frame(
  Representative_group = factor(
    selected_metabolites$GroupRep,
    levels = c(
      "CON",
      "ELC",
      "ORG"
    )
  )
)


rownames(
  annotation_row_fig2d
) <- rownames(
  mat_fig2d_z
)


# ============================================================
# 16.14 Annotation colors
# ============================================================

annotation_colors_fig2d <- list(
  Group = group_colors,
  Representative_group = group_colors
)


# ============================================================
# 16.15 Define row and column gaps
# ============================================================

n_CON_fig2d <- nrow(
  con_pool
)

n_ELC_fig2d <- nrow(
  elc_pool
)

n_ORG_fig2d <- nrow(
  org_pool
)


row_gaps_fig2d <- c(
  n_CON_fig2d,
  n_CON_fig2d + n_ELC_fig2d
)


column_groups_fig2d <- as.character(
  annotation_col_fig2d$Group
)


col_gaps_fig2d <- which(
  column_groups_fig2d[
    -length(
      column_groups_fig2d
    )
  ] !=
    column_groups_fig2d[
      -1
    ]
)


# ============================================================
# 16.16 Export Figure 2D matrices
# ============================================================

fig2d_abundance_export <- data.frame(
  metab_id = selected_metabolites$metab_id,
  Metabolite = metabolite_labels,
  GroupRep = selected_metabolites$GroupRep,
  mat_fig2d,
  check.names = FALSE
)


write.csv(
  fig2d_abundance_export,
  file.path(
    out_dir,
    "Figure2D_log2_abundance_matrix.csv"
  ),
  row.names = FALSE
)


fig2d_zscore_export <- data.frame(
  metab_id = selected_metabolites$metab_id,
  Metabolite = rownames(mat_fig2d_z),
  GroupRep = selected_metabolites$GroupRep,
  mat_fig2d_z,
  check.names = FALSE
)


write.csv(
  fig2d_zscore_export,
  file.path(
    out_dir,
    "Figure2D_zscore_matrix.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16.17 Draw Figure 2D in RStudio
# ============================================================

p_fig2d <- pheatmap(
  mat_fig2d_z,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  annotation_col = annotation_col_fig2d,
  annotation_row = annotation_row_fig2d,
  annotation_colors = annotation_colors_fig2d,
  
  gaps_row = row_gaps_fig2d,
  gaps_col = col_gaps_fig2d,
  
  border_color = NA,
  
  fontsize = 10,
  fontsize_row = 9,
  
  show_colnames = FALSE,
  show_rownames = TRUE
)
p_fig2d

# ============================================================
# 16.18 Save Figure 2D as PDF
# ============================================================

pheatmap(
  mat_fig2d_z,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  annotation_col = annotation_col_fig2d,
  annotation_row = annotation_row_fig2d,
  
  annotation_colors = annotation_colors_fig2d,
  
  gaps_row = row_gaps_fig2d,
  gaps_col = col_gaps_fig2d,
  
  border_color = NA,
  
  fontsize = 10,
  fontsize_row = 9,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  
  filename = file.path(
    out_dir,
    "Figure2D_heatmap.pdf"
  ),
  
  width = 9,
  height = 7
)


# ============================================================
# 16.19 Save Figure 2D as PNG
# ============================================================

pheatmap(
  mat_fig2d_z,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  annotation_col = annotation_col_fig2d,
  annotation_row = annotation_row_fig2d,
  
  annotation_colors = annotation_colors_fig2d,
  
  gaps_row = row_gaps_fig2d,
  gaps_col = col_gaps_fig2d,
  
  border_color = NA,
  
  fontsize = 10,
  fontsize_row = 9,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  main = NULL,
  
  filename = file.path(
    out_dir,
    "Figure2D_heatmap.png"
  ),
  
  width = 9,
  height = 7
)


# ============================================================
# 16.20 Figure 2D summary
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "Figure 2D completed\n"
)

cat(
  "CON:",
  n_CON_fig2d,
  "\n"
)

cat(
  "ELC:",
  n_ELC_fig2d,
  "\n"
)

cat(
  "ORG:",
  n_ORG_fig2d,
  "\n"
)

cat(
  "Total:",
  nrow(selected_metabolites),
  "\n"
)

cat(
  "====================================\n"
)


print(
  selected_metabolites %>%
    select(
      metab_id,
      Metabolite,
      GroupRep,
      support,
      best_FDR,
      effect_score,
      dominance
    )
)


# ============================================================
# 17. Export complete limma result tables
# ============================================================

anno_export <- metab_df %>%
  select(
    metab_id,
    Metabolite,
    MSI_level,
    Annotation_status,
    Annotation_evidence_type
  )


export_limma <- function(
    res,
    filename
) {
  
  res %>%
    left_join(
      anno_export,
      by = "metab_id"
    ) %>%
    write.csv(
      file.path(
        out_dir,
        filename
      ),
      row.names = FALSE
    )
}


# Group-only models

export_limma(
  res_CON_ELC,
  "limma_CON_vs_ELC.csv"
)

export_limma(
  res_CON_ORG,
  "limma_CON_vs_ORG.csv"
)

export_limma(
  res_ELC_ORG,
  "limma_ELC_vs_ORG.csv"
)


# Farm-adjusted models

export_limma(
  res_CON_ELC_farm,
  "limma_CON_vs_ELC_farm_adjusted.csv"
)

export_limma(
  res_CON_ORG_farm,
  "limma_CON_vs_ORG_farm_adjusted.csv"
)

export_limma(
  res_ELC_ORG_farm,
  "limma_ELC_vs_ORG_farm_adjusted.csv"
)


# ============================================================
# 18. Final output summary
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "Analysis completed.\n"
)

cat(
  "Output directory:\n",
  out_dir,
  "\n"
)

cat(
  "====================================\n"
)

