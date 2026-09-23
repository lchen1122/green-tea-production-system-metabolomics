############################################################
# Differential metabolite tables for Result 4
# Tea metabolomics: CON / ELC / ORG
############################################################

rm(list = ls())
gc()

#==============================
# 0. packages
#==============================
library(tidyverse)
library(limma)

#==============================
# 1. input / output paths
#==============================
data_file <- "F:/科研/有机生态低碳茶叶分析/代谢组结果/1.Preprocess/mix/all_metab_data.txt"
annotation_file <- "F:/科研/有机生态低碳茶叶分析/Paper/Result3/annotation_screened_full.csv"
out_dir <- "F:/科研/有机生态低碳茶叶分析/Paper/Result4/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#==============================
# 2. parameters
#==============================
# 你的数据已经是log后的，所以这里不再log
logfc_cutoff <- 0.58
p_cutoff <- 0.05
fdr_cutoff <- 0.10

comparisons <- list(
  c("CON", "ELC"),
  c("CON", "ORG"),
  c("ELC", "ORG")
)

#==============================
# 3. read data
#==============================
metab_df <- read.delim(data_file, check.names = FALSE)
anno_df  <- read.csv(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot("metab_id" %in% colnames(metab_df))
stopifnot(all(c("metab_id", "Metabolite") %in% colnames(anno_df)))

#==============================
# 4. sample info
#==============================
sample_cols <- grep("^(CON|ELC|ORG)_", colnames(metab_df), value = TRUE)

sample_info <- data.frame(
  Sample = sample_cols,
  Group = case_when(
    grepl("^CON_", sample_cols) ~ "CON",
    grepl("^ELC_", sample_cols) ~ "ELC",
    grepl("^ORG_", sample_cols) ~ "ORG"
  ),
  stringsAsFactors = FALSE
)

sample_info$Group <- factor(sample_info$Group, levels = c("CON", "ELC", "ORG"))

#==============================
# 5. expression matrix
#    rows = metabolites
#    cols = samples
#==============================
expr_mat <- as.matrix(metab_df[, sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- metab_df$metab_id

#==============================
# 6. annotation table
#==============================
# 去重，避免 join 出问题
anno_df <- anno_df %>%
  distinct(metab_id, .keep_all = TRUE)

# 只保留常用列；没有的列会自动跳过
anno_keep_cols <- intersect(
  c(
    "metab_id", "Metabolite", "DisplayName", "annotation_screen",
    "Formula", "Adducts", "compound_source", "Score", "level"
  ),
  colnames(anno_df)
)

anno_use <- anno_df %>%
  select(all_of(anno_keep_cols))

#==============================
# 7. helper function
#==============================
run_limma_pair <- function(expr_mat, sample_info, g1, g2, anno_use,
                           logfc_cutoff = 0.58,
                           p_cutoff = 0.05,
                           fdr_cutoff = 0.10) {
  
  # subset samples
  keep_samples <- sample_info$Sample[sample_info$Group %in% c(g1, g2)]
  info_sub <- sample_info %>%
    filter(Sample %in% keep_samples) %>%
    mutate(Group = factor(Group, levels = c(g1, g2)))
  
  expr_sub <- expr_mat[, info_sub$Sample, drop = FALSE]
  
  # design matrix
  design <- model.matrix(~ 0 + Group, data = info_sub)
  colnames(design) <- levels(info_sub$Group)
  
  # limma
  fit <- lmFit(expr_sub, design)
  contrast_mat <- makeContrasts(contrasts = paste0(g1, "-", g2), levels = design)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)
  
  tt <- topTable(fit2, number = Inf, sort.by = "none") %>%
    rownames_to_column("metab_id")
  
  # group means
  mean_g1 <- rowMeans(expr_sub[, info_sub$Group == g1, drop = FALSE], na.rm = TRUE)
  mean_g2 <- rowMeans(expr_sub[, info_sub$Group == g2, drop = FALSE], na.rm = TRUE)
  
  mean_df <- data.frame(
    metab_id = rownames(expr_sub),
    mean_g1 = mean_g1,
    mean_g2 = mean_g2,
    stringsAsFactors = FALSE
  )
  
  # merge and annotate
  res <- tt %>%
    left_join(mean_df, by = "metab_id") %>%
    left_join(anno_use, by = "metab_id") %>%
    mutate(
      Comparison = paste0(g1, "_vs_", g2),
      Group1 = g1,
      Group2 = g2,
      higher_in = case_when(
        logFC > 0  ~ g1,
        logFC < 0  ~ g2,
        TRUE       ~ "Tie"
      ),
      abs_logFC = abs(logFC),
      sig_rawp = (P.Value < p_cutoff & abs_logFC >= logfc_cutoff),
      sig_fdr  = (adj.P.Val < fdr_cutoff & abs_logFC >= logfc_cutoff),
      curated_annotation = case_when(
        !is.na(annotation_screen) &
          annotation_screen %in% c("recommended", "keep_with_caution") &
          !is.na(Metabolite) &
          Metabolite != "" &
          !grepl("^metab_", Metabolite, ignore.case = TRUE) ~ TRUE,
        TRUE ~ FALSE
      ),
      curated_sig_rawp = sig_rawp & curated_annotation,
      curated_sig_fdr  = sig_fdr  & curated_annotation
    ) %>%
    select(
      metab_id, Comparison, Group1, Group2, higher_in,
      logFC, abs_logFC, AveExpr, t, P.Value, adj.P.Val, B,
      mean_g1, mean_g2,
      everything()
    )
  
  res
}

#==============================
# 8. run all pairwise comparisons
#==============================
res_list <- list()
summary_list <- list()

for (cmp in comparisons) {
  g1 <- cmp[1]
  g2 <- cmp[2]
  cmp_name <- paste0(g1, "_vs_", g2)
  
  cat("Running:", cmp_name, "\n")
  
  res <- run_limma_pair(
    expr_mat = expr_mat,
    sample_info = sample_info,
    g1 = g1,
    g2 = g2,
    anno_use = anno_use,
    logfc_cutoff = logfc_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff
  )
  
  res_list[[cmp_name]] <- res
  
  # save full results
  write.csv(
    res,
    file.path(out_dir, paste0(cmp_name, "_all_features.csv")),
    row.names = FALSE
  )
  
  # save significant by raw p
  write.csv(
    res %>% filter(sig_rawp),
    file.path(out_dir, paste0(cmp_name, "_sig_rawp.csv")),
    row.names = FALSE
  )
  
  # save significant by FDR
  write.csv(
    res %>% filter(sig_fdr),
    file.path(out_dir, paste0(cmp_name, "_sig_fdr.csv")),
    row.names = FALSE
  )
  
  # save curated significant by raw p
  write.csv(
    res %>% filter(curated_sig_rawp),
    file.path(out_dir, paste0(cmp_name, "_curated_sig_rawp.csv")),
    row.names = FALSE
  )
  
  # save curated significant by FDR
  write.csv(
    res %>% filter(curated_sig_fdr),
    file.path(out_dir, paste0(cmp_name, "_curated_sig_fdr.csv")),
    row.names = FALSE
  )
  
  summary_list[[cmp_name]] <- data.frame(
    Comparison = cmp_name,
    n_all = nrow(res),
    n_sig_rawp = sum(res$sig_rawp, na.rm = TRUE),
    n_sig_fdr = sum(res$sig_fdr, na.rm = TRUE),
    n_curated_sig_rawp = sum(res$curated_sig_rawp, na.rm = TRUE),
    n_curated_sig_fdr = sum(res$curated_sig_fdr, na.rm = TRUE)
  )
}

#==============================
# 9. comparison summary
#==============================
summary_df <- bind_rows(summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "comparison_summary.csv"),
  row.names = FALSE
)

summary_df

#==============================
# 10. union table for Result 4
#     use curated significant metabolites across all comparisons
#==============================
curated_union_long <- bind_rows(res_list) %>%
  filter(curated_sig_rawp) %>%
  distinct(Comparison, metab_id, .keep_all = TRUE)

# membership across comparisons
membership_df <- curated_union_long %>%
  mutate(flag = 1) %>%
  select(metab_id, Comparison, flag) %>%
  pivot_wider(
    names_from = Comparison,
    values_from = flag,
    values_fill = 0
  )

# summary by metabolite
curated_union_summary <- curated_union_long %>%
  group_by(metab_id) %>%
  summarise(
    Metabolite = first(Metabolite),
    DisplayName = first(DisplayName),
    annotation_screen = first(annotation_screen),
    Formula = first(Formula),
    Adducts = first(Adducts),
    compound_source = first(compound_source),
    Score = first(Score),
    level = first(level),
    n_comparisons = n(),
    comparisons = paste(sort(unique(Comparison)), collapse = "; "),
    higher_in_summary = paste(paste0(Comparison, ":", higher_in), collapse = "; "),
    min_raw_p = min(P.Value, na.rm = TRUE),
    min_adj_p = min(adj.P.Val, na.rm = TRUE),
    max_abs_logFC = max(abs_logFC, na.rm = TRUE),
    mean_abs_logFC = mean(abs_logFC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(membership_df, by = "metab_id") %>%
  arrange(desc(n_comparisons), min_raw_p)

write.csv(
  curated_union_summary,
  file.path(out_dir, "curated_sig_union_summary.csv"),
  row.names = FALSE
)

#==============================
# 11. template for manual chemical class annotation
#     this is the file you can use for Result 4
#==============================
class_template <- curated_union_summary %>%
  transmute(
    metab_id,
    Metabolite,
    DisplayName,
    annotation_screen,
    comparisons,
    higher_in_summary,
    max_abs_logFC,
    min_raw_p,
    ChemClass = "",
    ChemSubclass = "",
    Note = ""
  )

write.csv(
  class_template,
  file.path(out_dir, "chemical_class_annotation_template.csv"),
  row.names = FALSE
)

#==============================
# 12. optional: one combined long table
#==============================
all_pairwise_results <- bind_rows(res_list)

write.csv(
  all_pairwise_results,
  file.path(out_dir, "all_pairwise_results_combined.csv"),
  row.names = FALSE
)

cat("\nDone.\n")
cat("Main output files:\n")
cat("1) comparison_summary.csv\n")
cat("2) *_curated_sig_rawp.csv for each comparison\n")
cat("3) curated_sig_union_summary.csv\n")
cat("4) chemical_class_annotation_template.csv\n")



############################################################
# Result 4. Chemical class analysis
# Tea metabolomics: CON / ELC / ORG
############################################################

rm(list = ls())
gc()

#==============================
# 0. packages
#==============================
library(tidyverse)

#==============================
# 1. input / output paths
#==============================
data_file <- "F:/科研/有机生态低碳茶叶分析/代谢组结果/1.Preprocess/mix/all_metab_data.txt"
union_file <- "F:/科研/有机生态低碳茶叶分析/Paper/Result4/curated_sig_union_summary.csv"
class_file <- "F:/科研/有机生态低碳茶叶分析/Paper/Result4/chemical_class_annotation_template.csv"

out_dir <- "F:/科研/有机生态低碳茶叶分析/Paper/Result4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#==============================
# 2. read data
#==============================
metab_df <- read.delim(data_file, check.names = FALSE)
union_df <- read.csv(union_file, check.names = FALSE, stringsAsFactors = FALSE)
class_df <- read.csv(class_file, check.names = FALSE, stringsAsFactors = FALSE)

#==============================
# 3. basic checks
#==============================
stopifnot("metab_id" %in% colnames(metab_df))
stopifnot("metab_id" %in% colnames(union_df))
stopifnot("metab_id" %in% colnames(class_df))
stopifnot("ChemClass" %in% colnames(class_df))

#==============================
# 4. sample info
#==============================
sample_cols <- grep("^(CON|ELC|ORG)_", colnames(metab_df), value = TRUE)

sample_info <- data.frame(
  Sample = sample_cols,
  Group = case_when(
    grepl("^CON_", sample_cols) ~ "CON",
    grepl("^ELC_", sample_cols) ~ "ELC",
    grepl("^ORG_", sample_cols) ~ "ORG"
  ),
  stringsAsFactors = FALSE
)

sample_info$Group <- factor(sample_info$Group, levels = c("CON", "ELC", "ORG"))

#==============================
# 5. clean class table
#==============================
class_df2 <- class_df %>%
  distinct(metab_id, .keep_all = TRUE) %>%
  mutate(
    ChemClass = trimws(ChemClass),
    ChemSubclass = if ("ChemSubclass" %in% colnames(.)) trimws(ChemSubclass) else NA_character_,
    Note = if ("Note" %in% colnames(.)) Note else NA_character_
  ) %>%
  filter(!is.na(ChemClass), ChemClass != "")

#==============================
# 6. merge differential pool
#==============================
diff_pool <- union_df %>%
  distinct(metab_id, .keep_all = TRUE) %>%
  left_join(
    class_df2 %>% select(metab_id, ChemClass, ChemSubclass, Note),
    by = "metab_id"
  ) %>%
  filter(!is.na(ChemClass), ChemClass != "")

# 保存一下主分析池
write.csv(
  diff_pool,
  file.path(out_dir, "Result4_diff_pool_with_classes.csv"),
  row.names = FALSE
)

#==============================
# 7. calculate group means and dominant group
#==============================
expr_long <- metab_df %>%
  filter(metab_id %in% diff_pool$metab_id) %>%
  select(metab_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(sample_info, by = "Sample")

group_mean_df <- expr_long %>%
  group_by(metab_id, Group) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_abundance,
    names_prefix = "mean_"
  ) %>%
  rowwise() %>%
  mutate(
    DominantGroup = c("CON", "ELC", "ORG")[which.max(c(mean_CON, mean_ELC, mean_ORG))]
  ) %>%
  ungroup()

diff_pool2 <- diff_pool %>%
  left_join(group_mean_df, by = "metab_id")

write.csv(
  diff_pool2,
  file.path(out_dir, "Result4_diff_pool_with_classes_and_group_means.csv"),
  row.names = FALSE
)

#==============================
# 8. Figure 4A data
#==============================
# 统计每个 dominant group 下各 chemical class 的数量

fig4a_df <- diff_pool2 %>%
  count(DominantGroup, ChemClass, name = "Count") %>%
  mutate(
    DominantGroup = factor(DominantGroup, levels = c("CON", "ELC", "ORG"))
  )

# 保存原始计数表
write.csv(
  fig4a_df,
  file.path(out_dir, "Figure4A_chemical_class_counts.csv"),
  row.names = FALSE
)

# 计算百分比
fig4a_prop_df <- fig4a_df %>%
  group_by(DominantGroup) %>%
  mutate(
    Total = sum(Count),
    Proportion = Count / Total
  ) %>%
  ungroup()

# 按总数排序 chemical class
chemclass_order <- fig4a_prop_df %>%
  group_by(ChemClass) %>%
  summarise(total_n = sum(Count), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  pull(ChemClass)

fig4a_prop_df <- fig4a_prop_df %>%
  mutate(
    ChemClass = factor(ChemClass, levels = chemclass_order)
  )

# 5. 配色
class_cols <- c(
  "Phenolics / flavonoids / glycosides" = "#E2AC7E",
  "Amino acid / peptide derivatives" = "#779A93",
  "Lipid-related / terpenoid-related metabolites" = "#E3E4D4",
  "Nucleosides / nucleotides / nitrogenous metabolites" = "#BD9095",
  "Other specialized / uncertain compounds" = "#EF8C39",
  "Organic acids / small polar metabolites" = "#A91C5E",
  "Tetrapyrroles / pigments" = "#2E4690",
  "Carbohydrates / polysaccharide-related" = "#BC5546",
  "Cofactors / vitamins" = "#72A4C5"
)

# 6. 作图
p_fig4a_prop <- ggplot(fig4a_prop_df, aes(x = DominantGroup, y = Proportion, fill = ChemClass)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proportion of curated differential metabolites",
    fill = "Chemical class"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_fig4a_prop



#==============================
# 9. class score calculation
#==============================

# 对每个代谢物在所有样本中做 z-score
expr_long2 <- expr_long %>%
  left_join(
    diff_pool2 %>% select(metab_id, DisplayName, ChemClass, ChemSubclass),
    by = "metab_id"
  ) %>%
  filter(!is.na(ChemClass), ChemClass != "") %>%
  group_by(metab_id) %>%
  mutate(
    z_abundance = as.numeric(scale(Abundance))
  ) %>%
  ungroup()

# 每个样本、每个类求平均 class score
class_score_df <- expr_long2 %>%
  group_by(Sample, Group, ChemClass) %>%
  summarise(
    class_score = mean(z_abundance, na.rm = TRUE),
    n_metabolites = n(),
    .groups = "drop"
  )

write.csv(
  class_score_df,
  file.path(out_dir, "Figure4B_class_scores_long.csv"),
  row.names = FALSE
)

#==============================
# 10. statistics for class scores
#==============================

overall_stats <- class_score_df %>%
  group_by(ChemClass) %>%
  summarise(
    overall_p = kruskal.test(class_score ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    overall_fdr = p.adjust(overall_p, method = "BH")
  )

pairwise_comparisons <- list(
  c("CON", "ELC"),
  c("CON", "ORG"),
  c("ELC", "ORG")
)

pairwise_stats_long <- class_score_df %>%
  group_by(ChemClass) %>%
  group_modify(~{
    df_now <- .x
    
    res <- lapply(pairwise_comparisons, function(comp) {
      g1 <- comp[1]
      g2 <- comp[2]
      
      x1 <- df_now$class_score[df_now$Group == g1]
      x2 <- df_now$class_score[df_now$Group == g2]
      
      wt <- wilcox.test(x1, x2, exact = FALSE)
      
      data.frame(
        comparison = paste0(g1, "_vs_", g2),
        group1 = g1,
        group2 = g2,
        pairwise_p = wt$p.value
      )
    })
    
    bind_rows(res) %>%
      mutate(pairwise_fdr = p.adjust(pairwise_p, method = "BH"))
  }) %>%
  ungroup()

pairwise_stats_wide <- pairwise_stats_long %>%
  pivot_wider(
    names_from = comparison,
    values_from = c(pairwise_p, pairwise_fdr)
  )

class_score_stats <- overall_stats %>%
  left_join(pairwise_stats_wide, by = "ChemClass") %>%
  arrange(overall_p)

write.csv(
  class_score_stats,
  file.path(out_dir, "Figure4B_class_score_stats.csv"),
  row.names = FALSE
)

#==============================
# 11. summary data for Figure 4B
#==============================

class_score_plot_df <- class_score_df %>%
  group_by(ChemClass, Group) %>%
  summarise(
    mean_score = mean(class_score, na.rm = TRUE),
    sd_score   = sd(class_score, na.rm = TRUE),
    se_score   = sd_score / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# 主图只保留这 5 类
main_classes <- c(
  "Phenolics / flavonoids / glycosides",
  "Amino acid / peptide derivatives",
  "Cofactors / vitamins",
  "Nucleosides / nucleotides / nitrogenous metabolites",
  "Organic acids / small polar metabolites"
)

class_score_plot_df2 <- class_score_plot_df %>%
  filter(ChemClass %in% main_classes)

overall_stats2 <- overall_stats %>%
  filter(ChemClass %in% main_classes) %>%
  mutate(
    label = ifelse(
      overall_p < 0.001,
      "P < 0.001",
      paste0("P = ", sprintf("%.3f", overall_p))
    )
  )

# 顺序：正文里想讲的顺序
chem_order <- rev(main_classes)

class_score_plot_df2 <- class_score_plot_df2 %>%
  mutate(
    ChemClass = factor(ChemClass, levels = chem_order),
    Group = factor(Group, levels = c("CON", "ELC", "ORG"))
  )

overall_stats2 <- overall_stats2 %>%
  mutate(
    ChemClass = factor(ChemClass, levels = chem_order)
  )

group_cols <- c(
  "CON" = "#E64B35",
  "ELC" = "#00A087",
  "ORG" = "#4DBBD5"
)

# P值标签位置
xmax <- max(class_score_plot_df2$mean_score + class_score_plot_df2$se_score, na.rm = TRUE)
xmin <- min(class_score_plot_df2$mean_score - class_score_plot_df2$se_score, na.rm = TRUE)
label_x <- xmax + 0.28 * (xmax - xmin + 0.1)

#==============================
# 12. plot Figure 4B
#==============================

p_fig4b <- ggplot(class_score_plot_df2,
                  aes(y = ChemClass, x = mean_score, fill = Group)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.62,
    color = "black",
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(
      xmin = mean_score - se_score,
      xmax = mean_score + se_score,
      color = Group
    ),
    position = position_dodge(width = 0.75),
    width = 0.18,
    linewidth = 0.55
  ) +
  geom_text(
    data = overall_stats2,
    aes(y = ChemClass, x = label_x, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.3
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  coord_cartesian(xlim = c(xmin - 0.1, label_x + 0.25)) +
  labs(
    x = "Class score (mean ± SE)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )

p_fig4b


#==============================
# 15. optional summary table for text
#==============================
class_score_summary <- class_score_df %>%
  group_by(ChemClass, Group) %>%
  summarise(
    mean_class_score = mean(class_score, na.rm = TRUE),
    sd_class_score = sd(class_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  class_score_summary,
  file.path(out_dir, "Figure4B_class_score_summary_by_group.csv"),
  row.names = FALSE
)

cat("\nDone.\n")
cat("Main outputs:\n")
cat("1) Figure4A_chemical_class_composition.pdf / .png\n")
cat("2) Figure4B_class_scores.pdf / .png\n")
cat("3) Figure4B_class_score_stats.csv\n")
cat("4) Result4_diff_pool_with_classes_and_group_means.csv\n")

