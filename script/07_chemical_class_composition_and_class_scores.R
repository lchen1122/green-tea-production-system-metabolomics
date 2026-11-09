############################################################
# 07. Chemical class composition and class-score analysis
# Tea metabolomics: CON / ELC / ORG
############################################################

rm(list = ls())
gc()

# ============================================================
# 0. Packages
# ============================================================
library(tidyverse)
library(limma)

# ============================================================
# 1. Paths
# ============================================================
data_file <- paste0(
  "F:/科研/有机生态低碳茶叶分析/",
  "代谢组结果/1.Preprocess/mix/all_metab_data.txt"
)

annotation_file <- paste0(
  "F:/科研/有机生态低碳茶叶分析/",
  "Paper/Result3/annotation_screened_full.csv"
)

# Original manually curated chemical-class file.
class_file <- paste0(
  "F:/科研/有机生态低碳茶叶分析/",
  "Paper/Result4/chemical_class_annotation_template.csv"
)

out_dir <- paste0(
  "F:/科研/有机生态低碳茶叶分析/",
  "Paper/Revised2/Result4/"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Parameters
# ============================================================
# Metabolomics abundance data are already log2-transformed upstream.
logfc_cutoff <- 0.58
p_cutoff <- 0.05
fdr_cutoff <- 0.10

group_levels <- c(
  "CON",
  "ELC",
  "ORG"
)

comparisons <- list(
  c("CON", "ELC"),
  c("CON", "ORG"),
  c("ELC", "ORG")
)

# ============================================================
# 3. Read data
# ============================================================
metab_df <- read.delim(
  data_file,
  check.names = FALSE
)

anno_df <- read.csv(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  distinct(
    metab_id,
    .keep_all = TRUE
  )

class_raw <- read.csv(
  class_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  "metab_id" %in% colnames(metab_df),
  all(
    c(
      "metab_id",
      "Metabolite"
    ) %in% colnames(anno_df)
  ),
  all(
    c(
      "metab_id",
      "ChemClass"
    ) %in% colnames(class_raw)
  )
)

# ============================================================
# 4. Sample information and expression matrix
# ============================================================
sample_cols <- grep(
  "^(CON|ELC|ORG)_",
  colnames(metab_df),
  value = TRUE
)

sample_info <- tibble(
  Sample = sample_cols,
  Group = case_when(
    grepl("^CON_", Sample) ~ "CON",
    grepl("^ELC_", Sample) ~ "ELC",
    grepl("^ORG_", Sample) ~ "ORG",
    TRUE ~ NA_character_
  )
) %>%
  filter(
    !is.na(Group)
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = group_levels
    )
  )

expr_mat <- as.matrix(
  metab_df[
    ,
    sample_info$Sample,
    drop = FALSE
  ]
)

mode(expr_mat) <- "numeric"
rownames(expr_mat) <- metab_df$metab_id

# ============================================================
# 5. Annotation and chemical-class tables
# ============================================================
anno_keep_cols <- intersect(
  c(
    "metab_id",
    "Metabolite",
    "DisplayName",
    "annotation_screen",
    "Formula",
    "Adducts",
    "compound_source",
    "Score",
    "level"
  ),
  colnames(anno_df)
)

anno_use <- anno_df %>%
  select(
    all_of(anno_keep_cols)
  )

# Normalize the original manual class table.
if (!"ChemSubclass" %in% colnames(class_raw)) {
  class_raw$ChemSubclass <- NA_character_
}

if (!"Note" %in% colnames(class_raw)) {
  class_raw$Note <- NA_character_
}

class_diff <- class_raw %>%
  transmute(
    metab_id = as.character(metab_id),
    ChemClass = na_if(
      trimws(
        as.character(ChemClass)
      ),
      ""
    ),
    ChemSubclass = na_if(
      trimws(
        as.character(ChemSubclass)
      ),
      ""
    ),
    Note = na_if(
      trimws(
        as.character(Note)
      ),
      ""
    )
  ) %>%
  distinct(
    metab_id,
    .keep_all = TRUE
  ) %>%
  filter(
    !is.na(ChemClass)
  )

cat(
  "\nManual chemical-class annotations available:",
  nrow(class_diff),
  "\n"
)

# ============================================================
# 6. Differential-metabolite analysis used for Figure 4A/4B
# ============================================================
run_limma_pair <- function(
  expr_mat,
  sample_info,
  g1,
  g2,
  anno_use,
  logfc_cutoff,
  p_cutoff,
  fdr_cutoff
) {

  info_sub <- sample_info %>%
    filter(
      Group %in% c(
        g1,
        g2
      )
    ) %>%
    mutate(
      Group = factor(
        Group,
        levels = c(
          g1,
          g2
        )
      )
    )

  expr_sub <- expr_mat[
    ,
    info_sub$Sample,
    drop = FALSE
  ]

  design <- model.matrix(
    ~ 0 + Group,
    data = info_sub
  )

  colnames(design) <- levels(
    info_sub$Group
  )

  fit <- lmFit(
    expr_sub,
    design
  )

  contrast_mat <- makeContrasts(
    contrasts = paste0(
      g1,
      "-",
      g2
    ),
    levels = design
  )

  fit2 <- contrasts.fit(
    fit,
    contrast_mat
  )

  fit2 <- eBayes(
    fit2
  )

  tt <- topTable(
    fit2,
    number = Inf,
    sort.by = "none"
  ) %>%
    rownames_to_column(
      "metab_id"
    )

  res <- tt %>%
    left_join(
      anno_use,
      by = "metab_id"
    ) %>%
    mutate(
      Comparison = paste0(
        g1,
        "_vs_",
        g2
      ),
      higher_in = case_when(
        logFC > 0 ~ g1,
        logFC < 0 ~ g2,
        TRUE ~ "Tie"
      ),
      abs_logFC = abs(
        logFC
      ),
      sig_rawp = (
        P.Value < p_cutoff &
          abs_logFC >= logfc_cutoff
      ),
      sig_fdr = (
        adj.P.Val < fdr_cutoff &
          abs_logFC >= logfc_cutoff
      ),
      curated_annotation = (
        !is.na(Metabolite) &
          Metabolite != "" &
          !grepl(
            "^metab_",
            Metabolite,
            ignore.case = TRUE
          )
      )
    )

  if ("annotation_screen" %in% colnames(res)) {
    res <- res %>%
      mutate(
        curated_annotation = curated_annotation &
          !is.na(annotation_screen) &
          annotation_screen %in% c(
            "recommended",
            "keep_with_caution"
          )
      )
  }

  res %>%
    mutate(
      curated_sig_rawp = sig_rawp &
        curated_annotation,
      curated_sig_fdr = sig_fdr &
        curated_annotation
    )
}

res_list <- list()

for (cmp in comparisons) {

  g1 <- cmp[1]
  g2 <- cmp[2]
  cmp_name <- paste0(
    g1,
    "_vs_",
    g2
  )

  res_list[[cmp_name]] <- run_limma_pair(
    expr_mat = expr_mat,
    sample_info = sample_info,
    g1 = g1,
    g2 = g2,
    anno_use = anno_use,
    logfc_cutoff = logfc_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff
  )
}

# Original exploratory differential-metabolite definition:
# raw P < 0.05 and |log2FC| >= 0.58 among curated annotations.
diff_union_long <- bind_rows(
  res_list
) %>%
  filter(
    curated_sig_rawp
  ) %>%
  distinct(
    Comparison,
    metab_id,
    .keep_all = TRUE
  )

# Unique curated differential metabolites with manual class annotation.
diff_union <- diff_union_long %>%
  distinct(
    metab_id,
    .keep_all = TRUE
  ) %>%
  left_join(
    class_diff,
    by = "metab_id"
  ) %>%
  filter(
    !is.na(ChemClass)
  )

cat(
  "\nDifferential-metabolite pool:\n",
  "Unique curated differential metabolites: ",
  n_distinct(
    diff_union_long$metab_id
  ),
  "\nWith ChemClass annotation: ",
  nrow(
    diff_union
  ),
  "\n",
  sep = ""
)

if (nrow(diff_union) == 0) {
  stop(
    paste0(
      "No curated differential metabolite matched the manual ",
      "chemical-class annotation table."
    )
  )
}

write.csv(
  diff_union,
  file.path(
    out_dir,
    "curated_differential_metabolites_with_classes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. Figure 4A
#    Chemical-class composition of curated differential metabolites
# ============================================================

# Group mean abundance for each differential metabolite.
group_mean_long <- metab_df %>%
  filter(
    metab_id %in% diff_union$metab_id
  ) %>%
  select(
    metab_id,
    all_of(
      sample_info$Sample
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      sample_info$Sample
    ),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(
    sample_info,
    by = "Sample"
  ) %>%
  group_by(
    metab_id,
    Group
  ) %>%
  summarise(
    mean_abundance = mean(
      Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Dominant group = group with the highest mean abundance.
dominant_group_df <- group_mean_long %>%
  filter(
    is.finite(
      mean_abundance
    )
  ) %>%
  group_by(
    metab_id
  ) %>%
  slice_max(
    order_by = mean_abundance,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    metab_id,
    DominantGroup = as.character(
      Group
    )
  )

fig4a_pool <- diff_union %>%
  left_join(
    dominant_group_df,
    by = "metab_id"
  ) %>%
  filter(
    !is.na(
      DominantGroup
    )
  )

fig4a_df <- fig4a_pool %>%
  count(
    DominantGroup,
    ChemClass,
    name = "Count"
  ) %>%
  mutate(
    DominantGroup = factor(
      DominantGroup,
      levels = group_levels
    )
  )

fig4a_prop_df <- fig4a_df %>%
  group_by(
    DominantGroup
  ) %>%
  mutate(
    Total = sum(
      Count
    ),
    Proportion = Count /
      Total
  ) %>%
  ungroup()

write.csv(
  fig4a_prop_df,
  file.path(
    out_dir,
    "Figure4A_chemical_class_counts_and_proportions.csv"
  ),
  row.names = FALSE
)

chemclass_order <- fig4a_prop_df %>%
  group_by(
    ChemClass
  ) %>%
  summarise(
    total_n = sum(
      Count
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      total_n
    )
  ) %>%
  pull(
    ChemClass
  )

fig4a_prop_df <- fig4a_prop_df %>%
  mutate(
    ChemClass = factor(
      ChemClass,
      levels = chemclass_order
    )
  )

class_cols <- c(
  "Phenolics / flavonoids / glycosides" = "#E2AC7E",
  "Amino acid / peptide derivatives" = "#779A93",
  "Lipid-related / terpenoid-related metabolites" = "#E3E4D4",
  "Nucleosides / nucleotides / nitrogenous metabolites" = "#BD9095",
  "Other specialized / uncertain compounds" = "#EF8C39",
  "Organic acids / small polar metabolites" = "#A91C5E",
  "Tetrapyrroles / pigments" = "#2E4690",
  "Carbohydrates / polysaccharide-related" = "#BC5546",
  "Steroids / hormone-related metabolites" = "#A775B4",
  "Cofactors / vitamins" = "#72A4C5"
)

missing_palette <- setdiff(
  unique(
    as.character(
      fig4a_prop_df$ChemClass
    )
  ),
  names(
    class_cols
  )
)

if (length(missing_palette) > 0) {
  stop(
    paste0(
      "Figure 4A contains ChemClass labels missing from class_cols: ",
      paste(
        missing_palette,
        collapse = "; "
      )
    )
  )
}

p_fig4a <- ggplot(
  fig4a_prop_df,
  aes(
    x = DominantGroup,
    y = Proportion,
    fill = ChemClass
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = Count
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 3
  ) +
  scale_fill_manual(
    values = class_cols
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion of curated differential metabolites",
    fill = "Chemical class"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      face = "bold"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

print(
  p_fig4a
)

ggsave(
  file.path(
    out_dir,
    "Figure4A_chemical_class_composition.pdf"
  ),
  p_fig4a,
  width = 8.2,
  height = 5.5
)

ggsave(
  file.path(
    out_dir,
    "Figure4A_chemical_class_composition.png"
  ),
  p_fig4a,
  width = 8.2,
  height = 5.5,
  dpi = 600
)

# ============================================================
# 8. Figure 4B
#    Class-level quantitative summary of the curated
#    differential-metabolite set
# ============================================================

# Long-format abundance table for the SAME metabolites used above.
expr_long_class <- metab_df %>%
  filter(
    metab_id %in% diff_union$metab_id
  ) %>%
  select(
    metab_id,
    all_of(
      sample_info$Sample
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      sample_info$Sample
    ),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(
    sample_info,
    by = "Sample"
  ) %>%
  left_join(
    diff_union %>%
      select(
        metab_id,
        ChemClass
      ) %>%
      distinct(),
    by = "metab_id"
  ) %>%
  filter(
    !is.na(
      ChemClass
    )
  )

# Standardize each metabolite across all biological samples.
expr_long_class_z <- expr_long_class %>%
  group_by(
    metab_id
  ) %>%
  mutate(
    z_abundance = as.numeric(
      scale(
        Abundance
      )
    )
  ) %>%
  ungroup() %>%
  filter(
    is.finite(
      z_abundance
    )
  )

# Number of metabolites contributing to each class.
class_n_df <- expr_long_class_z %>%
  distinct(
    metab_id,
    ChemClass
  ) %>%
  count(
    ChemClass,
    name = "n_metabolites"
  ) %>%
  arrange(
    desc(
      n_metabolites
    )
  )

write.csv(
  class_n_df,
  file.path(
    out_dir,
    "Figure4B_metabolite_counts_by_class.csv"
  ),
  row.names = FALSE
)

# Sample-level class score:
# mean standardized abundance of all metabolites in the class.
class_score_df <- expr_long_class_z %>%
  group_by(
    Sample,
    Group,
    ChemClass
  ) %>%
  summarise(
    class_score = mean(
      z_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  class_score_df,
  file.path(
    out_dir,
    "Figure4B_class_scores_long.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8A. Overall statistics
#     Test ALL represented classes, then BH-correct across ALL
#     tested classes.
# ------------------------------------------------------------
overall_stats <- class_score_df %>%
  group_by(
    ChemClass
  ) %>%
  summarise(
    P_raw = kruskal.test(
      class_score ~ Group
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    P_BH = p.adjust(
      P_raw,
      method = "BH"
    )
  ) %>%
  left_join(
    class_n_df,
    by = "ChemClass"
  ) %>%
  arrange(
    P_BH,
    P_raw
  )

write.csv(
  overall_stats,
  file.path(
    out_dir,
    "Figure4B_overall_statistics_all_classes.csv"
  ),
  row.names = FALSE
)

print(
  overall_stats
)

sig_classes <- overall_stats %>%
  filter(
    P_BH < 0.1
  ) %>%
  pull(
    ChemClass
  )

write.csv(
  overall_stats,
  file.path(
    out_dir,
    "Figure4B_significant_classes_BH.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFigure 4B statistics:\n",
  "Tested chemical classes: ",
  nrow(
    overall_stats
  ),
  "\nBH-significant classes: ",
  length(
    sig_classes
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8B. Pairwise comparisons
#     Secondary interpretation of significant classes only.
# ------------------------------------------------------------
pairwise_stats <- tibble()

if (length(sig_classes) > 0) {

  pairwise_stats <- class_score_df %>%
    filter(
      ChemClass %in% sig_classes
    ) %>%
    group_by(
      ChemClass
    ) %>%
    group_modify(
      ~ {

        df_now <- .x

        bind_rows(
          lapply(
            comparisons,
            function(comp) {

              g1 <- comp[1]
              g2 <- comp[2]

              wt <- wilcox.test(
                df_now$class_score[
                  df_now$Group == g1
                ],
                df_now$class_score[
                  df_now$Group == g2
                ],
                exact = FALSE
              )

              tibble(
                comparison = paste0(
                  g1,
                  "_vs_",
                  g2
                ),
                P_raw = wt$p.value
              )
            }
          )
        ) %>%
          mutate(
            P_BH_within_class = p.adjust(
              P_raw,
              method = "BH"
            )
          )
      }
    ) %>%
    ungroup()

  write.csv(
    pairwise_stats,
    file.path(
      out_dir,
      "Figure4B_pairwise_statistics_significant_classes.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 8C. Figure 4B
# ------------------------------------------------------------
if (length(sig_classes) > 0) {

  class_score_plot_df <- class_score_df %>%
    filter(
      ChemClass %in% sig_classes
    ) %>%
    group_by(
      ChemClass,
      Group
    ) %>%
    summarise(
      mean_score = mean(
        class_score,
        na.rm = TRUE
      ),
      sd_score = sd(
        class_score,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  class_label_df <- overall_stats %>%
    filter(
      ChemClass %in% sig_classes
    ) %>%
    mutate(
      ChemClass_label = paste0(
        ChemClass,
        " (n=",
        n_metabolites,
        ")"
      ),
      P_label = case_when(
        P_BH < 0.1 ~ "BH-adjusted P < 0.001",
        TRUE ~ paste0(
          "BH-adjusted P = ",
          sprintf(
            "%.3f",
            P_BH
          )
        )
      )
    )

  class_score_plot_df <- class_score_plot_df %>%
    left_join(
      class_label_df %>%
        select(
          ChemClass,
          ChemClass_label,
          P_raw,
          P_BH,
          P_label
        ),
      by = "ChemClass"
    )

  # Keep a stable, result-driven order.
  chem_order <- class_label_df %>%
    arrange(
      desc(
        P_BH
      ),
      desc(
        P_raw
      )
    ) %>%
    pull(
      ChemClass_label
    )

  class_score_plot_df <- class_score_plot_df %>%
    mutate(
      ChemClass_label = factor(
        ChemClass_label,
        levels = chem_order
      ),
      Group = factor(
        Group,
        levels = group_levels
      )
    )

  stats_plot_df <- class_label_df %>%
    mutate(
      ChemClass_label = factor(
        ChemClass_label,
        levels = chem_order
      )
    )

  group_cols <- c(
    "CON" = "#E64B35",
    "ELC" = "#00A087",
    "ORG" = "#4DBBD5"
  )

  xmax <- max(
    class_score_plot_df$mean_score +
      class_score_plot_df$sd_score,
    na.rm = TRUE
  )

  xmin <- min(
    class_score_plot_df$mean_score -
      class_score_plot_df$sd_score,
    na.rm = TRUE
  )

  label_x <- xmax +
    0.25 * (
      xmax -
        xmin +
        0.1
    )

  p_fig4b <- ggplot(
    class_score_plot_df,
    aes(
      y = ChemClass_label,
      x = mean_score,
      fill = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      color = "grey50",
      linewidth = 0.5
    ) +
    geom_col(
      position = position_dodge(
        width = 0.75
      ),
      width = 0.62,
      color = "black",
      linewidth = 0.3
    ) +
    geom_errorbar(
      aes(
        xmin = mean_score -
          sd_score,
        xmax = mean_score +
          sd_score,
        color = Group
      ),
      position = position_dodge(
        width = 0.75
      ),
      width = 0.18,
      linewidth = 0.55
    ) +
    geom_text(
      data = stats_plot_df,
      aes(
        y = ChemClass_label,
        x = label_x,
        label = P_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.2
    ) +
    scale_fill_manual(
      values = group_cols
    ) +
    scale_color_manual(
      values = group_cols
    ) +
    coord_cartesian(
      xlim = c(
        xmin - 0.1,
        label_x + 0.4
      )
    ) +
    labs(
      x = "Class score (mean ± SD)",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.text = element_text(
        face = "bold"
      ),
      axis.text.y = element_text(
        face = "bold"
      ),
      axis.title.x = element_text(
        face = "bold"
      )
    )

  print(
    p_fig4b
  )

  ggsave(
    file.path(
      out_dir,
      "Figure4B_class_scores.pdf"
    ),
    p_fig4b,
    width = 9.6,
    height = 5.8
  )

  ggsave(
    file.path(
      out_dir,
      "Figure4B_class_scores.png"
    ),
    p_fig4b,
    width = 9.6,
    height = 5.8,
    dpi = 600
  )

} else {

  message(
    "\nNo chemical class remained significant after BH correction; Figure 4B was not generated."
  )
}

# ------------------------------------------------------------
# 8D. Complete class-score summary
# ------------------------------------------------------------
class_score_summary <- class_score_df %>%
  group_by(
    ChemClass,
    Group
  ) %>%
  summarise(
    mean_score = mean(
      class_score,
      na.rm = TRUE
    ),
    sd_score = sd(
      class_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    overall_stats,
    by = "ChemClass"
  ) %>%
  mutate(
    Shown_in_Figure4B = ChemClass %in%
      sig_classes
  ) %>%
  arrange(
    P_BH,
    P_raw,
    ChemClass,
    Group
  )

write.csv(
  class_score_summary,
  file.path(
    out_dir,
    "Figure4B_class_score_summary_all_classes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. Finish
# ============================================================
cat(
  "\n07 analysis completed.\n"
)

