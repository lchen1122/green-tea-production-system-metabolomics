# ============================================================
# 05_result3_append_refined_logic.R
#
# Append this block after 05_result3_serum_metabolome_redesign.R
#
# Refined Result 3 logic:
# Hormone replacement reveals sex-divergent modes of metabolic
# reorganization, rather than a universal difference in
# "metabolic reversibility".
#
# Fig3C: biologically curated representative metabolites
# Fig3D: biological-module trajectories
# Fig3E: directional pathway-level support
# ============================================================

required_objects <- c(
  "female_patterns","male_patterns","serum_mat","serum_meta",
  "group_mean_z","pathway_all","PLOT_MAIN_DIR",
  "TABLE_PATTERN_DIR","TABLE_PATHWAY_DIR","TABLE_CANDIDATE_DIR"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Run the main Result 3 script first. Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ============================================================
# 1. Fig3C — biologically curated representative metabolites
# ============================================================

preferred_metabolites <- tribble(
  ~BiologicalAxis, ~PreferredName, ~PreferredSex,
  "Steroid-related metabolism", "Cholesterol", "Female",
  "Steroid-related metabolism", "Cortisol", "Male",
  "Steroid-related metabolism", "Dihydrocortisol", "Male",
  "Arachidonic-acid lipid mediators", "5-HETE", "Male",
  "Arachidonic-acid lipid mediators", "20-HETE", "Male",
  "Bile-acid-related metabolism", "Lithocholate 3-O-Glucuronide", "Male",
  "Amino-acid / immune metabolism", "Kynurenine", "Female",
  "Amino-acid / immune metabolism", "L-Lysine", "Female",
  "Amino-acid / immune metabolism", "Pipecolic Acid", "Female",
  "Microbiome-associated co-metabolism", "Trimethylamine N-Oxide", "Female",
  "Sphingolipid metabolism", "Sphinganine", "Female"
)

normalize_name <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "")
}

candidate_pool <- bind_rows(
  female_patterns %>% mutate(Sex = "Female"),
  male_patterns %>% mutate(Sex = "Male")
) %>%
  mutate(NameKey = normalize_name(DisplayName))

preferred_lookup <- preferred_metabolites %>%
  mutate(NameKey = normalize_name(PreferredName))

curated_hits <- preferred_lookup %>%
  left_join(
    candidate_pool,
    by = c("NameKey", "PreferredSex" = "Sex")
  ) %>%
  filter(!is.na(metab_id)) %>%
  mutate(
    BestFDR = pmin(
      Loss_FDR, Replacement_FDR, Residual_FDR,
      na.rm = TRUE
    ),
    MaxAbsEffect = pmax(
      abs(Loss_effect),
      abs(Replacement_effect),
      abs(Residual_effect),
      na.rm = TRUE
    )
  ) %>%
  group_by(BiologicalAxis, PreferredName, PreferredSex) %>%
  arrange(BestFDR, desc(MaxAbsEffect), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

write.csv(
  curated_hits,
  file.path(
    TABLE_CANDIDATE_DIR,
    "serum_biologically_curated_metabolites_refined.csv"
  ),
  row.names = FALSE
)

make_metabolite_trajectory <- function(
    metab_id, display_name, sex_label, pattern_label, biological_axis
) {
  
  groups_use <- if (sex_label == "Female") {
    c("Sham_F","Ovx_F","OvxE_F")
  } else {
    c("Sham_M","Orx_M","OrxT_M")
  }
  
  stage_labels <- c("Intact","Gonadectomy","Replacement")
  
  df <- tibble(
    SampleID = colnames(serum_mat),
    Intensity = as.numeric(
      serum_mat[metab_id, , drop = TRUE]
    )
  ) %>%
    left_join(
      serum_meta %>% select(SampleID, Group),
      by = "SampleID"
    ) %>%
    filter(as.character(Group) %in% groups_use) %>%
    mutate(
      Group = factor(as.character(Group), levels = groups_use),
      Stage = factor(
        stage_labels[match(as.character(Group), groups_use)],
        levels = stage_labels
      )
    )
  
  ggplot(df, aes(x = Stage, y = Intensity)) +
    geom_jitter(width = 0.07, size = 1.2, alpha = 0.42) +
    stat_summary(
      aes(group = 1),
      fun = mean,
      geom = "line",
      linewidth = 0.85
    ) +
    stat_summary(fun = mean, geom = "point", size = 2.7) +
    labs(
      x = NULL,
      y = "Intensity",
      title = display_name,
      subtitle = paste0(
        biological_axis, " | ", sex_label, " | ", pattern_label
      )
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "italic"),
      plot.subtitle = element_text(size = 7.6)
    )
}

axis_priority <- c(
  "Steroid-related metabolism",
  "Arachidonic-acid lipid mediators",
  "Bile-acid-related metabolism",
  "Amino-acid / immune metabolism",
  "Microbiome-associated co-metabolism",
  "Sphingolipid metabolism"
)

curated_plot_set <- curated_hits %>%
  mutate(
    BiologicalAxis = factor(BiologicalAxis, levels = axis_priority)
  ) %>%
  arrange(BiologicalAxis, PreferredSex, BestFDR) %>%
  group_by(BiologicalAxis) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  slice_head(n = 8)

if (nrow(curated_plot_set) > 0) {
  curated_plots <- pmap(
    curated_plot_set,
    function(
    BiologicalAxis, PreferredName, PreferredSex, NameKey,
    ..., metab_id, DisplayName, Pattern
    ) {
      make_metabolite_trajectory(
        metab_id = metab_id,
        display_name = DisplayName,
        sex_label = PreferredSex,
        pattern_label = Pattern,
        biological_axis = as.character(BiologicalAxis)
      )
    }
  )
  
  p_fig3c <- wrap_plots(curated_plots, ncol = 4)
  
  ggsave(
    file.path(
      PLOT_MAIN_DIR,
      "Fig3C_biologically_curated_metabolite_trajectories.pdf"
    ),
    p_fig3c,
    width = 13.2,
    height = 7.2
  )
}

# ============================================================
# 2. Fig3D — biological-module trajectories
#
# NOTE:
# Curated summary only; not a formal enrichment result.
# ============================================================

module_members <- curated_hits %>%
  select(
    BiologicalAxis,
    metab_id,
    PreferredSex,
    DisplayName,
    Pattern
  ) %>%
  distinct()

make_module_trajectory <- function(module_name, sex_label) {
  
  groups_use <- if (sex_label == "Female") {
    c("Sham_F","Ovx_F","OvxE_F")
  } else {
    c("Sham_M","Orx_M","OrxT_M")
  }
  
  ids <- module_members %>%
    filter(
      BiologicalAxis == module_name,
      PreferredSex == sex_label
    ) %>%
    pull(metab_id) %>%
    unique()
  
  ids <- intersect(ids, rownames(group_mean_z))
  
  if (length(ids) == 0) return(tibble())
  
  tibble(
    BiologicalAxis = module_name,
    Sex = sex_label,
    Stage = factor(
      c("Intact","Gonadectomy","Replacement"),
      levels = c("Intact","Gonadectomy","Replacement")
    ),
    MeanZ = colMeans(
      group_mean_z[ids, groups_use, drop = FALSE],
      na.rm = TRUE
    ),
    N_metabolites = length(ids)
  )
}

module_names <- unique(
  as.character(module_members$BiologicalAxis)
)

module_trajectory <- map_dfr(
  module_names,
  ~ bind_rows(
    make_module_trajectory(.x, "Female"),
    make_module_trajectory(.x, "Male")
  )
)

write.csv(
  module_trajectory,
  file.path(
    TABLE_CANDIDATE_DIR,
    "serum_biological_module_trajectories.csv"
  ),
  row.names = FALSE
)

if (nrow(module_trajectory) > 0) {
  p_fig3d <- ggplot(
    module_trajectory,
    aes(
      x = Stage,
      y = MeanZ,
      group = Sex,
      linetype = Sex
    )
  ) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 0.35) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.4) +
    facet_wrap(~ BiologicalAxis, scales = "free_y") +
    labs(
      x = NULL,
      y = "Mean metabolite z-score",
      linetype = NULL,
      title = "Biological modules show sex-divergent modes of metabolic reorganization",
      subtitle = "Curated module trajectories summarize biologically related endogenous metabolites"
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(
      PLOT_MAIN_DIR,
      "Fig3D_biological_module_trajectories.pdf"
    ),
    p_fig3d,
    width = 10,
    height = 7
  )
}

# ============================================================
# 3. Fig3E — directional pathway-level support
# ============================================================

BIOLOGICAL_PATHWAYS <- c(
  "Steroid hormone biosynthesis",
  "Arachidonic acid metabolism",
  "Glycerophospholipid metabolism",
  "Linoleic acid metabolism",
  "alpha-Linolenic acid metabolism",
  "Bile secretion",
  "Tryptophan metabolism",
  "Arginine and proline metabolism",
  "D-Amino acid metabolism",
  "Phenylalanine metabolism"
)

PATHWAY_CONTRASTS <- c(
  "Sex_intact",
  "Female_loss",
  "Female_replacement",
  "Male_loss",
  "Male_replacement"
)

PATHWAY_LABELS <- c(
  Sex_intact = "Sex difference\nIntact",
  Female_loss = "Female\nGonadectomy",
  Female_replacement = "Female\nReplacement",
  Male_loss = "Male\nGonadectomy",
  Male_replacement = "Male\nReplacement"
)

pathway_bio <- pathway_all %>%
  mutate(
    PathwayName = if_else(
      !is.na(description) & description != "",
      description,
      pathway_id
    )
  ) %>%
  filter(
    Contrast %in% PATHWAY_CONTRASTS,
    PathwayName %in% BIOLOGICAL_PATHWAYS
  ) %>%
  complete(
    Contrast = PATHWAY_CONTRASTS,
    PathwayName = BIOLOGICAL_PATHWAYS,
    fill = list(
      HitN = 0,
      UpN = 0,
      DownN = 0,
      MeanEffect = 0,
      DirectionScore = 0,
      FDR_enrichment = NA_real_
    )
  ) %>%
  mutate(
    Contrast = factor(Contrast, levels = PATHWAY_CONTRASTS),
    PathwayName = factor(
      PathwayName,
      levels = rev(BIOLOGICAL_PATHWAYS)
    )
  )

write.csv(
  pathway_bio,
  file.path(
    TABLE_PATHWAY_DIR,
    "serum_biologically_curated_pathways_refined.csv"
  ),
  row.names = FALSE
)

p_fig3e <- ggplot(
  pathway_bio,
  aes(x = MeanEffect, y = PathwayName)
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
  geom_segment(
    aes(
      x = 0,
      xend = MeanEffect,
      yend = PathwayName
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(size = HitN),
    shape = 21,
    stroke = 0.35
  ) +
  facet_wrap(
    ~ Contrast,
    nrow = 1,
    labeller = as_labeller(PATHWAY_LABELS)
  ) +
  labs(
    x = "Mean effect of significant mapped metabolites",
    y = NULL,
    size = "Mapped hits",
    title = "Directional pathway-level support for sex-divergent metabolic remodeling",
    subtitle = "Positive and negative values indicate the average direction of significant mapped metabolites"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    PLOT_MAIN_DIR,
    "Fig3E_directional_pathway_support.pdf"
  ),
  p_fig3e,
  width = 13,
  height = 6.5
)

# ============================================================
# 4. Compact response-mode summary for writing
# ============================================================

response_mode_summary <- bind_rows(
  female_patterns %>% mutate(Sex = "Female"),
  male_patterns %>% mutate(Sex = "Male")
) %>%
  count(
    Sex,
    Pattern,
    name = "N"
  ) %>%
  group_by(Sex) %>%
  mutate(
    Percent = N / sum(N) * 100
  ) %>%
  ungroup() %>%
  mutate(
    Mode = case_when(
      Pattern %in% c(
        "Strong_toward_intact",
        "Partial_toward_intact"
      ) ~ "Toward-intact",
      Pattern == "Replacement_specific" ~ "Replacement-responsive",
      Pattern == "Persistent_after_replacement" ~ "Persistent",
      Pattern == "Complex_other" ~ "Complex",
      Pattern == "Loss_only" ~ "Loss-only",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(
    Sex,
    Mode
  ) %>%
  summarise(
    N = sum(N),
    Percent = sum(Percent),
    .groups = "drop"
  )

write.csv(
  response_mode_summary,
  file.path(
    TABLE_PATTERN_DIR,
    "serum_response_modes_compact.csv"
  ),
  row.names = FALSE
)

cat(
  "\n\n========================================\n",
  "REFINED RESULT 3 LOGIC COMPLETED\n",
  "========================================\n",
  "New figures:\n",
  "- Fig3C_biologically_curated_metabolite_trajectories.pdf\n",
  "- Fig3D_biological_module_trajectories.pdf\n",
  "- Fig3E_directional_pathway_support.pdf\n\n",
  "Interpretation:\n",
  "Hormone replacement reveals sex-divergent modes of metabolic reorganization.\n",
  "Male responses contain more toward-intact and replacement-responsive patterns,\n",
  "whereas female responses contain a larger persistent/complex component.\n",
  "Do NOT generalize this as intrinsic male metabolic reversibility.\n",
  sep = ""
)

