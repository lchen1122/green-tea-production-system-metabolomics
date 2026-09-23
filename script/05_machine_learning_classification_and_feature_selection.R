# ============================================================
# APPEND TO END OF:
#   05_result3_serum_metabolome_redesign.R
#
# FINAL 5-PANEL RESULT 3 LAYOUT
#
# A. Fate of baseline sex-biased metabolites after gonadectomy
# B. Simplified female/male response architecture
# C. Biologically interpretable representative metabolites
# D. Chemical-class directional remodeling
# E. Directional KEGG pathway remodeling
#
# This block assumes the main Result 3 script has already created:
#   female_patterns
#   male_patterns
#   serum_mat
#   serum_meta
#   annotation_use
#   pathway_all
#   PLOT_MAIN_DIR
#   PLOT_SUPP_DIR
#   TABLE_PATTERN_DIR
#   TABLE_CANDIDATE_DIR
#   TABLE_PATHWAY_DIR
# ============================================================


# ============================================================
# 0. Checks
# ============================================================

required_objects <- c(
  "female_patterns",
  "male_patterns",
  "serum_mat",
  "serum_meta",
  "annotation_use",
  "pathway_all",
  "PLOT_MAIN_DIR",
  "PLOT_SUPP_DIR",
  "TABLE_PATTERN_DIR",
  "TABLE_CANDIDATE_DIR",
  "TABLE_PATHWAY_DIR"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Please run the main Result 3 script first. Missing objects: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


# ============================================================
# 1. Fig3B — simplified response architecture
#
# Collapse the detailed six response classes into three
# biologically interpretable modes:
#
#   1. Toward-intact / replacement-responsive
#   2. Persistent / complex
#   3. Loss-only
#
# Keep the full six-class plot in supplementary.
# ============================================================

response_mode_summary <- bind_rows(
  female_patterns %>%
    mutate(
      Sex = "Female"
    ),
  male_patterns %>%
    mutate(
      Sex = "Male"
    )
) %>%
  filter(
    Pattern != "No_significant_change"
  ) %>%
  mutate(
    Mode = case_when(
      Pattern %in%
        c(
          "Strong_toward_intact",
          "Partial_toward_intact",
          "Replacement_specific"
        ) ~
        "Toward-intact / replacement-responsive",
      
      Pattern %in%
        c(
          "Persistent_after_replacement",
          "Complex_other"
        ) ~
        "Persistent / complex",
      
      Pattern ==
        "Loss_only" ~
        "Loss-only",
      
      TRUE ~
        "Other"
    )
  ) %>%
  count(
    Sex,
    Mode,
    name = "N"
  ) %>%
  group_by(
    Sex
  ) %>%
  mutate(
    Percent =
      N /
      sum(N) *
      100
  ) %>%
  ungroup()

mode_order <- c(
  "Toward-intact / replacement-responsive",
  "Persistent / complex",
  "Loss-only"
)

response_mode_summary <- response_mode_summary %>%
  mutate(
    Mode = factor(
      Mode,
      levels = mode_order
    )
  )

write.csv(
  response_mode_summary,
  file.path(
    TABLE_PATTERN_DIR,
    "serum_response_modes_simplified.csv"
  ),
  row.names = FALSE
)


p_fig3b_simple <- ggplot(
  response_mode_summary,
  aes(
    x = Sex,
    y = Percent,
    fill = Mode
  )
) +
  geom_col(
    width = 0.62
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      scale = 1
    )
  ) +
  labs(
    x = NULL,
    y = "Composition of responsive metabolites",
    fill = NULL,
    title = "Female and male serum metabolomes adopt distinct response modes"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "right"
  )


ggsave(
  file.path(
    PLOT_MAIN_DIR,
    "Fig3B_simplified_response_modes.pdf"
  ),
  p_fig3b_simple,
  width = 7.6,
  height = 5.2
)


# ============================================================
# 2. Fig3C — biologically curated metabolite trajectories
# ============================================================

preferred_metabolites <- tribble(
  ~BiologicalAxis, ~PreferredName, ~PreferredSex,
  
  "Steroid-related metabolism",
  "Cholesterol",
  "Female",
  
  "Steroid-related metabolism",
  "Cortisol",
  "Male",
  
  "Steroid-related metabolism",
  "Dihydrocortisol",
  "Male",
  
  "Arachidonic-acid lipid mediators",
  "5-HETE",
  "Male",
  
  "Arachidonic-acid lipid mediators",
  "20-HETE",
  "Male",
  
  "Bile-acid-related metabolism",
  "Lithocholate 3-O-Glucuronide",
  "Male",
  
  "Amino-acid / immune metabolism",
  "Kynurenine",
  "Female",
  
  "Amino-acid / immune metabolism",
  "L-Lysine",
  "Female",
  
  "Amino-acid / immune metabolism",
  "Pipecolic Acid",
  "Female",
  
  "Microbiome-associated co-metabolism",
  "Trimethylamine N-Oxide",
  "Female",
  
  "Sphingolipid metabolism",
  "Sphinganine",
  "Female"
)


normalize_name <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all(
      "[^a-z0-9]+",
      ""
    )
}


candidate_pool <- bind_rows(
  female_patterns %>%
    mutate(
      Sex = "Female"
    ),
  
  male_patterns %>%
    mutate(
      Sex = "Male"
    )
) %>%
  mutate(
    NameKey =
      normalize_name(
        DisplayName
      )
  )


preferred_lookup <- preferred_metabolites %>%
  mutate(
    NameKey =
      normalize_name(
        PreferredName
      )
  )


curated_hits <- preferred_lookup %>%
  left_join(
    candidate_pool,
    by = c(
      "NameKey",
      "PreferredSex" = "Sex"
    )
  ) %>%
  filter(
    !is.na(
      metab_id
    )
  ) %>%
  mutate(
    BestFDR =
      pmin(
        Loss_FDR,
        Replacement_FDR,
        Residual_FDR,
        na.rm = TRUE
      ),
    
    MaxAbsEffect =
      pmax(
        abs(
          Loss_effect
        ),
        abs(
          Replacement_effect
        ),
        abs(
          Residual_effect
        ),
        na.rm = TRUE
      )
  ) %>%
  group_by(
    BiologicalAxis,
    PreferredName,
    PreferredSex
  ) %>%
  arrange(
    BestFDR,
    desc(
      MaxAbsEffect
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup()


axis_priority <- c(
  "Steroid-related metabolism",
  "Arachidonic-acid lipid mediators",
  "Bile-acid-related metabolism",
  "Amino-acid / immune metabolism",
  "Microbiome-associated co-metabolism",
  "Sphingolipid metabolism"
)


curated_hits <- curated_hits %>%
  mutate(
    BiologicalAxis = factor(
      BiologicalAxis,
      levels = axis_priority
    )
  ) %>%
  arrange(
    BiologicalAxis,
    PreferredSex,
    BestFDR
  )


write.csv(
  curated_hits,
  file.path(
    TABLE_CANDIDATE_DIR,
    "serum_biologically_curated_metabolites_final.csv"
  ),
  row.names = FALSE
)


make_metabolite_trajectory <- function(
    metab_id,
    display_name,
    sex_label,
    pattern_label,
    biological_axis
) {
  
  groups_use <- if (
    sex_label ==
    "Female"
  ) {
    c(
      "Sham_F",
      "Ovx_F",
      "OvxE_F"
    )
  } else {
    c(
      "Sham_M",
      "Orx_M",
      "OrxT_M"
    )
  }
  
  stage_labels <- c(
    "Intact",
    "Gonadectomy",
    "Replacement"
  )
  
  df <- tibble(
    SampleID =
      colnames(
        serum_mat
      ),
    Intensity =
      as.numeric(
        serum_mat[
          metab_id,
          ,
          drop = TRUE
        ]
      )
  ) %>%
    left_join(
      serum_meta %>%
        select(
          SampleID,
          Group
        ),
      by = "SampleID"
    ) %>%
    filter(
      as.character(
        Group
      ) %in%
        groups_use
    ) %>%
    mutate(
      Group = factor(
        as.character(
          Group
        ),
        levels = groups_use
      ),
      
      Stage = factor(
        stage_labels[
          match(
            as.character(
              Group
            ),
            groups_use
          )
        ],
        levels = stage_labels
      )
    )
  
  ggplot(
    df,
    aes(
      x = Stage,
      y = Intensity
    )
  ) +
    geom_jitter(
      width = 0.07,
      size = 1.2,
      alpha = 0.42
    ) +
    stat_summary(
      aes(
        group = 1
      ),
      fun = mean,
      geom = "line",
      linewidth = 0.85
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2.7
    ) +
    labs(
      x = NULL,
      y = "Intensity",
      title = display_name,
      subtitle = paste0(
        biological_axis,
        " | ",
        sex_label,
        " | ",
        pattern_label
      )
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        face = "italic"
      ),
      plot.subtitle = element_text(
        size = 7.5
      )
    )
}


curated_plot_set <- curated_hits %>%
  group_by(
    BiologicalAxis
  ) %>%
  slice_head(
    n = 2
  ) %>%
  ungroup() %>%
  slice_head(
    n = 8
  )


if (
  nrow(
    curated_plot_set
  ) >
  0
) {
  
  curated_plots <- pmap(
    curated_plot_set,
    function(
    BiologicalAxis,
    PreferredName,
    PreferredSex,
    NameKey,
    ...,
    metab_id,
    DisplayName,
    Pattern
    ) {
      
      make_metabolite_trajectory(
        metab_id = metab_id,
        display_name = DisplayName,
        sex_label = PreferredSex,
        pattern_label = Pattern,
        biological_axis =
          as.character(
            BiologicalAxis
          )
      )
    }
  )
  
  p_fig3c <- wrap_plots(
    curated_plots,
    ncol = 4
  )
  
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
# 3. Fig3D — chemical-class directional remodeling
#
# This provides an independent middle layer between individual
# metabolites and KEGG pathways.
#
# We summarize ALL significant annotated metabolites within
# each chemical superclass:
#
#   HitN            = significant metabolites in the class
#   UpN             = Effect > 0
#   DownN           = Effect < 0
#   DirectionScore  = (UpN - DownN) / HitN
#   MeanEffect      = mean Effect
#
# This is not pathway enrichment.
# ============================================================

chemical_contrasts <- c(
  "Sex_intact",
  "Female_loss",
  "Female_replacement",
  "Male_loss",
  "Male_replacement"
)


chemical_contrast_labels <- c(
  Sex_intact = "Sex difference\nIntact",
  Female_loss = "Female\nGonadectomy",
  Female_replacement = "Female\nReplacement",
  Male_loss = "Male\nGonadectomy",
  Male_replacement = "Male\nReplacement"
)


chemical_da <- da_all %>%
  filter(
    Contrast %in%
      chemical_contrasts,
    Significant,
    !is.na(
      SuperClass
    ),
    SuperClass != "",
    SuperClass != "-"
  ) %>%
  group_by(
    Contrast,
    SuperClass
  ) %>%
  summarise(
    HitN =
      n_distinct(
        metab_id
      ),
    
    UpN =
      n_distinct(
        metab_id[
          Effect > 0
        ]
      ),
    
    DownN =
      n_distinct(
        metab_id[
          Effect < 0
        ]
      ),
    
    MeanEffect =
      mean(
        Effect,
        na.rm = TRUE
      ),
    
    DirectionScore =
      (
        UpN -
          DownN
      ) /
      HitN,
    
    .groups = "drop"
  )


# retain the most informative chemical classes:
# classes must contribute at least 5 significant metabolites
# in at least one contrast
chemical_keep <- chemical_da %>%
  group_by(
    SuperClass
  ) %>%
  summarise(
    MaxHitN =
      max(
        HitN,
        na.rm = TRUE
      ),
    MeanAbsDirection =
      mean(
        abs(
          DirectionScore
        ),
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  filter(
    MaxHitN >= 5
  ) %>%
  arrange(
    desc(
      MaxHitN
    ),
    desc(
      MeanAbsDirection
    )
  ) %>%
  slice_head(
    n = 10
  ) %>%
  pull(
    SuperClass
  )


chemical_plot_df <- chemical_da %>%
  filter(
    SuperClass %in%
      chemical_keep
  ) %>%
  complete(
    Contrast =
      chemical_contrasts,
    SuperClass =
      chemical_keep,
    fill = list(
      HitN = 0,
      UpN = 0,
      DownN = 0,
      MeanEffect = 0,
      DirectionScore = 0
    )
  ) %>%
  mutate(
    Contrast = factor(
      Contrast,
      levels =
        chemical_contrasts
    ),
    
    SuperClass = factor(
      SuperClass,
      levels = rev(
        chemical_keep
      )
    )
  )


write.csv(
  chemical_plot_df,
  file.path(
    TABLE_CANDIDATE_DIR,
    "serum_chemical_superclass_directional_summary.csv"
  ),
  row.names = FALSE
)


p_fig3d <- ggplot(
  chemical_plot_df,
  aes(
    x = Contrast,
    y = SuperClass
  )
) +
  geom_point(
    aes(
      size = HitN,
      fill = DirectionScore
    ),
    shape = 21,
    stroke = 0.35
  ) +
  scale_x_discrete(
    labels =
      chemical_contrast_labels,
    drop = FALSE
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Direction\nscore"
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Significant\nmetabolites",
    title = "Chemical classes show distinct directional remodeling across hormone contrasts",
    subtitle = "Direction score = (upregulated - downregulated) / significant metabolites"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "right"
  )


ggsave(
  file.path(
    PLOT_MAIN_DIR,
    "Fig3D_chemical_superclass_directional_remodeling.pdf"
  ),
  p_fig3d,
  width = 10.5,
  height = 6.8
)


# ============================================================
# 4. Fig3E — directional KEGG pathway remodeling
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
      !is.na(
        description
      ) &
        description != "",
      description,
      pathway_id
    )
  ) %>%
  filter(
    Contrast %in%
      PATHWAY_CONTRASTS,
    PathwayName %in%
      BIOLOGICAL_PATHWAYS
  ) %>%
  complete(
    Contrast =
      PATHWAY_CONTRASTS,
    PathwayName =
      BIOLOGICAL_PATHWAYS,
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
    Contrast = factor(
      Contrast,
      levels =
        PATHWAY_CONTRASTS
    ),
    
    PathwayName = factor(
      PathwayName,
      levels = rev(
        BIOLOGICAL_PATHWAYS
      )
    )
  )


write.csv(
  pathway_bio,
  file.path(
    TABLE_PATHWAY_DIR,
    "serum_biologically_curated_pathways_final.csv"
  ),
  row.names = FALSE
)


p_fig3e <- ggplot(
  pathway_bio,
  aes(
    x = MeanEffect,
    y = PathwayName
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.45
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = MeanEffect,
      yend = PathwayName
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      size = HitN
    ),
    shape = 21,
    stroke = 0.35
  ) +
  facet_wrap(
    ~ Contrast,
    nrow = 1,
    labeller =
      as_labeller(
        PATHWAY_LABELS
      )
  ) +
  labs(
    x = "Mean effect of significant mapped metabolites",
    y = NULL,
    size = "Mapped hits",
    title = "Directional pathway-level support for sex-divergent metabolic remodeling",
    subtitle = "Positive and negative values indicate the average direction of significant mapped metabolites"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
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
# 5. Final summary
# ============================================================

cat(
  "\n\n========================================\n",
  "FINAL RESULT 3 FIVE-PANEL BLOCK DONE\n",
  "========================================\n",
  sep = ""
)

cat(
  "\nMain figure panels:\n",
  "A. Keep Fig3A alluvial from the main script\n",
  "B. Fig3B_simplified_response_modes.pdf\n",
  "C. Fig3C_biologically_curated_metabolite_trajectories.pdf\n",
  "D. Fig3D_chemical_superclass_directional_remodeling.pdf\n",
  "E. Fig3E_directional_pathway_support.pdf\n",
  sep = ""
)

cat(
  "\nSuggested supplementary figures:\n",
  "- Fig3A sex-effect scatter\n",
  "- full six-class response architecture\n",
  "- previous curated module trajectory plot\n",
  sep = ""
)

cat(
  "\nDone.\n"
)
