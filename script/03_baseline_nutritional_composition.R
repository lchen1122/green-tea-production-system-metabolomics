# ============================================================
# Untargeted metabolomics vs HPLC targets
# Targets: EGCG, ECG, EGC, EC, Caffeine
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Read untargeted metabolomics data
# ------------------------------------------------------------

df <- fread(
  "F:/科研/有机生态低碳茶叶分析/代谢组结果/1.Preprocess/mix/all_metab_data.txt",
  data.table = FALSE
)

# ------------------------------------------------------------
# 2. Define metabolites corresponding to HPLC targets
# ------------------------------------------------------------

target_map <- data.frame(
  Metabolite = c(
    # EGCG
    "Epigallocatechin-3-Gallate",
    "Epigallocatechin Gallate",
    
    # ECG
    "Epicatechin-3-Gallate",
    "(-)-Epicatechin 3-O-Gallate",
    
    # EGC
    "Epigallocatechin",
    
    # EC
    "Epicatechin",
    
    # Caffeine
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

# ------------------------------------------------------------
# 3. Extract matched metabolites
# ------------------------------------------------------------

target_df <- df %>%
  inner_join(
    target_map,
    by = "Metabolite"
  )

# Check matched features
target_info <- target_df %>%
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
    `Mass Error (ppm)`
  ) %>%
  arrange(
    factor(
      HPLC_target,
      levels = c("EGCG", "ECG", "EGC", "EC", "Caffeine")
    )
  )

print(target_info)

# ------------------------------------------------------------
# 4. Identify biological sample columns
# ------------------------------------------------------------

sample_cols <- grep(
  "^(CON|ELC|ORG)_",
  colnames(target_df),
  value = TRUE
)

# ------------------------------------------------------------
# 5. Convert to long format
# ------------------------------------------------------------

target_long <- target_df %>%
  select(
    HPLC_target,
    metab_id,
    Metabolite,
    all_of(sample_cols)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

# ------------------------------------------------------------
# 6. Merge duplicate annotations
# If one HPLC target has multiple untargeted features,
# use the mean abundance within each sample
# ------------------------------------------------------------

target_combined <- target_long %>%
  group_by(
    HPLC_target,
    Sample
  ) %>%
  summarise(
    Abundance = mean(
      as.numeric(Abundance),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Assign production system
# ------------------------------------------------------------

target_combined <- target_combined %>%
  mutate(
    Group = case_when(
      grepl("^CON_", Sample) ~ "CON",
      grepl("^ELC_", Sample) ~ "ELC",
      grepl("^ORG_", Sample) ~ "ORG"
    ),
    
    Group = factor(
      Group,
      levels = c("CON", "ELC", "ORG")
    ),
    
    HPLC_target = factor(
      HPLC_target,
      levels = c(
        "EGCG",
        "ECG",
        "EGC",
        "EC",
        "Caffeine"
      )
    )
  )

# ------------------------------------------------------------
# 8. Scale within each metabolite
# 0-1 scaling for trend comparison
# ------------------------------------------------------------

target_scaled <- target_combined %>%
  group_by(HPLC_target) %>%
  mutate(
    min_value = min(Abundance, na.rm = TRUE),
    max_value = max(Abundance, na.rm = TRUE),
    
    Scaled_abundance = ifelse(
      max_value == min_value,
      0.5,
      (Abundance - min_value) /
        (max_value - min_value)
    )
  ) %>%
  ungroup() %>%
  select(
    -min_value,
    -max_value
  )

# ------------------------------------------------------------
# 9. Calculate mean and SD for CON / ELC / ORG
# ------------------------------------------------------------

target_summary <- target_scaled %>%
  group_by(
    HPLC_target,
    Group
  ) %>%
  summarise(
    Mean = mean(
      Scaled_abundance,
      na.rm = TRUE
    ),
    SD = sd(
      Scaled_abundance,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

print(target_summary)

# ------------------------------------------------------------
# 10. Plot faceted bar figure
# ------------------------------------------------------------

p_untargeted <- ggplot(
  target_summary,
  aes(
    x = Group,
    y = Mean
  )
) +
  geom_col(
    width = 0.68,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(Mean - SD, 0),
      ymax = pmin(Mean + SD, 1)
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ HPLC_target,
    nrow = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Scaled relative abundance"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.title.y = element_text(
      size = 11
    ),
    panel.spacing = unit(
      1,
      "lines"
    )
  )

print(p_untargeted)

# ============================================================
# HPLC results
# Faceted bar plot: EGCG, ECG, EGC, EC, Caffeine
# Unit converted from % to mg/g
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Read HPLC data
# ------------------------------------------------------------

hplc <- read_excel(
  "F:/科研/有机生态低碳茶叶分析/Paper/HPLC.xlsx"
)

# ------------------------------------------------------------
# 2. Convert unit and define groups
# 1% = 10 mg/g
# ------------------------------------------------------------

hplc <- hplc %>%
  mutate(
    Value_mg_g = Value * 10,
    
    Treatment = case_when(
      grepl("^CON_", Group) ~ "CON",
      grepl("^ELC_", Group) ~ "ELC",
      grepl("^ORG_", Group) ~ "ORG"
    ),
    
    Treatment = factor(
      Treatment,
      levels = c("CON", "ELC", "ORG")
    ),
    
    Class = factor(
      Class,
      levels = c(
        "EGCG",
        "ECG",
        "EGC",
        "EC",
        "Caffeine"
      )
    )
  )

# ------------------------------------------------------------
# 3. Check data
# ------------------------------------------------------------

head(hplc)

table(
  hplc$Class,
  hplc$Treatment,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 4. Calculate mean ± SD
# ------------------------------------------------------------

hplc_summary <- hplc %>%
  group_by(
    Class,
    Treatment
  ) %>%
  summarise(
    Mean = mean(
      Value_mg_g,
      na.rm = TRUE
    ),
    SD = sd(
      Value_mg_g,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  )

hplc_summary

# ------------------------------------------------------------
# 5. Faceted bar plot
# ------------------------------------------------------------

p_hplc <- ggplot(
  hplc_summary,
  aes(
    x = Treatment,
    y = Mean
  )
) +
  geom_col(
    width = 0.65,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = Mean - SD,
      ymax = Mean + SD
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ Class,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "Content (mg/g)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.title.y = element_text(
      size = 11
    ),
    panel.spacing = unit(
      0.8,
      "lines"
    )
  )

# ------------------------------------------------------------
# 6. Show plot
# ------------------------------------------------------------

p_hplc

# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "HPLC_faceted_bar_mg_g.pdf",
  plot = p_hplc,
  width = 10,
  height = 3.6
)

ggsave(
  filename = "HPLC_faceted_bar_mg_g.png",
  plot = p_hplc,
  width = 10,
  height = 3.6,
  dpi = 600
)


library(patchwork)
library(ggplot2)

# ------------------------------------------------------------
# 1. combine
# ------------------------------------------------------------

p_combined <- p_untargeted / p_hplc +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

# 显示
p_combined

# ------------------------------------------------------------
# 2. export
# ------------------------------------------------------------

ggsave(
  filename = "Untarget_HPLC_combined.pdf",
  plot = p_combined,
  width = 10,
  height = 7.2
)

ggsave(
  filename = "Untarget_HPLC_combined.png",
  plot = p_combined,
  width = 10,
  height = 7.2,
  dpi = 600
)




# ------------------------------------------------------------
# 8. Export summary table
# ------------------------------------------------------------

write.csv(
  hplc_summary,
  file = "HPLC_group_summary_mg_g.csv",
  row.names = FALSE
)
