############################################################
# Result 5. Biochemical assays under HFD conditions
# Long-format input table:
#   Item | Value | Class
# Statistics revised according to Reviewer 3:
#   overall test: ordinary one-way ANOVA
#   post hoc test: Dunnett's multiple-comparison test
#                  using HFD as the reference group
#   plot labels: significance versus HFD only
############################################################

rm(list = ls())
gc()

#==============================
# 0. packages
#==============================
library(tidyverse)
library(readxl)
library(multcomp)

#==============================
# 1. input / output paths
#==============================
bio_file   <- "F:/科研/有机生态低碳茶叶分析/生化实验汇总.xlsx"
sheet_name <- 1
out_dir    <- "F:/科研/有机生态低碳茶叶分析/Paper/Result5"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#==============================
# 2. read data
#==============================
bio_raw <- read_excel(bio_file, sheet = sheet_name) %>%
  rename(
    Item  = 1,
    Value = 2,
    Class = 3
  ) %>%
  mutate(
    Item  = as.character(Item),
    Class = as.character(Class),
    Value = as.numeric(Value)
  ) %>%
  filter(!is.na(Item), !is.na(Value), !is.na(Class))

# Check original marker names
sort(unique(bio_raw$Class))

#==============================
# 3. parse groups + map markers
#==============================
bio_df <- bio_raw %>%
  mutate(
    Group = case_when(
      grepl("^HFD_ORG_", Item) ~ "ORG",
      grepl("^HFD_CON_", Item) ~ "CON",
      grepl("^HFD_ELC_", Item) ~ "ELC",
      grepl("^ND_", Item)      ~ "ND",
      grepl("^HFD_", Item)     ~ "HFD",
      TRUE ~ NA_character_
    ),
    Rep = readr::parse_number(Item)
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(
    # Plotting order is kept as ND, HFD, CON, ELC, ORG.
    # HFD will be re-leveled as the reference only inside the statistical model.
    Group = factor(Group, levels = c("ND", "HFD", "CON", "ELC", "ORG")),
    
    DisplayName = case_when(
      str_detect(Class, regex("^ROS|Reactive oxygen|ROS", ignore_case = TRUE)) ~ "ROS",
      str_detect(Class, regex("^MDA|MDA", ignore_case = TRUE)) ~ "MDA",
      str_detect(Class, regex("^GSSG|GSSG", ignore_case = TRUE)) ~ "GSSG",
      str_detect(Class, regex("GSH[-– ]?Px|GPX", ignore_case = TRUE)) ~ "GPX",
      str_detect(Class, regex("^CAT|CAT", ignore_case = TRUE)) ~ "CAT",
      str_detect(Class, regex("T[-– ]?AOC", ignore_case = TRUE)) ~ "T-AOC",
      str_detect(Class, regex("^GSH\\s*\\(|T[-– ]?GSH", ignore_case = TRUE)) ~ "T-GSH",
      str_detect(Class, regex("T[-– ]?SOD", ignore_case = TRUE)) ~ "T-SOD",
      str_detect(Class, regex("NEFA", ignore_case = TRUE)) ~ "NEFA",
      TRUE ~ NA_character_
    ),
    
    Panel = case_when(
      DisplayName %in% c("ROS", "MDA", "GSSG") ~ "A",
      DisplayName %in% c("GPX", "CAT", "T-AOC", "T-GSH", "T-SOD") ~ "B",
      TRUE ~ NA_character_
    ),
    
    MarkerOrder = case_when(
      DisplayName == "ROS"   ~ 1,
      DisplayName == "MDA"   ~ 2,
      DisplayName == "GSSG"  ~ 3,
      DisplayName == "GPX"   ~ 1,
      DisplayName == "CAT"   ~ 2,
      DisplayName == "T-AOC" ~ 3,
      DisplayName == "T-GSH" ~ 4,
      DisplayName == "T-SOD" ~ 5,
      DisplayName == "NEFA"  ~ 6,
      TRUE ~ NA_real_
    )
  )

# Check unmapped Class names
unknown_classes <- bio_df %>%
  filter(is.na(DisplayName)) %>%
  distinct(Class)

if (nrow(unknown_classes) > 0) {
  cat("\nUnmapped Class names:\n")
  print(unknown_classes)
}

#==============================
# 4. keep markers for main analysis
#    NEFA is excluded from the main figure
#==============================
bio_use <- bio_df %>%
  filter(!is.na(Panel)) %>%
  mutate(
    DisplayName = factor(
      DisplayName,
      levels = c("ROS", "MDA", "GSSG", "GPX", "CAT", "T-AOC", "T-GSH", "T-SOD")
    )
  )

write.csv(
  bio_use,
  file.path(out_dir, "Result5_source_data_long.csv"),
  row.names = FALSE
)

#==============================
# 5. summary statistics by group
#==============================
summary_df <- bio_use %>%
  group_by(DisplayName, Panel, MarkerOrder, Group) %>%
  summarise(
    n = n(),
    mean_value = mean(Value, na.rm = TRUE),
    sd_value   = sd(Value, na.rm = TRUE),
    se_value   = sd_value / sqrt(n),
    .groups = "drop"
  )

write.csv(
  summary_df,
  file.path(out_dir, "Result5_summary_by_group.csv"),
  row.names = FALSE
)

#==============================
# 6. one-way ANOVA + Dunnett test
#    HFD is the reference group
#==============================

# Overall ordinary one-way ANOVA
anova_stats <- bio_use %>%
  group_by(DisplayName, Panel, MarkerOrder) %>%
  group_modify(~{
    df_now <- .x %>% drop_na(Value, Group)
    
    fit <- aov(Value ~ Group, data = df_now)
    p_val <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    tibble(anova_p = p_val)
  }) %>%
  ungroup()

# Dunnett's multiple-comparison test:
# ND, CON, ELC, ORG versus HFD
# HFD is re-leveled as the control only for the model.
dunnett_stats <- bio_use %>%
  group_by(DisplayName, Panel, MarkerOrder) %>%
  group_modify(~{
    df_now <- .x %>%
      drop_na(Value, Group) %>%
      mutate(
        Group = droplevels(Group),
        Group = relevel(Group, ref = "HFD")
      )
    
    fit <- aov(Value ~ Group, data = df_now)
    
    dun_fit <- multcomp::glht(
      fit,
      linfct = multcomp::mcp(Group = "Dunnett")
    )
    
    dun_sum <- summary(dun_fit)
    dun_ci  <- confint(dun_fit)$confint
    
    comparison_names <- names(dun_sum$test$coefficients)
    
    tibble(
      comparison = comparison_names,
      Group = sub(" - HFD$", "", comparison_names),
      estimate = unname(dun_sum$test$coefficients),
      SE = unname(dun_sum$test$sigma),
      t_value = unname(dun_sum$test$tstat),
      p_dunnett = unname(dun_sum$test$pvalues),
      CI_low = dun_ci[, "lwr"],
      CI_high = dun_ci[, "upr"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    Group = factor(Group, levels = c("ND", "CON", "ELC", "ORG")),
    Significance = case_when(
      p_dunnett < 0.001 ~ "***",
      p_dunnett < 0.01  ~ "**",
      p_dunnett < 0.05  ~ "*",
      TRUE              ~ "ns"
    )
  )

write.csv(
  anova_stats,
  file.path(out_dir, "Result5_oneway_ANOVA.csv"),
  row.names = FALSE
)

write.csv(
  dunnett_stats,
  file.path(out_dir, "Result5_Dunnett_vs_HFD.csv"),
  row.names = FALSE
)

#==============================
# 7. combined statistics table
#==============================
result5_stats_table <- dunnett_stats %>%
  left_join(
    anova_stats,
    by = c("DisplayName", "Panel", "MarkerOrder")
  ) %>%
  arrange(Panel, MarkerOrder, Group)

write.csv(
  result5_stats_table,
  file.path(out_dir, "Result5_stats_ANOVA_Dunnett.csv"),
  row.names = FALSE
)

#==============================
# 8. supplementary-style summary table
#==============================
summary_wide <- summary_df %>%
  mutate(
    mean_se = sprintf("%.3f ± %.3f", mean_value, se_value)
  ) %>%
  dplyr::select(
    DisplayName,
    Panel,
    MarkerOrder,
    Group,
    mean_se
  ) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = mean_se
  )

dunnett_wide <- dunnett_stats %>%
  dplyr::select(
    DisplayName,
    Panel,
    MarkerOrder,
    Group,
    p_dunnett
  ) %>%
  dplyr::mutate(
    p_name = paste0(
      "Dunnett_",
      as.character(Group),
      "_vs_HFD"
    )
  ) %>%
  dplyr::select(
    -Group
  ) %>%
  tidyr::pivot_wider(
    names_from = p_name,
    values_from = p_dunnett
  )

supp_table_s7 <- summary_wide %>%
  left_join(
    anova_stats,
    by = c("DisplayName", "Panel", "MarkerOrder")
  ) %>%
  left_join(
    dunnett_wide,
    by = c("DisplayName", "Panel", "MarkerOrder")
  ) %>%
  arrange(Panel, MarkerOrder)

write.csv(
  supp_table_s7,
  file.path(out_dir, "Supplementary_Table_S7_Result5_biochemistry.csv"),
  row.names = FALSE
)

#==============================
# 9. label positions for Dunnett significance
#==============================

# Only ND, CON, ELC and ORG receive labels because all are compared with HFD.
# HFD itself is the reference and therefore has no significance label.
dunnett_label_df <- summary_df %>%
  dplyr::inner_join(
    dunnett_stats %>%
      dplyr::select(
        DisplayName,
        Panel,
        MarkerOrder,
        Group,
        p_dunnett,
        Significance
      ),
    by = c(
      "DisplayName",
      "Panel",
      "MarkerOrder",
      "Group"
    )
  ) %>%
  dplyr::group_by(
    DisplayName,
    Panel,
    MarkerOrder
  ) %>%
  dplyr::mutate(
    y_range =
      max(mean_value + se_value, na.rm = TRUE) -
      min(mean_value - se_value, na.rm = TRUE),
    
    y_label =
      mean_value +
      se_value +
      0.07 * (y_range + 0.01)
  ) %>%
  dplyr::ungroup()

#==============================
# 10. colors
#==============================
group_cols <- c(
  "ND"  = "#BDBDBD",
  "HFD" = "#7A7A7A",
  "CON" = "#E64B35",
  "ELC" = "#00A087",
  "ORG" = "#4DBBD5"
)

#==============================
# 11. selected markers for plotting
#==============================
selected_markers <- c("ROS", "MDA", "T-SOD", "GPX", "CAT", "T-AOC")

#==============================
# 12. plot function
#     6 markers, 2 rows x 3 cols
#     significance labels are Dunnett-adjusted P values vs HFD
#==============================
plot_result5_selected <- function(summary_df, dunnett_label_df) {
  
  summary_sub <- summary_df %>%
    filter(as.character(DisplayName) %in% selected_markers) %>%
    mutate(
      DisplayName = factor(as.character(DisplayName), levels = selected_markers)
    )
  
  label_sub <- dunnett_label_df %>%
    filter(as.character(DisplayName) %in% selected_markers) %>%
    mutate(
      DisplayName = factor(as.character(DisplayName), levels = selected_markers)
    )
  
  p <- ggplot(summary_sub, aes(x = Group, y = mean_value, fill = Group)) +
    geom_col(width = 0.68, color = "black", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
      width = 0.16,
      linewidth = 0.55
    ) +
    geom_text(
      data = label_sub,
      aes(x = Group, y = y_label, label = Significance),
      inherit.aes = FALSE,
      size = 4.0,
      fontface = "bold",
      vjust = 0
    ) +
    facet_wrap(~ DisplayName, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_cols) +
    labs(
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.22))
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5)
    )
  
  return(p)
}

#==============================
# 13. Figure 5
#==============================
p_fig5 <- plot_result5_selected(
  summary_df = summary_df,
  dunnett_label_df = dunnett_label_df
)

p_fig5

ggsave(
  file.path(out_dir, "Figure5_selected6_markers_Dunnett_vs_HFD.pdf"),
  p_fig5,
  width = 10,
  height = 6.5
)

ggsave(
  file.path(out_dir, "Figure5_selected6_markers_Dunnett_vs_HFD.png"),
  p_fig5,
  width = 10,
  height = 6.5,
  dpi = 300
)

#==============================
# 14. print selected Dunnett results
#==============================
selected_dunnett <- dunnett_stats %>%
  filter(as.character(DisplayName) %in% selected_markers) %>%
  arrange(
    factor(as.character(DisplayName), levels = selected_markers),
    Group
  )

print(selected_dunnett)

#==============================
# 15. done
#==============================
cat("\nDone.\n")
cat("Main outputs:\n")
cat("1) Figure5_selected6_markers_Dunnett_vs_HFD.pdf / .png\n")
cat("2) Result5_oneway_ANOVA.csv\n")
cat("3) Result5_Dunnett_vs_HFD.csv\n")
cat("4) Result5_stats_ANOVA_Dunnett.csv\n")
cat("5) Supplementary_Table_S7_Result5_biochemistry.csv\n")
