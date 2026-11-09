# ============================================================
# 03_baseline_nutritional_composition.R
# Conventional nutritional composition of CON / ELC / ORG teas
# ============================================================

rm(list = ls())
gc()

library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(pheatmap)

# ============================================================
# 1. Paths
# ============================================================

project_dir <- "F:/科研/有机生态低碳茶叶分析"

data_file <- file.path(
  project_dir,
  "常规营养检测结果.xlsx"
)

out_dir <- file.path(
  project_dir,
  "Paper/Revised2/03_baseline_nutrition"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read and split data
# ============================================================

nutrient <- read_excel(data_file)

nutrient <- nutrient %>%
  mutate(
    Class = case_when(
      Item %in% c(
        "Ca", "Mn", "Fe", "Cu", "Zn", "Se", "Sr"
      ) ~ "Micro",
      
      Item %in% c(
        "Fat",
        "Protein",
        "Moisture",
        "Ash",
        "Dietary fiber",
        "Carbohydrates"
      ) ~ "Macro",
      
      TRUE ~ "Protein"
    )
  )

macronutrient <- nutrient %>% filter(Class == "Macro")
microelement  <- nutrient %>% filter(Class == "Micro")
protein       <- nutrient %>% filter(Class == "Protein")

group_colors <- c(
  "CON" = "#E64B35",
  "ELC" = "#00A087",
  "ORG" = "#4DBBD5"
)

# ============================================================
# 3. Helper: wide -> long
# ============================================================

to_long <- function(dat) {
  dat %>%
    pivot_longer(
      cols = -c(Item, Unit, Class),
      names_to = c("Group", "Rep"),
      names_sep = "_",
      values_to = "Value"
    ) %>%
    mutate(
      Rep = as.integer(Rep),
      Group = factor(Group, levels = c("CON", "ELC", "ORG"))
    )
}

# ============================================================
# 4. Protein / hydrolyzed amino-acid related measurements
# ============================================================

protein_long <- to_long(protein)

protein_summary <- protein_long %>%
  group_by(Item, Unit, Class, Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    .groups = "drop"
  )

protein_test <- protein_long %>%
  group_by(Item) %>%
  anova_test(dv = Value, between = Group) %>%
  ungroup()

write.csv(
  protein_summary,
  file.path(out_dir, "protein_amino_acid_summary.csv"),
  row.names = FALSE
)

write.csv(
  protein_test,
  file.path(out_dir, "protein_amino_acid_ANOVA.csv"),
  row.names = FALSE
)

p_protein <- ggplot(
  protein_summary,
  aes(x = Group, y = mean, fill = Group)
) +
  geom_col(
    width = 0.65,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.15,
    linewidth = 0.4
  ) +
  geom_point(
    data = protein_long,
    aes(x = Group, y = Value),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.10, height = 0),
    size = 1.8
  ) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ Item, scales = "free_y") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )
p_protein

ggsave(
  file.path(out_dir, "baseline_protein_amino_acids.pdf"),
  p_protein,
  width = 10,
  height = 8
)

# ============================================================
# 5. Mineral elements
# ============================================================

micro_long <- to_long(microelement)

micro_summary <- micro_long %>%
  group_by(Item, Unit, Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    .groups = "drop"
  )

write.csv(
  micro_summary,
  file.path(out_dir, "mineral_element_summary.csv"),
  row.names = FALSE
)

micro_mat <- microelement %>%
  select(
    Item,
    matches("^(CON|ELC|ORG)_")
  ) %>%
  as.data.frame()

rownames(micro_mat) <- micro_mat$Item
micro_mat$Item <- NULL

micro_mat <- micro_mat[
  ,
  c(
    sort(grep("^CON_", names(micro_mat), value = TRUE)),
    sort(grep("^ELC_", names(micro_mat), value = TRUE)),
    sort(grep("^ORG_", names(micro_mat), value = TRUE))
  ),
  drop = FALSE
]

micro_mat <- as.matrix(micro_mat)
mode(micro_mat) <- "numeric"

micro_log <- log10(micro_mat + 1e-6)
micro_z <- t(scale(t(micro_log)))
micro_z[is.na(micro_z)] <- 0

annotation_col <- data.frame(
  Group = factor(
    sub("_\\d+$", "", colnames(micro_z)),
    levels = c("CON", "ELC", "ORG")
  )
)

rownames(annotation_col) <- colnames(micro_z)

annotation_colors <- list(
  Group = group_colors
)


pheatmap(
  micro_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  filename = file.path(
    out_dir,
    "baseline_mineral_heatmap2.pdf"
  ),
  width = 9,
  height = 7
)

# ============================================================
# 6. Proximate / macronutrient composition
# ============================================================

macro_long <- to_long(macronutrient) %>%
  mutate(
    Item = fct_reorder(Item, Value, .fun = mean)
  )

macro_summary <- macro_long %>%
  group_by(Item, Unit, Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    .groups = "drop"
  )

macro_test <- macro_long %>%
  group_by(Item) %>%
  anova_test(dv = Value, between = Group) %>%
  ungroup()

write.csv(
  macro_summary,
  file.path(out_dir, "macronutrient_summary.csv"),
  row.names = FALSE
)

write.csv(
  macro_test,
  file.path(out_dir, "macronutrient_ANOVA.csv"),
  row.names = FALSE
)

pd <- position_dodge(width = 0.6)

p_macro <- ggplot(
  macro_summary,
  aes(x = mean, y = Item, color = Group)
) +
  geom_errorbar(
    aes(
      xmin = mean - sd,
      xmax = mean + sd
    ),
    position = pd,
    width = 0.18,
    linewidth = 0.6,
    orientation = "y"
  ) +
  geom_point(
    position = pd,
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  geom_point(
    data = macro_long,
    aes(
      x = Value,
      y = Item,
      color = Group
    ),
    position = position_jitterdodge(
      dodge.width = 0.6,
      jitter.width = 0,
      jitter.height = 0.12
    ),
    size = 1.8,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = "Content",
    y = NULL,
    color = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal"
  )
p_macro

ggsave(
  file.path(out_dir, "baseline_macronutrients.pdf"),
  p_macro,
  width = 8,
  height = 6
)

cat("\nDone.\nOutput directory:\n", out_dir, "\n")
