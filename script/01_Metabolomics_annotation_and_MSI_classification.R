# ============================================================
# 01_metabolomics_annotation_confidence.R
# Annotation confidence review for untargeted metabolomics
#
# Revision logic:
# - MS1 mass error: |ppm| < 10
# - B(i): experimental MS/MS, Fragmentation Score > 35
#         -> MSI Level 2 (putatively annotated compound)
# - B(ii): theoretical fragmentation, Theoretical Score > 40
#          -> retained as lower-confidence theoretical annotation,
#             NOT counted as MSI Level 2
# - No authentic-standard confirmation -> no MSI Level 1 assignment
# ============================================================

rm(list = ls())
gc()

library(data.table)
library(dplyr)

# ============================================================
# 1. Paths
# ============================================================

project_dir <- "F:/科研/有机生态低碳茶叶分析"

data_file <- file.path(
  project_dir,
  "代谢组结果/1.Preprocess/mix/all_metab_data.txt"
)

anno_file <- file.path(
  project_dir,
  "代谢组结果/3.Anno/04.AnnoOverview/anno.xls"
)

out_dir <- file.path(
  project_dir,
  "Paper/Revised2/01_annotation"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read data
# ============================================================

df <- fread(
  data_file,
  data.table = FALSE,
  check.names = FALSE
)

anno <- fread(
  anno_file,
  data.table = FALSE,
  check.names = FALSE
)

anno_class <- anno %>%
  select(
    metab_id,
    hmdb_name,
    compound_name,
    compound_first_category,
    compound_second_category,
    kingdom,
    super_class,
    class,
    sub_class
  ) %>%
  distinct(metab_id, .keep_all = TRUE)

msi_df <- df %>%
  left_join(anno_class, by = "metab_id") %>%
  mutate(
    `Fragmentation Score` =
      suppressWarnings(as.numeric(`Fragmentation Score`)),
    `Theoretical Fragmentation Score` =
      suppressWarnings(as.numeric(`Theoretical Fragmentation Score`)),
    `Mass Error (ppm)` =
      suppressWarnings(as.numeric(`Mass Error (ppm)`))
  )

# ============================================================
# 3. Annotation evidence
# ============================================================

msi_df <- msi_df %>%
  mutate(
    has_specific_annotation =
      !is.na(Metabolite) &
      trimws(Metabolite) != "" &
      !grepl("^metab_", Metabolite, ignore.case = TRUE),

    has_class_annotation =
      (
        (!is.na(super_class) & trimws(super_class) != "" & super_class != "-") |
        (!is.na(class)       & trimws(class)       != "" & class       != "-") |
        (!is.na(sub_class)   & trimws(sub_class)   != "" & sub_class   != "-")
      ),

    mass_error_ok =
      !is.na(`Mass Error (ppm)`) &
      abs(`Mass Error (ppm)`) < 10,

    experimental_msms_ok =
      level == "B(i)" &
      !is.na(`Fragmentation Score`) &
      `Fragmentation Score` > 35,

    theoretical_fragmentation_ok =
      level == "B(ii)" &
      !is.na(`Theoretical Fragmentation Score`) &
      `Theoretical Fragmentation Score` > 40
  )

# ============================================================
# 4. Revised annotation confidence classification
# ============================================================

msi_df <- msi_df %>%
  mutate(
    MSI_level = case_when(
      has_specific_annotation &
        mass_error_ok &
        experimental_msms_ok ~ "Level 2",

      !has_specific_annotation &
        has_class_annotation ~ "Level 3",

      !has_specific_annotation &
        !has_class_annotation ~ "Level 4",

      TRUE ~ NA_character_
    ),

    Annotation_status = case_when(
      has_specific_annotation &
        mass_error_ok &
        experimental_msms_ok ~
        "MSI Level 2: accurate mass + experimental MS/MS spectral matching",

      has_specific_annotation &
        mass_error_ok &
        theoretical_fragmentation_ok ~
        "Tentative theoretical annotation: accurate mass + theoretical fragmentation matching",

      has_specific_annotation &
        !mass_error_ok ~
        "Specific database annotation, but mass error >=10 ppm",

      has_specific_annotation &
        mass_error_ok ~
        "Specific database annotation without accepted fragmentation support",

      !has_specific_annotation &
        has_class_annotation ~
        "MSI Level 3: compound-class annotation",

      TRUE ~
        "MSI Level 4: unknown feature"
    ),

    Accepted_fragmentation_score = case_when(
      experimental_msms_ok ~ `Fragmentation Score`,
      theoretical_fragmentation_ok ~ `Theoretical Fragmentation Score`,
      TRUE ~ NA_real_
    ),

    Annotation_evidence_type = case_when(
      experimental_msms_ok ~ "Experimental MS/MS",
      theoretical_fragmentation_ok ~ "Theoretical fragmentation",
      TRUE ~ "Other / insufficient fragmentation evidence"
    )
  )

# ============================================================
# 5. Summary counts
# ============================================================

annotation_summary <- msi_df %>%
  count(
    Annotation_status,
    MSI_level,
    name = "n_features"
  ) %>%
  arrange(desc(n_features))

print(annotation_summary)

cat("\nMSI Level 2 features:",
    sum(msi_df$MSI_level == "Level 2", na.rm = TRUE), "\n")

cat("Accepted theoretical-only annotations:",
    sum(
      msi_df$theoretical_fragmentation_ok &
        msi_df$mass_error_ok &
        msi_df$has_specific_annotation,
      na.rm = TRUE
    ),
    "\n"
)

cat("No feature is assigned MSI Level 1 because no authentic-standard confirmation was performed.\n")

# ============================================================
# 7. Export
# ============================================================

helper_cols <- c(
  "has_specific_annotation",
  "has_class_annotation",
  "mass_error_ok",
  "experimental_msms_ok",
  "theoretical_fragmentation_ok"
)

Supplementary_Table_S2_revised <- msi_df %>%
  select(-any_of(helper_cols))

fwrite(
  Supplementary_Table_S2_revised,
  file.path(out_dir, "Supplementary_Table_S2_revised.csv"),
  na = ""
)

fwrite(
  annotation_summary,
  file.path(out_dir, "annotation_confidence_summary.csv"),
  na = ""
)


cat("\nDone.\nOutput directory:\n", out_dir, "\n")

# ============================================================
# 7. 检查所有原始 all_metab_data 列是否完整保留
# ============================================================

stopifnot(
  all(names(df) %in% names(Supplementary_Table_S2_revised))
)

# 查看最终附表列名
print(names(Supplementary_Table_S2_revised))

cat(
  "\nOriginal features:", nrow(df),
  "\nFinal S2 features:", nrow(Supplementary_Table_S2_revised),
  "\nOriginal columns:", ncol(df),
  "\nFinal S2 columns:", ncol(Supplementary_Table_S2_revised),
  "\n"
)
