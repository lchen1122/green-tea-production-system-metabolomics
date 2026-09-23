# ============================================================
# MSI annotation confidence classification
# create a NEW table instead of modifying df
# ============================================================

library(data.table)
library(dplyr)

# ------------------------------------------------------------
# 1. 读取
# ------------------------------------------------------------

df <- fread(
  "F:/科研/有机生态低碳茶叶分析/代谢组结果/1.Preprocess/mix/all_metab_data.txt",
  data.table = FALSE
)

# ------------------------------------------------------------
# 2. 读取 AnnoOverview：保持原表不变
# ------------------------------------------------------------

anno <- fread(
  "F:/科研/有机生态低碳茶叶分析/代谢组结果/3.Anno/04.AnnoOverview/anno.xls",
  data.table = FALSE
)

# ------------------------------------------------------------
# 3. 从 anno 中提取 MSI 判定需要的分类信息
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 4. 新建 MSI 专用数据表
# 不修改原始 df
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 5. 定义注释证据
# ------------------------------------------------------------

msi_df <- msi_df %>%
  mutate(
    
    # 有明确具体代谢物名称
    has_specific_annotation =
      !is.na(Metabolite) &
      trimws(Metabolite) != "" &
      !grepl("^metab_", Metabolite, ignore.case = TRUE),
    
    # 至少有 compound-class 层面的信息
    has_class_annotation =
      (
        (!is.na(super_class) & trimws(super_class) != "" & super_class != "-") |
          (!is.na(class)       & trimws(class)       != "" & class       != "-") |
          (!is.na(sub_class)   & trimws(sub_class)   != "" & sub_class   != "-")
      ),
    
    # B(i): experimental MS/MS spectral matching
    has_experimental_fragmentation =
      level == "B(i)" &
      !is.na(`Fragmentation Score`) &
      `Fragmentation Score` > 0,
    
    # B(ii): fragmentation pattern matching
    has_theoretical_fragmentation =
      level == "B(ii)" &
      !is.na(`Theoretical Fragmentation Score`) &
      `Theoretical Fragmentation Score` > 0,
    
    # mass accuracy criterion
    mass_error_ok =
      !is.na(`Mass Error (ppm)`) &
      abs(`Mass Error (ppm)`) < 10
  )

# ------------------------------------------------------------
# 6. MSI level assignment
# ------------------------------------------------------------

msi_df <- msi_df %>%
  mutate(
    MSI_level = case_when(
      
      # Level 2:
      # specific metabolite annotation + mass error < 10 ppm
      # + fragmentation-based evidence
      has_specific_annotation &
        mass_error_ok &
        (
          has_experimental_fragmentation |
            has_theoretical_fragmentation
        ) ~ "Level 2",
      
      # Level 3:
      # compound class identified, but no specific metabolite identity
      !has_specific_annotation &
        has_class_annotation ~ "Level 3",
      
      # Level 4:
      # unannotated feature
      !has_specific_annotation &
        !has_class_annotation ~ "Level 4",
      
      # specific metabolite annotation available,
      # but mass error exceeds the predefined 10 ppm threshold
      has_specific_annotation &
        !mass_error_ok ~ "Excluded",
      
      TRUE ~ "Review"
    )
  )

# ------------------------------------------------------------
# 7. Annotation evidence
# ------------------------------------------------------------

msi_df <- msi_df %>%
  mutate(
    Annotation_basis = case_when(
      
      MSI_level == "Level 2" &
        has_experimental_fragmentation ~
        "Accurate mass + experimental MS/MS spectral matching",
      
      MSI_level == "Level 2" &
        has_theoretical_fragmentation ~
        "Accurate mass + theoretical fragmentation matching",
      
      MSI_level == "Level 3" ~
        "Compound-class annotation",
      
      MSI_level == "Level 4" ~
        "Unannotated feature",
      
      MSI_level == "Excluded" ~
        "Specific annotation available; mass error >10 ppm",
      
      MSI_level == "Review" ~
        "Manual review required",
      
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 8. 检查最终分类
# ------------------------------------------------------------

table(msi_df$MSI_level, useNA = "ifany")

table(
  Original_level = msi_df$level,
  MSI_level = msi_df$MSI_level,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 9. 生成最终 Supplementary Table S2
# 保留原始信息、样本丰度值和新增 MSI 注释
# ------------------------------------------------------------

Supplementary_Table_S2 <- msi_df %>%
  select(
    -has_specific_annotation,
    -has_class_annotation,
    -has_experimental_fragmentation,
    -has_theoretical_fragmentation,
    -mass_error_ok
  )

# ------------------------------------------------------------
# 10. 导出 CSV
# ------------------------------------------------------------

write.csv(
  Supplementary_Table_S2,
  file = "Supplementary_Table_S2.csv",
  row.names = FALSE,
  na = ""
)





