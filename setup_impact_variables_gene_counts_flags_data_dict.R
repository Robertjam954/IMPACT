#!/usr/bin/env Rscript
# Setup variables for IMPACT merged dataset
# - Reads merged patient-level IMPACT CSV
# - Reads sample mutation data and aggregates per sample
# - Recodes variables as requested and writes cleaned dataset + data dictionary

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(tidyr)
})

# Paths (edit if needed)
merged_path <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Data", "merged_impact.csv")
mutation_path <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Data", "sample mutation data.xlsx")
outdir <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

message('Reading merged dataset: ', merged_path)
df <- read_csv(merged_path, guess_max = 10000)
message('Rows:', nrow(df), 'Cols:', ncol(df))

# Read mutation file header and sample of columns
message('Reading mutation file: ', mutation_path)
mut <- read_excel(mutation_path)
message('Mutation rows:', nrow(mut), 'Cols:', ncol(mut))

# --- Merge notes ---
# merged dataset sample id column: 'SAMPLE_ID'
# mutation sample id column: 'Tumor_Sample_Barcode'

# Aggregate mutation-level data to sample-level
mut_summary <- mut %>%
  mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
         gene_name = as.character(gene_name)) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    mutation_count_all_sites_sum = n(),
    t_alt_count_max = ifelse(all(is.na(t_alt_count)), NA_real_, max(as.numeric(t_alt_count), na.rm = TRUE)),
    genes = paste(unique(gene_name), collapse = ';')
  ) %>%
  ungroup()

# Merge summaries into patient-level merged df by SAMPLE_ID -> Tumor_Sample_Barcode
# Ensure SAMPLE_ID and Tumor_Sample_Barcode are character and comparable
if ('SAMPLE_ID' %in% names(df)) {
  df <- df %>% mutate(SAMPLE_ID = as.character(SAMPLE_ID))
  df <- df %>% left_join(mut_summary, by = c('SAMPLE_ID' = 'Tumor_Sample_Barcode'))
} else {
  stop('SAMPLE_ID not found in merged dataset; update the script to use correct sample ID column')
}

# --- Variable recoding per instructions ---
# 1) sample type binary: 1 Primary, 2 Metastasis
if ('SAMPLE_TYPE' %in% names(df)) {
  df <- df %>% mutate(
    sample_type_bin = case_when(
      str_to_lower(SAMPLE_TYPE) == 'primary' ~ 1L,
      str_to_lower(SAMPLE_TYPE) == 'metastasis' ~ 2L,
      TRUE ~ NA_integer_
    )
  )
}

# 2) sample site categorical (unordered)
if ('SAMPLE_SITE' %in% names(df)) {
  df <- df %>% mutate(sample_site = as.character(SAMPLE_SITE)) %>%
    mutate(sample_site = ifelse(trimws(sample_site)=="", NA_character_, sample_site)) %>%
    mutate(sample_site = factor(sample_site, ordered = FALSE))
}

# 3) tumor sample histology categorical (unordered)
if ('TUMOR_SAMPLE_HISTOLOGY' %in% names(df)) {
  df <- df %>% mutate(tumor_sample_histology = factor(TUMOR_SAMPLE_HISTOLOGY, ordered = FALSE))
}

# 4) Keep continuous invasive carcinoma age AND create age_cat (<50, 50-70, >70)
if ('INVASIVE_CARCINOMA_DX_AGE' %in% names(df)) {
  df <- df %>% mutate(Invasive_carcinoma_dx_age = as.numeric(INVASIVE_CARCINOMA_DX_AGE)) %>%
    mutate(age_cat = case_when(
      !is.na(Invasive_carcinoma_dx_age) & Invasive_carcinoma_dx_age < 50 ~ '<50',
      !is.na(Invasive_carcinoma_dx_age) & Invasive_carcinoma_dx_age >= 50 & Invasive_carcinoma_dx_age <= 70 ~ '50-70',
      !is.na(Invasive_carcinoma_dx_age) & Invasive_carcinoma_dx_age > 70 ~ '>70',
      TRUE ~ NA_character_
    )) %>%
    mutate(age_cat = factor(age_cat, levels = c('<50','50-70','>70'), ordered = TRUE))
} else {
  warning('INVASIVE_CARCINOMA_DX_AGE not found; age variables not created')
}

# 5) overall tumor grade ordered: 1 low, 2 intermediate, 3 high
if ('OVERALL_TUMOR_GRADE' %in% names(df)) {
  df <- df %>% mutate(
    overall_tumor_grade_ord = case_when(
      # Grade I: look for roman numeral I or the word 'Low'
      str_detect(OVERALL_TUMOR_GRADE, regex('\\bI\\b|Low', ignore_case = TRUE)) ~ 1L,
      # Grade II: roman numeral II or 'Intermediate'
      str_detect(OVERALL_TUMOR_GRADE, regex('\\bII\\b|Intermediate', ignore_case = TRUE)) ~ 2L,
      # Grade III: roman numeral III or 'High'
      str_detect(OVERALL_TUMOR_GRADE, regex('\\bIII\\b|High', ignore_case = TRUE)) ~ 3L,
      TRUE ~ NA_integer_
    )
  )
  df$overall_tumor_grade_ord <- factor(df$overall_tumor_grade_ord, levels = c(1,2,3), labels = c('Low','Intermediate','High'), ordered = TRUE)
}

# 6) Stage diagnosis ordered mapping
if ('STAGE_AT_DIAGNOSIS' %in% names(df)) {
  df <- df %>% mutate(
    stage_diag_group = case_when(
      str_to_upper(STAGE_AT_DIAGNOSIS) %in% c('IA','IB') ~ 1L,
      str_to_upper(STAGE_AT_DIAGNOSIS) %in% c('IIA','IIB') ~ 2L,
      str_to_upper(STAGE_AT_DIAGNOSIS) %in% c('IIIA','IIIB','IIIC') ~ 3L,
      str_to_upper(STAGE_AT_DIAGNOSIS) %in% c('IV') ~ 4L,
      TRUE ~ NA_integer_
    )
  )
  df$stage_diag_group <- factor(df$stage_diag_group, levels = c(1,2,3,4), labels = c('Stage I','Stage II','Stage III','Stage IV'), ordered = TRUE)
}

# 7) receptor status primary mapping
if ('RECEPTOR_STATUS_PRIMARY' %in% names(df)) {
  df <- df %>% mutate(
    receptor_primary_cat = case_when(
      RECEPTOR_STATUS_PRIMARY %in% c('HR+/HER2-','HR+/HER2_Equivocal') ~ 1L,
      RECEPTOR_STATUS_PRIMARY == 'HR+/HER2+' ~ 2L,
      RECEPTOR_STATUS_PRIMARY %in% c('HR-/HER2+','HR-/HER2_Equivocal') ~ 3L,
      RECEPTOR_STATUS_PRIMARY %in% c('Triple Negative','HR-/HER2-') ~ 4L,
      RECEPTOR_STATUS_PRIMARY %in% c('HR+/HER2_Unknown','HR-/HER2_Unknown') ~ as.integer(NA),
      TRUE ~ as.integer(NA)
    )
  )
  df$receptor_primary_cat <- factor(df$receptor_primary_cat, levels = c(1,2,3,4), labels = c('HR+/HER2-','HR+/HER2+','HR-/HER2+','Triple Negative'), ordered = TRUE)
}

# 8) tmb quartiles (TMB_NONSYNONYMOUS)
if ('TMB_NONSYNONYMOUS' %in% names(df) || 'TMB (nonsynonymous)' %in% names(df)) {
  tmb_col <- ifelse('TMB_NONSYNONYMOUS' %in% names(df), 'TMB_NONSYNONYMOUS', 'TMB (nonsynonymous)')
  df <- df %>% mutate(tmb_val = as.numeric(.data[[tmb_col]]))
  df <- df %>% mutate(tmb_quartile = ntile(tmb_val, 4))
  df$tmb_quartile <- factor(df$tmb_quartile, levels = 1:4, labels = c('Q1','Q2','Q3','Q4'), ordered = TRUE)
}

# 9) FGA (Fraction Genome Altered) quartiles
if ('Fraction Genome Altered' %in% names(df)) {
  df <- df %>% mutate(fga_val = as.numeric(`Fraction Genome Altered`))
  df <- df %>% mutate(fga_quartile = ntile(fga_val, 4))
  df$fga_quartile <- factor(df$fga_quartile, levels = 1:4, labels = c('Q1','Q2','Q3','Q4'), ordered = TRUE)
}

# 10) MSI Type unordered categorical
if ('MSI Type' %in% names(df)) {
  df <- df %>% mutate(MSI_type = case_when(
    str_to_lower(`MSI Type`) == 'stable' ~ 'Stable',
    str_to_lower(`MSI Type`) == 'instable' ~ 'Instable',
    TRUE ~ NA_character_
  )) %>% mutate(MSI_type = factor(MSI_type, ordered = FALSE))
}

# 11) mutation_count_sum -> mutation_count_all_sites_sum (ordered quartiles)
if ('mutation_count_all_sites_sum' %in% names(df)) {
  df <- df %>% mutate(mutation_count_all_sites_sum = as.numeric(mutation_count_all_sites_sum))
} else if ('Mutation Count' %in% names(df)) {
  df <- df %>% mutate(mutation_count_all_sites_sum = as.numeric(`Mutation Count`))
}
if ('mutation_count_all_sites_sum' %in% names(df)) {
  df <- df %>% mutate(mutation_count_q = ntile(mutation_count_all_sites_sum, 4))
  df$mutation_count_q <- factor(df$mutation_count_q, levels=1:4, labels=c('Q1','Q2','Q3','Q4'), ordered=TRUE)
}

# 12) tmb_nonsynonmous
if (exists('tmb_val')) {
  df <- df %>% mutate(tmb_nonsynonmous = tmb_val)
  df <- df %>% mutate(tmb_nonsynonmous_q = ntile(tmb_nonsynonmous, 4))
  df$tmb_nonsynonmous_q <- factor(df$tmb_nonsynonmous_q, levels=1:4, labels=c('Q1','Q2','Q3','Q4'), ordered=TRUE)
}

# 13) t_alt_count quartiles from mutation summary
if ('t_alt_count_max' %in% names(df)) {
  df <- df %>% mutate(t_alt_count = as.numeric(t_alt_count_max))
  df <- df %>% mutate(t_alt_count_q = ntile(t_alt_count, 4))
  df$t_alt_count_q <- factor(df$t_alt_count_q, levels=1:4, labels=c('Q1','Q2','Q3','Q4'), ordered=TRUE)
}

# 14) Race mapping
if ('Race' %in% names(df)) {
  df <- df %>% mutate(Race_clean = case_when(
    Race %in% c('WHITE','White','Middle Eastern') ~ 'White or Caucasian',
    Race %in% c('BLACK OR AFRICAN AMERICAN','Black or African American') ~ 'Black or African American',
    Race %in% c('ASIAN-FAR EAST/INDIAN SUBCONT','Asian-Far East/Indian Subcont','Chinese') ~ 'Asian',
    Race %in% c('NATIVE AMERICAN-AM IND/ALASKA') ~ 'Native American or Alaska Native',
    TRUE ~ NA_character_
  )) %>% mutate(Race = factor(Race_clean, levels=c('White or Caucasian','Black or African American','Asian','Native American or Alaska Native')))
}

# 15) Ethnicity mapping
if ('Ethnicity' %in% names(df)) {
  df <- df %>% mutate(Ethnicity_clean = case_when(
    Ethnicity %in% c('Non-Spanish; Non-Hispanic') ~ 'Non-Hispanic',
    Ethnicity %in% c('Cuban','Dominican Republic','Mexican (includes Chicano)','Puerto Rican','South/Central America (except Brazil)','Spanish  NOS; Hispanic NOS, Latino NOS') ~ 'Hispanic',
    TRUE ~ NA_character_
  )) %>% mutate(Ethnicity = factor(Ethnicity_clean))
}

# 16) Overall HER2 status binary unordered factor
if ('OVERALL_HER2_STATUS' %in% names(df) || 'OVERALL_HER2_STATUS_PATIENT' %in% names(df)) {
  her2col <- ifelse('OVERALL_HER2_STATUS' %in% names(df), 'OVERALL_HER2_STATUS', 'OVERALL_HER2_STATUS_PATIENT')
  df <- df %>% mutate(OVERALL_her2_status = case_when(
    str_detect(as.character(.data[[her2col]]), regex('pos|positive|1', ignore_case=TRUE)) ~ 1L,
    str_detect(as.character(.data[[her2col]]), regex('neg|negative|0', ignore_case=TRUE)) ~ 0L,
    TRUE ~ NA_integer_
  )) %>% mutate(OVERALL_her2_status = factor(OVERALL_her2_status, levels=c(0,1), labels=c('Negative','Positive')))
}

# 17) OS_STATUS mapping
if ('OS_STATUS' %in% names(df)) {
  df <- df %>% mutate(OS_STATUS_f = case_when(
    str_starts(as.character(OS_STATUS), '0') ~ 'Alive',
    str_starts(as.character(OS_STATUS), '1') ~ 'Deceased',
    TRUE ~ NA_character_
  )) %>% mutate(OS_STATUS_f = factor(OS_STATUS_f, levels=c('Alive','Deceased')))
}

# 18) DFS_EVENT mapping
if ('DFS_EVENT' %in% names(df)) {
  df <- df %>% mutate(DFS_EVENT_f = case_when(
    DFS_EVENT == 1 ~ 'Event',
    DFS_EVENT == 0 ~ 'No event',
    TRUE ~ NA_character_
  )) %>% mutate(DFS_EVENT_f = factor(DFS_EVENT_f, levels=c('No event','Event')))
}

# 19) metastatic_recurrence_status
if ('METASTATIC_RECURRENCE_TIME_MONTHS' %in% names(df)) {
  df <- df %>% mutate(metastatic_recurrence_status = case_when(
    is.na(METASTATIC_RECURRENCE_TIME_MONTHS) | METASTATIC_RECURRENCE_TIME_MONTHS == '' ~ 0L,
    as.numeric(METASTATIC_RECURRENCE_TIME_MONTHS) > 0 ~ 1L,
    TRUE ~ 0L
  ))
}

# --- Top 5 genes analysis ---
if (nrow(mut_summary) > 0) {
  long_genes <- mut %>% mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode), gene_name = as.character(gene_name)) %>% distinct(Tumor_Sample_Barcode, gene_name)
  gene_prev <- long_genes %>% count(gene_name, name='n_samples') %>% mutate(n_total = n_distinct(df$SAMPLE_ID), prevalence = n_samples / n_total) %>% arrange(desc(prevalence))
  top5 <- gene_prev %>% slice_head(n=5) %>% pull(gene_name)
  write_csv(gene_prev, file.path(outdir, 'gene_prevalence.csv'))
  writeLines(top5, file.path(outdir, 'top5_genes.txt'))

  for (g in top5) {
    flagcol <- paste0('G__', make.names(g))
    pres <- long_genes %>% filter(gene_name == g) %>% pull(Tumor_Sample_Barcode)
    df[[flagcol]] <- as.integer(df$SAMPLE_ID %in% pres)
  }

  df <- df %>% mutate(met_loc = case_when(
    str_to_lower(SAMPLE_TYPE) == 'metastasis' & str_detect(str_to_lower(coalesce(SAMPLE_SITE, '')), 'brain') ~ 'Brain',
    str_to_lower(SAMPLE_TYPE) == 'metastasis' ~ 'Other',
    TRUE ~ 'No'
  ))
  df$met_loc <- factor(df$met_loc, levels=c('Brain','Other','No'))

  # produce multi-panel stacked plots (PDF)
  library(ggplot2)
  pdf(file.path(outdir, 'top5_gene_metloc_stacked_panels.pdf'), width = 14, height = 6)
  for (g in top5) {
    flagcol <- paste0('G__', make.names(g))
    flag_sym <- rlang::sym(flagcol)
    tb <- df %>% group_by(!!flag_sym, met_loc) %>% summarise(n = n(), .groups='drop') %>%
      tidyr::complete(!!flag_sym := c(0,1), met_loc = c('Brain','Other','No'), fill = list(n=0)) %>%
      group_by(!!flag_sym) %>% mutate(prop = n / sum(n))
    tbp <- tb %>% filter(!!flag_sym == 1)
    if (nrow(tbp)==0) next
    gg <- ggplot(tbp, aes(x=factor(1), y=prop, fill=met_loc)) + geom_col() + coord_flip() + labs(title=g, x=NULL, y='Proportion') + theme_minimal()
    print(gg)
  }
  dev.off()
}

# --- Data dictionary creation ---
make_map <- function(varname, original_vals, mapped_vals) {
  tibble::tibble(variable = varname, original = as.character(original_vals), mapped = as.character(mapped_vals))
}

maps <- list()
if ('SAMPLE_TYPE' %in% names(df)) maps[[length(maps)+1]] <- make_map('sample_type_bin', c('Primary','Metastasis'), c('1','2'))
if ('OVERALL_TUMOR_GRADE' %in% names(df)) maps[[length(maps)+1]] <- make_map('overall_tumor_grade_ord', c('I  Low Grade (Well Differentiated)','II  Intermediate Grade (Moderately Differentiated)','III High Grade (Poorly Differentiated)','Unknown'), c('Low','Intermediate','High','NA'))
if ('RECEPTOR_STATUS_PRIMARY' %in% names(df)) maps[[length(maps)+1]] <- make_map('receptor_primary_cat', c('HR+/HER2-','HR+/HER2_Equivocal','HR+/HER2+','HR-/HER2+','HR-/HER2_Equivocal','Triple Negative'), c('1','1','2','3','3','4'))

if (length(maps)>0) {
  dict <- dplyr::bind_rows(maps)
  readr::write_csv(dict, file.path(outdir, 'data_dictionary_mappings.csv'))
}

# Save cleaned merged dataset
out_clean <- file.path(outdir, 'merged_impact_patient_level.csv')
readr::write_csv(df, out_clean)
message('Wrote cleaned dataset to: ', out_clean)
message('Done. Review output folder for derived variables, gene lists and data dictionary.')
