#!/usr/bin/env Rscript
## Diagnose which predictors trigger glm warnings (separation/convergence)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(rlang)
  library(stringr)
})

outdir <- file.path("C:", "Users", "jamesr4", "OneDrive - Memorial Sloan Kettering Cancer Center", "Documents", "Research", "Projects", "impact", "Output")
infile <- file.path(outdir, 'merged_impact_patient_level.csv')
if (!file.exists(infile)) stop('Merged input not found: ', infile)
df <- read_csv(infile, show_col_types = FALSE)

# ensure .dfs exists
if (!('.dfs' %in% names(df))){
  if ('DFS_EVENT' %in% names(df)) df <- df %>% mutate(.dfs = as.integer(DFS_EVENT)) else if ('DFS_EVENT_f' %in% names(df)) df <- df %>% mutate(.dfs = as.integer(DFS_EVENT_f == 'Event')) else stop('No DFS_EVENT found; run generate_tables_1_2.R first')
}

# predictors same logic as generate script
dem_vars <- c('age_cat','Race_clean','Ethnicity_clean','sample_site','sample_type_bin')
classical <- c('stage_diag_group','overall_tumor_grade_ord','receptor_primary_cat')
possible_mut_vars <- c('tmb_quartile','fga_quartile','MSI_type','mutation_count_q','tmb_nonsynonmous_q','t_alt_count_q','mutation_count_all_sites_sum','tmb_val','t_alt_count')
gene_flags <- grep('^G__', names(df), value = TRUE)
mut_vars_present <- intersect(possible_mut_vars, names(df))
non_classical <- c(mut_vars_present, gene_flags)

preds <- unique(c(classical, non_classical))
preds <- preds[preds %in% names(df)]

diag_rows <- list()
warnings_all <- list()

for (v in preds){
  vec <- df[[v]]
  n_total <- length(vec)
  n_na <- sum(is.na(vec))
  n_event <- sum(df$.dfs==1, na.rm = TRUE)
  n_non <- sum(df$.dfs==0, na.rm = TRUE)

  # level summary
  if (is.numeric(vec) || is.integer(vec)){
    type <- 'numeric'
    med_all <- tryCatch(median(as.numeric(vec), na.rm=TRUE), error=function(e) NA_real_)
    mean_all <- tryCatch(mean(as.numeric(vec), na.rm=TRUE), error=function(e) NA_real_)
    sd_all <- tryCatch(sd(as.numeric(vec), na.rm=TRUE), error=function(e) NA_real_)
    levels_info <- paste0('mean=', signif(mean_all,4), '; sd=', signif(sd_all,4), '; median=', signif(med_all,4))
  } else {
    type <- 'categorical'
    tabs <- table(as.character(vec), useNA='ifany')
    levels_info <- paste(names(tabs), tabs, sep=':', collapse='; ')
  }

  warning_msgs <- character()
  error_msg <- NA_character_

  # fit univariate glm and capture warnings/errors
  formula_text <- paste('.dfs ~', paste0('`', v, '`'))
  form <- as.formula(formula_text)
  wcollector <- character()
  res <- withCallingHandlers({
    m <- tryCatch(glm(form, data = df %>% select(.dfs, all_of(v)) %>% filter(!is.na(.dfs)), family=binomial()), error=function(e) e)
    m
  }, warning = function(w){
    wcollector <<- c(wcollector, conditionMessage(w))
    invokeRestart('muffleWarning')
  })

  if (inherits(res, 'error')){
    error_msg <- conditionMessage(res)
  }
  if (length(wcollector)>0) warning_msgs <- unique(wcollector)
  warnings_all[[v]] <- warning_msgs

  diag_rows[[length(diag_rows)+1]] <- tibble::tibble(variable = v, type = type, n_total = n_total, n_na = n_na, n_event = n_event, n_non_event = n_non, levels_info = levels_info, warnings = paste(warning_msgs, collapse=' | '), error = ifelse(is.na(error_msg),'', error_msg))
}

diag_df <- dplyr::bind_rows(diag_rows)
readr::write_csv(diag_df, file.path(outdir, 'glm_diagnostics.csv'))
if (length(warnings_all)>0) writeLines(unlist(lapply(names(warnings_all), function(n) paste0(n, ': ', paste(warnings_all[[n]], collapse=' | ')))), con = file.path(outdir, 'glm_diagnostics_warnings.txt'))

message('Wrote glm_diagnostics.csv and glm_diagnostics_warnings.txt to ', outdir)
