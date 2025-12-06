 # Analyszes + visualizes missingness then orrelation Analaysis, followed by exploratory analysis
 
  "  library(optparse)\n",
    "  library(readxl)\n",
    "  library(readr)\n",
    "  library(dplyr)\n",
    "  library(stringr)\n",
    "  library(janitor)\n",
    "  library(skimr)\n",
    "  library(purrr)\n",
    "  library(forcats)\n",
    "  library(tidyr)\n",
    "  library(tidyverse)\n",
    "  library(gtsummary)\n",
    "  library(knitr)\n",
    "  library(kableExtra)\n",
    "  library(tfrmt)\n",
    "  library(naniar)      # For missing data visualization\n",
    "  library(GGally)      # For correlation matrix\n",
    "  library(ggplot2)     # For plots\n",
    "  library(broom)       # For tidy model outputs\n",
    "  library(viridisLite) # For color palettes\n",
    "  library(scales)      # For formatting\n",
    "})\n",
    "\n",
    "# ===== Configuration and Setup =====\n",
    "message(\"📂 Working directory: \", getwd())\n",
    "\n",
    "# Create output directory\n",
    "outdir <- \"output_descriptive_vars\"\n",
    "if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)\n",
    "\n",
    "# ===== Data Loading =====\n",
    "# Load main dataset\n",
    "df_raw <- read.csv(\"~/Research/Projects/genomics_brain_mets_genie_bpc/tcga_gene_alteration_dfs_os_xgb/merged_tcga_df.csv\", na = c('NA', ''))\n",
    "colnames_original <- colnames(df_raw)\n",
    "message(\"📊 Loaded dataset with \", nrow(df_raw), \" rows and \", ncol(df_raw), \" columns\")\n",
    "\n",
    "# Load breast cancer dataset if available\n",
    "if (file.exists(\"complete_genie_bpc_dataset_deid.xlsx\")) {\n",
    "  bc_data <- read_xlsx(\"complete_genie_bpc_dataset_deid.xlsx\")\n",
    "  message(\"📊 Loaded breast cancer dataset with \", nrow(bc_data), \" rows and \", ncol(bc_data), \" columns\")\n",
    "} else {\n",
    "  message(\"⚠️ Breast cancer dataset not found, skipping BC-specific analysis\")\n",
    "  bc_data <- NULL\n",
    "}\n",
    "\n",
    
    "# ===== Missingness Analysis =====\n",
    "missing_summary <- df %>%\n",
    "  summarise(across(everything(), ~ sum(is.na(.)))) %>%\n",
    "  pivot_longer(everything(), names_to = 'variable', values_to = 'n_missing') %>%\n",
    "  mutate(n_total = nrow(df), prop_missing = n_missing / n_total)\n",
    "\n",
    "write.csv(missing_summary, file = file.path(outdir, 'missing_summary.csv'), row.names = FALSE)\n",
    "message('Wrote missing_summary.csv')\n",
    "\n",
    "# ===== Skim Summary =====\n",
    "sk <- skim(df)\n",
    "write.csv(as.data.frame(sk), file = file.path(outdir, 'skim_summary.csv'), row.names = FALSE)\n",
    "message('Wrote skim summary (csv)')\n",
    "\n",
    "# ===== Visualization Function/theme =====\n",
    "\n",
    "# Theme for plots\n",
    "theme_msk <- function(base = 10) {\n",
    "  theme_minimal(base_size = base) +\n",
    "    theme(\n",
    "      panel.grid.minor = element_blank(),\n",
    "      panel.spacing = unit(0.6, \"lines\"),\n",
    "      strip.text = element_text(face = \"bold\"),\n",
    "      legend.position = \"right\"\n",
    "    )\n",
    "}\n",
    "\n",
    "# Safe factor with explicit unknown level \".\"\n",
    "as_factor_with_dot <- function(x, lvl = NULL, ordered = FALSE) {\n",
    "  x <- as.character(x)\n",
    "  x[is.na(x) | trimws(x) == \"\"] <- \".\"\n",
    "  if (is.null(lvl)) return(factor(x, ordered = ordered))\n",
    "  factor(x, levels = c(lvl, if (!\".\" %in% lvl) \".\"), ordered = ordered)\n",
   
    # Visualize missingness\n",
    "png(file.path(outdir, 'missingness_plot.png'), width = 800, height = 600)\n",
    "print(gg_miss_var(df))\n",
    "dev.off()\n",
    "\n",
    "# Missingness pattern\n",
    "png(file.path(outdir, 'missingness_pattern.png'), width = 800, height = 600)\n",
    "print(vis_miss(df))\n",
    "dev.off()\n",
    "\n",

    # ===== Exploratory Analysis =====\n",
    "# Clean numeric subset for correlation analysis\n",
    "numeric_vars <- df %>% select(where(is.numeric))\n",
    "\n",
    "# Correlation matrix plot\n",
 def improved_scatter(df, outfile_png, outfile_pdf=None):
    group_col = '-top_5_gene_status'
    if group_col not in df.columns:
        print('Scatter plot skipped: top_5_gene_status not available in dataset.')

   
    for cand in ['SUBTYPE', 'CANCER_TYPE_DETAILED', 'CANCER_TYPE', 'TUMOR_TYPE', 'DFS_event']:
        if cand in df.columns:
            subtype_col = cand
            break
   
  
    metric_col = None
    for cand in ['TMB_NONSYNONYMOUS', 'ANEUPLOIDY_SCORE', 'BUFFA_HYPOXIA_SCORE', 'MSI_SENSOR_SCORE', 'MANTIS', 'TBL_SCORE']:
        if cand in df.columns:
            metric_col = cand
            break
    if metric_col is None:
        print('Scatter plot skipped: no numeric burden column available for plotting.')
        return

    plot_df = df[[subtype_col, metric_col, group_col]].copy()
    plot_df[group_col] = pd.to_numeric(plot_df[group_col], errors='coerce')
    plot_df[metric_col] = pd.to_numeric(plot_df[metric_col], errors='coerce')
    plot_df = plot_df.dropna(subset=[subtype_col, metric_col, group_col])
    if plot_df.empty:
        print('Scatter plot skipped: insufficient data after filtering.')
        return

    plot_df[subtype_col] = plot_df[subtype_col].astype(str)
    order = plot_df[subtype_col].value_counts().index.tolist()

    sns.set(style='whitegrid')
    plt.figure(figsize=(12, 6))
    sns.violinplot(x=subtype_col, y=metric_col, data=plot_df, order=order, inner=None, color='.9')
    sns.stripplot(x=subtype_col, y=metric_col, hue=group_col, data=plot_df, order=order, jitter=0.25, dodge=True, alpha=0.7, palette='Set1', size=4)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel(subtype_col.replace('_', ' '))
    plt.ylabel(metric_col.replace('_', ' '))
    plt.title(f"{metric_col.replace('_', ' ')} by {subtype_col.replace('_', ' ')}")
    plt.legend(title=group_col)
    plt.tight_layout()
    plt.savefig(outfile_png, dpi=300)
    if outfile_pdf:
        plt.savefig(outfile_pdf)
    plt.close()


    #demographics tables and plots

def table1_by_dfs(df, out_dir=os.path.join('output','descriptive'), group_col='top_5_gene_status'):
    """Table 1: patient characteristics and non-genetic tumor variables by DFS group and Total.
    Drops NA/Unknown/blank per covariate and writes CSV/LaTeX and a dropped-counts CSV."""
    os.makedirs(out_dir, exist_ok=True)
    label_map = {0: 'Disease Free', 1: 'Recurred or Progressed'}

    numeric_summary_vars = ['AGE','DFS_MONTHS','PFS_MONTHS','OS_MONTHS','DSS_MONTHS','DAYS_LAST_FOLLOWUP','TMB_NONSYNONYMOUS','ANEUPLOIDY_SCORE','MSI_SENSOR_SCORE','BUFFA_HYPOXIA_SCORE','MANTIS','TBL_SCORE']
    categorical_summary_vars = ['AGE_CAT','SEX','RACE','ETHNICITY','SUBTYPE','CANCER_TYPE','CANCER_TYPE_DETAILED','AJCC_PATHOLOGIC_TUMOR_STAGE','MANTIS_BIN','HISTORY_NEOADJUVANT_TRTYN','SOMATIC_STATUS','SAMPLE_TYPE','TUMOR_TYPE','TUMOR_TISSUE_SITE','TISSUE_SOURCE_SITE','TISSUE_SOURCE_SITE_CODE']

    # First, collapse to patient-level using merge_key or patient_ID if present
    id_col = None
    for cand in ['merge_key', 'patient_ID', 'patient_id']:
        if cand in df.columns:
            id_col = cand
            break
    if id_col is None:
        patient_df = df.copy()
        print('Warning: no merge_key/patient_ID found; Table1 will run on original rows (variant-level).')
    else:
        patient_cols = list(dict.fromkeys(numeric_summary_vars + categorical_summary_vars + [group_col]))

        def _agg_numeric(series):
            s = pd.to_numeric(series, errors='coerce')
            return float(s.mean()) if s.notna().any() else pd.NA

        def _agg_first(series):
            s = series.dropna()
            return s.iloc[0] if not s.empty else pd.NA

        agg_map = {}
        for c in patient_cols:
            if c not in df.columns:
                continue
            if c in numeric_summary_vars:
                agg_map[c] = _agg_numeric
            elif c == group_col:
                agg_map[c] = _agg_first
            else:
                agg_map[c] = _agg_first

        if agg_map:
            patient_df = df.groupby(id_col).agg(agg_map).reset_index()
        else:
            patient_df = df.groupby(id_col).size().reset_index(name='count')

        patient_df[group_col] = pd.to_numeric(patient_df.get(group_col), errors='coerce')
        patient_df[group_col] = patient_df[group_col].where(patient_df[group_col].isin([0, 1]), pd.NA)

    working_df = patient_df
    total_n = len(working_df)

    rows = []
    dropped = []
    for v in numeric_summary_vars + categorical_summary_vars:
        if v not in working_df.columns:
            dropped.append({'variable': v, 'n_dropped_total': total_n, 'n_dropped_group_0': total_n, 'n_dropped_group_1': total_n})
            continue

        ser = working_df[v]
        mask_drop = _drop_mask_for_ser(working_df, ser)
        n_drop_total = int(mask_drop.sum())
        if group_col in working_df.columns:
            n_drop_g0 = int(mask_drop[working_df[group_col] == 0].sum())
            n_drop_g1 = int(mask_drop[working_df[group_col] == 1].sum())
        else:
            n_drop_g0 = 0
            n_drop_g1 = 0
        dropped.append({'variable': v, 'n_dropped_total': n_drop_total, 'n_dropped_group_0': n_drop_g0, 'n_dropped_group_1': n_drop_g1})

        ser_clean = ser[~mask_drop]
        if v in numeric_summary_vars:
            for g in [0, 1]:
                if group_col not in working_df.columns:
                    val = ''
                else:
                    vals = pd.to_numeric(working_df.loc[(working_df[group_col] == g) & (~mask_drop), v], errors='coerce').dropna()
                    if vals.empty:
                        val = ''
                    else:
                        med = vals.median()
                        iqr = vals.quantile(0.75) - vals.quantile(0.25)
                        val = f"{med:.2f} (IQR {iqr:.2f})"
                rows.append({'variable': v, 'level': '', 'group': label_map.get(g), 'value': val})
            vals_all = pd.to_numeric(ser_clean, errors='coerce').dropna()
            if vals_all.empty:
                allval = ''
            else:
                med = vals_all.median()
                iqr = vals_all.quantile(0.75) - vals_all.quantile(0.25)
                allval = f"{med:.2f} (IQR {iqr:.2f})"
            rows.append({'variable': v, 'level': '', 'group': 'Total', 'value': allval})
        else:
            if ser_clean.empty:
                continue
            levels = list(pd.Series(ser_clean.astype(str)).value_counts().index)
            for lev in levels:
                for g in [0, 1]:
                    if group_col not in working_df.columns:
                        cnt = 0
                        pct = 0.0
                    else:
                        sub = working_df.loc[(working_df[group_col] == g) & (~mask_drop), v].dropna().astype(str)
                        cnt = int((sub == str(lev)).sum())
                        pct = cnt / len(sub) * 100 if len(sub) > 0 else 0.0
                    rows.append({'variable': v, 'level': str(lev), 'group': label_map.get(g), 'value': f"{cnt} ({pct:.1f}%)"})
                cnt_all = int((ser_clean.astype(str) == str(lev)).sum())
                pct_all = cnt_all / len(ser_clean) * 100 if len(ser_clean) > 0 else 0.0
                rows.append({'variable': v, 'level': str(lev), 'group': 'Total', 'value': f"{cnt_all} ({pct_all:.1f}%)"})

    out = pd.DataFrame(rows)
    dropped_df = pd.DataFrame(dropped)
    dropped_df.to_csv(os.path.join(out_dir, 'table1_dropped_counts.csv'), index=False)

    if out.empty:
        print('Table1: no non-genetic variables available after dropping NA/Unknown/blank.')
        return out, dropped_df

    pivot = out.pivot_table(index=['variable', 'level'], columns='group', values='value', aggfunc=lambda x: ' | '.join(x.astype(str))).reset_index()
    for col in [label_map[0], label_map[1], 'Total']:
        if col not in pivot.columns:
            pivot[col] = ''
    pivot.to_csv(os.path.join(out_dir, 'table1_by_dfs.csv'), index=False)

    n0 = int(working_df[working_df[group_col] == 0].shape[0]) if group_col in working_df.columns else 0
    n1 = int(working_df[working_df[group_col] == 1].shape[0]) if group_col in working_df.columns else 0
    ntot = int(working_df.shape[0])
    tex_lines = []
    tex_lines.append('\\begin{table}[htb]')
    tex_lines.append('\\centering')
    tex_lines.append('\\small')
    tex_lines.append('\\begin{tabular}{llccc}')
    tex_lines.append('\\hline')
    header = "Variable & Level & Disease Free (n={}) & Recurred/Progressed (n={}) & Total (n={})".format(n0, n1, ntot)
    tex_lines.append(header + ' \\\\')
    for _, r in pivot.iterrows():
        var = r['variable']
        lev = r['level'] if r['level'] and str(r['level']) != 'nan' else ''
        a = r.get(label_map[0], '')
        b = r.get(label_map[1], '')
        t = r.get('Total', '')
        tex_lines.append(f"{var} & {lev} & {a} & {b} & {t} " + ' \\\\')
    tex_lines.append('\\hline')
    tex_lines.append('\\end{tabular}')
    tex_lines.append('\\caption{Patient characteristics by DFS status. Categorical variables are reported as n (\\%), continuous variables as median (IQR).}')
    tex_lines.append('\\end{table}')
    with open(os.path.join(out_dir, 'table1_by_dfs.tex'), 'w') as fh:
        fh.write('\n'.join(tex_lines))

    return pivot, dropped_df
    "# ===== Genomic Visualizations =====\n",
    
# 8) tmb quartiles (TMB_NONSYNONYMOUS)


# 9) FGA (Fraction Genome Altered) quartiles

# 10) MSI Type unordered categorical


# 11) mutation_count_sum -> mutation_count_all_sites_sum (ordered quartiles)


# 12) tmb_nonsynonmous


# 13) t_alt_count quartiles from mutation summary

    # Variant Allele Count Analysis\n",
    "if (!is.null(var_acol) && var_acol %in% names(df)) {\n",
    # Histogram of variant allele count\n",
    "  p_hist_vac <- ggplot(df, aes(x = .data[[var_acol]])) +\n",
    "    geom_histogram(bins = 30, alpha = 0.8, fill = \"skyblue\") +\n",
    "    labs(x = \"Variant Allele Count\", y = \"Frequency\", title = \"Distribution of Variant Allele Count\") +\n",
    "    theme_msk()\n",
    "  \n",
    "  ggsave(file.path(outdir, 'plot_vac_hist.png'), p_hist_vac, width = 8, height = 6, dpi = 300)\n",
    "  \n",
    # Boxplot by DFS status if available\n",
    "  if ('dfs_f' %in% names(df)) {\n",
    "    p_box_vac <- ggplot(df, aes(x = dfs_f, y = .data[[var_acol]], fill = dfs_f)) +\n",
    "      geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +\n",
    "      labs(x = \"DFS Status\", y = \"Variant Allele Count\", title = \"Variant Allele Count by DFS Status\") +\n",
    "      theme_msk()\n",
    "    \n",
    "    ggsave(file.path(outdir, 'plot_vac_by_dfs.png'), p_box_vac, width = 8, height = 6, dpi = 300)\n",
    "  }\n",
    "}\n",
    "\n",
"# Gene prevalence plots\n",
    "if (exists('gene_prev') && nrow(gene_prev) > 0 && 'dfs_f' %in% names(df)) {\n",
    "  # Create long format for gene analysis\n",
    "  if (!is.null(hugo_col) && exists('long_genes')) {\n",
    "    gene_dfs_data <- long_genes %>%\n",
    "      left_join(df %>% select(merge_key = row_number(), dfs_f), by = \"merge_key\") %>%\n",
    "      group_by(gene, dfs_f) %>%\n",
    "      summarise(n_with_gene = n(), .groups = 'drop') %>%\n",
    "      left_join(\n",
    "        df %>% group_by(dfs_f) %>% summarise(total_n = n(), .groups = 'drop'),\n",
    "        by = \"dfs_f\"\n",
    "      ) %>%\n",
    "      mutate(prevalence = n_with_gene / total_n) %>%\n",
    "      filter(gene %in% top10_genes)\n",
    "    \n",
    "    # Plot by DFS status\n",
    "    for (dfs_level in unique(gene_dfs_data$dfs_f)) {\n",
    "      if (is.na(dfs_level)) next\n",
    "      \n",
    "      p <- gene_dfs_data %>%\n",
    "        filter(dfs_f == dfs_level) %>%\n",
    "        arrange(prevalence) %>%\n",
    "        mutate(gene = factor(gene, levels = gene)) %>%\n",
    "        ggplot(aes(x = gene, y = prevalence)) +\n",
    "        geom_segment(aes(xend = gene, yend = 0), color = \"gray50\") +\n",
    "        geom_point(size = 3, color = \"orange\") +\n",
    "        coord_flip() +\n",
    "        labs(x = NULL, y = \"Prevalence\", title = paste(\"Gene Mutation Prevalence —\", dfs_level)) +\n",
    "        theme_msk()\n",
    "      \n",
    "      ggsave(file.path(outdir, paste0('plot_gene_prev_', gsub(\"[^A-Za-z0-9]\", \"_\", dfs_level), '.png')), \n",
    "             p, width = 8, height = 6, dpi = 300)\n",
    "    }\n",
    "  }\n",
    "}\n",
    "\n",


def table2_by_dfs(df, out_dir=os.path.join('output', 'descriptive'), group_col='DFS_event', top_genes=20):
    ""Table 2:Molecular burden and outcome metrics by DFS group and Total."""
    os.makedirs(out_dir, exist_ok=True)
    label_map = {0: 'Disease Free', 1: Died or Recurred'}

    numeric_vars = ['TMB_NONSYNONYMOUS', 'ANEUPLOIDY_SCORE', 'MSI_SENSOR_SCORE', 'MANTIS', 'BUFFA_HYPOXIA_SCORE', 'TBL_SCORE', 'DFS_MONTHS', 'PFS_MONTHS', 'OS_MONTHS', 'DSS_MONTHS', 'DAYS_LAST_FOLLOWUP']
    categorical_vars = ['MANTIS_BIN', 'HISTORY_NEOADJUVANT_TRTYN', 'SOMATIC_STATUS', 'SAMPLE_TYPE', 'TUMOR_TYPE', 'TUMOR_TISSUE_SITE', 'TISSUE_SOURCE_SITE', 'TISSUE_SOURCE_SITE_CODE', 'OS_STATUS', 'PFS_STATUS', 'DSS_STATUS']

    df[group_col] = pd.to_numeric(df.get(group_col), errors='coerce')
    df[group_col] = df[group_col].where(df[group_col].isin([0, 1]), pd.NA)

    pid_col = None
    for cand in ['merge_key', 'PATIENT_ID', 'patient_ID', 'patient_id']:
        if cand in df.columns:
            pid_col = cand
            break
    if pid_col is None:
        print('Warning: no patient identifier found; Table2 will run on original rows (variant-level).')
        patient_df = df.copy()
    else:
        patient_cols = list(dict.fromkeys(numeric_vars + categorical_vars + [group_col]))

        def _agg_numeric(series):
            s = pd.to_numeric(series, errors='coerce')
            return float(s.mean()) if s.notna().any() else pd.NA

        def _agg_first(series):
            s = series.dropna()
            return s.iloc[0] if not s.empty else pd.NA

        agg_map = {}
        for c in patient_cols:
            if c not in df.columns:
                continue
            if c in numeric_vars:
                agg_map[c] = _agg_numeric
            elif c == group_col:
                agg_map[c] = _agg_first
            else:
                agg_map[c] = _agg_first

        if agg_map:
            patient_df = df.groupby(pid_col).agg(agg_map).reset_index()
        else:
            patient_df = df.groupby(pid_col).size().reset_index(name='count')

    total_n = len(patient_df)
    rows = []
    dropped = []

    available_numeric = [v for v in numeric_vars if v in patient_df.columns]
    available_categorical = [v for v in categorical_vars if v in patient_df.columns]

    for v in available_numeric + available_categorical:
        ser = patient_df[v]
        mask_drop = _drop_mask_for_ser(patient_df, ser)
        n_drop_total = int(mask_drop.sum())
        if group_col in patient_df.columns:
            n_drop_g0 = int(mask_drop[patient_df[group_col] == 0].sum())
            n_drop_g1 = int(mask_drop[patient_df[group_col] == 1].sum())
        else:
            n_drop_g0 = 0
            n_drop_g1 = 0
        dropped.append({'variable': v, 'n_dropped_total': n_drop_total, 'n_dropped_group_0': n_drop_g0, 'n_dropped_group_1': n_drop_g1})

        ser_clean = ser[~mask_drop]
        if v in available_numeric:
            for g in [0, 1]:
                if group_col not in patient_df.columns:
                    val = ''
                else:
                    vals = pd.to_numeric(patient_df.loc[(patient_df[group_col] == g) & (~mask_drop), v], errors='coerce').dropna()
                    if vals.empty:
                        val = ''
                    else:
                        med = vals.median()
                        iqr = vals.quantile(0.75) - vals.quantile(0.25)
                        val = f"{med:.2f} (IQR {iqr:.2f})"
                rows.append({'variable': v, 'level': '', 'group': label_map.get(g), 'value': val})
            vals_all = pd.to_numeric(ser_clean, errors='coerce').dropna()
            if vals_all.empty:
                allval = ''
            else:
                med = vals_all.median()
                iqr = vals_all.quantile(0.75) - vals_all.quantile(0.25)
                allval = f"{med:.2f} (IQR {iqr:.2f})"
            rows.append({'variable': v, 'level': '', 'group': 'Total', 'value': allval})
        else:
            if ser_clean.empty:
                continue
            levels = list(pd.Series(ser_clean.astype(str)).value_counts().index)
            for lev in levels:
                for g in [0, 1]:
                    if group_col not in patient_df.columns:
                        cnt = 0
                        pct = 0.0
                    else:
                        sub = patient_df.loc[(patient_df[group_col] == g) & (~mask_drop), v].dropna().astype(str)
                        cnt = int((sub == str(lev)).sum())
                        pct = cnt / len(sub) * 100 if len(sub) > 0 else 0.0
                    rows.append({'variable': v, 'level': str(lev), 'group': label_map.get(g), 'value': f"{cnt} ({pct:.1f}%)"})
                cnt_all = int((ser_clean.astype(str) == str(lev)).sum())
                pct_all = cnt_all / len(ser_clean) * 100 if len(ser_clean) > 0 else 0.0
                rows.append({'variable': v, 'level': str(lev), 'group': 'Total', 'value': f"{cnt_all} ({pct_all:.1f}%)"})

    out = pd.DataFrame(rows)
    dropped_df = pd.DataFrame(dropped)
    dropped_df.to_csv(os.path.join(out_dir, 'table2_dropped_counts.csv'), index=False)

    if out.empty:
        print('Table2: no molecular/outcome variables available after dropping NA/Unknown/blank.')
        return out, dropped_df

    pivot = out.pivot_table(index=['variable', 'level'], columns='group', values='value', aggfunc=lambda x: ' | '.join(x.astype(str))).reset_index()
    for col in [label_map[0], label_map[1], 'Total']:
        if col not in pivot.columns:
            pivot[col] = ''
    pivot.to_csv(os.path.join(out_dir, 'table2_by_dfs.csv'), index=False)

    n0 = int(patient_df[patient_df[group_col] == 0].shape[0]) if group_col in patient_df.columns else 0
    n1 = int(patient_df[patient_df[group_col] == 1].shape[0]) if group_col in patient_df.columns else 0
    ntot = int(patient_df.shape[0])
    tex_lines = []
    tex_lines.append('\\begin{table}[htb]')
    tex_lines.append('\\centering')
    tex_lines.append('\\small')
    tex_lines.append('\\begin{tabular}{llccc}')
    tex_lines.append('\\hline')
    header = "Variable & Level & Disease Free (n={}) & Recurred/Death (n={}) & Total (n={})".format(n0, n1, ntot)
    tex_lines.append(header + ' \\\\')
    for _, r in pivot.iterrows():
        var = r['variable']
        lev = r['level'] if r['level'] and str(r['level']) != 'nan' else ''
        a = r.get(label_map[0], '')
        b = r.get(label_map[1], '')
        t = r.get('Total', '')
        tex_lines.append(f"{var} & {lev} & {a} & {b} & {t} " + ' \\\\')
    tex_lines.append('\\hline')
    tex_lines.append('\\end{tabular}')
    tex_lines.append('\\caption{Molecular burden and outcome metrics by DFS status. Continuous variables are reported as median (IQR); categorical variables are reported as n (\\%).}')
    tex_lines.append('\\end{table}')
    with open(os.path.join(out_dir, 'table2_by_dfs.tex'), 'w') as fh:
        fh.write('\n'.join(tex_lines))

    return pivot, dropped_df

#proportion of variablle amongst outcome group
    "plot_grouped_prop <- function(data, group_var, out_file) {\n",
    "  if (!'dfs_bin' %in% names(data) || !group_var %in% names(data)) return()\n",
    "  \n",
    "  p <- data %>%\n",
    "    group_by(.data[[group_var]]) %>%\n",
    "    summarise(\n",
    "      n = n(),\n",
    "      prop_dfs1 = mean(dfs_bin == 1, na.rm = TRUE),\n",
    "      .groups = 'drop'\n",
    "    ) %>%\n",
    "    filter(!is.na(.data[[group_var]])) %>%\n",
    "    mutate(group_factor = fct_reorder(as.factor(.data[[group_var]]), prop_dfs1, .desc = TRUE)) %>%\n",
    "    ggplot(aes(x = group_factor, y = prop_dfs1)) +\n",
    "    geom_col(fill = \"steelblue\", alpha = 0.7) +\n",
    "    geom_text(aes(label = percent(prop_dfs1, accuracy = 0.1)), vjust = -0.3, size = 3) +\n",
    "    labs(x = NULL, y = \"P(DFS = 1)\", title = paste(\"Probability of DFS = 1 by\", group_var)) +\n",
    "    scale_y_continuous(labels = percent_format()) +\n",
    "    theme_msk() +\n",
    "    theme(axis.text.x = element_text(angle = 45, hjust = 1))\n",
    "  \n",
    "  ggsave(file.path(outdir, out_file), p, width = 8, height = 6, dpi = 300)\n",
    "}\n",
    "\n",
    "# Create demographic plots if variables exist\n",
    "if ('dfs_bin' %in% names(df)) {\n",
    "  demo_vars <- intersect(names(df), c('race', 'ethnicity', 'subtype', 'sex'))\n",
    "  for (var in demo_vars) {\n",
    "    plot_grouped_prop(df, var, paste0('plot_prob_dfs1_by_', var, '.png'))\n",
    "  }\n",
    "}\n",
    "\n",
    "# ===== Statistical Analysis =====\n",
    "\n",
    "if (!is.null(dfs_status_col) && 'top_5_gene_status' %in% names(df)) {\n",
    "  \n",
    "  # Univariate linear regression\n",
    "  fit_uni <- lm(dfs_bin ~ top_5_gene_status, data = df)\n",
    "  summary(fit_uni)\n",
    "  tidy(fit_uni) %>% write.csv(file.path(outdir, 'univariate_lm.csv'), row.names = FALSE)\n",
    "  \n",
    "  # Multivariate regression with possible confounders\n",
    "  confounders <- intersect(names(df), c('sex_bin', 'age_cat', 'ajcc_pathologic_tumor_stage', \n",
    "                                        'buffa_hypoxia_quartile', 'aneuploidy_score_quartile'))\n",
    "  \n",
    "  if (length(confounders) > 0) {\n",
    "    formula_str <- paste(\"dfs_bin ~ top_5_gene_status +\", paste(confounders, collapse = \" + \"))\n",
    "    fit_multi <- lm(as.formula(formula_str), data = df)\n",
    "    summary(fit_multi)\n",
    "    tidy(fit_multi) %>% write.csv(file.path(outdir, 'multivariable_lm.csv'), row.names = FALSE)\n",
    "    \n",
    "    # Residual diagnostics\n",
    "    png(file.path(outdir, 'residuals_diagnostics.png'), width = 800, height = 800)\n",
    "    par(mfrow = c(2, 2))\n",
    "    plot(fit_multi)\n",
    "    dev.off()\n",
    "  }\n",
    "}\n",
    "\n",