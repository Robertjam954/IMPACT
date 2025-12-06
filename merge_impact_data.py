"""Merge MSK IMPACT patient/sample data sheets and mutation data into a single CSV/Excel file.

Usage:
  python scripts/merge_genie_sheets.py --patient-file path/to/patient_data.xlsx --mutation-file path/to/mutation_data.xlsx --key patient_id

Performs left joins to keep all patients from the base dataset.
If --key is not provided the script will attempt to find a common key column.
"""
import argparse
import os
import re
import pandas as pd
import shutil
import tempfile


def list_sheets(path):
    """Return sheet names for the given workbook path."""
    safe = _ensure_readable_copy(path)
    xls = pd.ExcelFile(safe)
    return xls.sheet_names


def _ensure_readable_copy(path):
    """Return a path that can be safely opened by pandas.

    If the original file cannot be opened due to permission issues (OneDrive lock
    or Excel open), copy it to a temp file and return the temp path.
    """
    try:
        # quick test open
        with open(path, 'rb'):
            return path
    except PermissionError:
        tmpdir = tempfile.gettempdir()
        base = os.path.basename(path)
        tmp = os.path.join(tmpdir, f"copy_{base}")
        shutil.copy2(path, tmp)
        return tmp


def preview_sheet(path, sheet, n=5):
    safe = _ensure_readable_copy(path)
    df = pd.read_excel(safe, sheet_name=sheet)
    return df.head(n), df.columns.tolist()


def guess_common_key(dfs):
    # Return the column name that is common to all dataframes (case-insensitive)
    sets = [set([c.lower() for c in df.columns]) for df in dfs]
    common = set.intersection(*sets)
    # Prefer columns with 'id' or 'sample' keywords
    for pref in ["id", "sample", "patient", "case"]:
        for c in common:
            if pref in c:
                return c
    return next(iter(common)) if common else None


def main():
    parser = argparse.ArgumentParser()
    default_patient_file = r"C:\Users\jamesr4\OneDrive - Memorial Sloan Kettering Cancer Center\Documents\Research\Projects\impact\Data\msk impact patient and sample data.xlsx"
    default_mutation_file = r"C:\Users\jamesr4\OneDrive - Memorial Sloan Kettering Cancer Center\Documents\Research\Projects\impact\Data\msk impact mutation data.xlsx"
    
    parser.add_argument("--patient-file", required=False, default=default_patient_file, help="Path to patient and sample data Excel file")
    parser.add_argument("--mutation-file", required=False, default=default_mutation_file, help="Path to mutation data Excel file")
    parser.add_argument("--key", default=None, help="Merge key column name (case-insensitive), e.g. 'patient_id'")
    parser.add_argument("--out", default="output/merged_impact.csv", help="Output file path")
    args = parser.parse_args()

    patient_file = args.patient_file
    mutation_file = args.mutation_file

    if not os.path.exists(patient_file):
        raise FileNotFoundError(f"Patient file not found: {patient_file}")
    if not os.path.exists(mutation_file):
        raise FileNotFoundError(f"Mutation file not found: {mutation_file}")

    print("Patient file sheets:", list_sheets(patient_file))
    print("Mutation file sheets:", list_sheets(mutation_file))

    # Read only the first two sheets from patient file (user requested two-sheet merge)
    print(f"Reading patient file: {patient_file}")
    patient_sheets = list_sheets(patient_file)
    if len(patient_sheets) < 2:
        raise ValueError(f"Expected at least 2 sheets in patient file, found {len(patient_sheets)}")
    patient_sheets = patient_sheets[:2]
    dfs = []
    for s in patient_sheets:
        print(f"Reading sheet '{s}' from patient file...")
        df = pd.read_excel(patient_file, sheet_name=s)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {df.columns.tolist()[:10]}")
        dfs.append(df)
    
    # Read mutation file (assume first/only sheet)
    print(f"Reading mutation file: {mutation_file}")
    mutation_sheets = list_sheets(mutation_file)
    mutation_sheet = mutation_sheets[0]  # Use first sheet
    print(f"Reading sheet '{mutation_sheet}' from mutation file...")
    mutation_df = pd.read_excel(mutation_file, sheet_name=mutation_sheet)
    print(f"  Shape: {mutation_df.shape}")
    print(f"  Columns: {mutation_df.columns.tolist()[:10]}")

    key = args.key
    if key is None:
        guessed = guess_common_key(dfs)
        if guessed:
            print("Guessed common key (lowercased):", guessed)
            # find actual column name matching guessed key (original case)
            for c in dfs[0].columns:
                if c.lower() == guessed:
                    key = c
                    break
        else:
            raise ValueError("Could not guess a common key. Provide --key explicitly.")

    print("Using merge key:", key)

    # Step 1: Merge patient file sheets using left joins
    print(f"\nStep 1: Merging {len(dfs)} sheets from patient file...")
    merged = dfs[0].rename(columns={c: c.strip() for c in dfs[0].columns})
    def _norm_col(s):
        if s is None:
            return ''
        return re.sub(r"[^a-z0-9]", "", str(s).lower())

    # Find key column in first dataframe using normalized names
    key_col = None
    target_norm = _norm_col(key)
    for c in merged.columns:
        if _norm_col(c) == target_norm:
            key_col = c
            break
    if key_col is None:
        raise ValueError(f"Key '{key}' not found in first sheet")
    
    merged = merged.rename(columns={key_col: "merge_key"})
    print(f"Using key column: {key_col} -> merge_key")
    print(f"Base dataset shape: {merged.shape}")
    
    # Merge remaining patient sheets with left joins
    for i, df in enumerate(dfs[1:], 1):
        df2 = df.rename(columns={c: c.strip() for c in df.columns})
        # find matching column name using normalized names (handles spaces/underscores/case)
        match = None
        for c in df2.columns:
            if _norm_col(c) == target_norm:
                match = c
                break
        if match is None:
            print(f"Warning: Key '{key}' not found in sheet {i+1}, skipping...")
            continue
        
        df2 = df2.rename(columns={match: "merge_key"})
        before_shape = merged.shape
        merged = pd.merge(merged, df2, on="merge_key", how="left", suffixes=("", "_r"))
        print(f"After merging sheet {i+1}: {before_shape} -> {merged.shape}")
    
    # Step 2: Merge with mutation data using left join
    print(f"\nStep 2: Merging with mutation data...")
    mutation_df_clean = mutation_df.rename(columns={c: c.strip() for c in mutation_df.columns})
    
    # Find key column in mutation data using normalized names
    mutation_key = None
    for c in mutation_df_clean.columns:
        if _norm_col(c) == target_norm:
            mutation_key = c
            break
    
    if mutation_key is None:
        print(f"Warning: Key '{key}' not found in mutation data, skipping mutation merge...")
    else:
        mutation_df_clean = mutation_df_clean.rename(columns={mutation_key: "merge_key"})
        before_shape = merged.shape
        merged = pd.merge(merged, mutation_df_clean, on="merge_key", how="left", suffixes=("", "_mut"))
        print(f"After merging mutation data: {before_shape} -> {merged.shape}")

    # save variant-level merged output (current behavior)
    out = args.out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    merged.to_csv(out, index=False)
    try:
        merged.to_excel(out.replace('.csv', '.xlsx'), index=False)
    except Exception:
        pass

    # Keep a copy of the patient-level merged base (before mutation join)
    base_merged = merged.copy()

    # helper to find column by normalized keywords
    def _find_col(df, keywords):
        def norm(s):
            return re.sub(r"[^a-z0-9]", "", str(s).lower())
        for c in df.columns:
            nc = norm(c)
            for kw in keywords:
                if kw in nc:
                    return c
        return None

    # determine sample and site columns from base_merged
    sample_id_col = _find_col(base_merged, ['sampleid', 'sample_id', 'sample'])
    sample_site_col = _find_col(base_merged, ['sample_site', 'site', 'metastaticsite', 'tumor_site'])
    sample_type_col = _find_col(base_merged, ['sampletype', 'sample_type', 'sampleclass'])

    # create sample-level table (one row per patient-sample) preserving locations
    samples_df = None
    sample_cols = ['merge_key']
    if sample_id_col and sample_id_col in base_merged.columns:
        sample_cols.append(sample_id_col)
    if sample_type_col and sample_type_col in base_merged.columns:
        sample_cols.append(sample_type_col)
    if sample_site_col and sample_site_col in base_merged.columns:
        sample_cols.append(sample_site_col)

    if len(sample_cols) > 1:
        samples_df = base_merged.loc[:, [c for c in sample_cols if c in base_merged.columns]].drop_duplicates()
        samples_out = os.path.join(os.path.dirname(out), 'samples_per_patient.csv')
        samples_df.to_csv(samples_out, index=False)
        print('Wrote sample-level per-patient file:', samples_out)

    # aggregate mutation metrics from mutation_df_clean if present
    mut_agg = None
    if 'mutation_df_clean' in locals() and 'merge_key' in mutation_df_clean.columns:
        md = mutation_df_clean
        aggs = {}
        # detect mutation count column if present
        mut_count_col = _find_col(md, ['mutationcount', 'mutation_count', 'mutationcount'])
        if mut_count_col and mut_count_col in md.columns:
            aggs['mutation_count_sum'] = md.groupby('merge_key')[mut_count_col].sum()
        # count mutation rows per patient
        aggs['mutation_row_count'] = md.groupby('merge_key').size()
        # detect TMB column
        tmb_col = _find_col(md, ['tmb', 'tmbnonsynonymous', 'tmb_nonsynonymous'])
        if tmb_col and tmb_col in md.columns:
            aggs['tmb_mean'] = md.groupby('merge_key')[tmb_col].mean()

        if aggs:
            mut_agg = pd.concat(aggs, axis=1).reset_index()

    # Build patient-level table from base_merged: pick representative patient-level values
    patient_level_cols = [c for c in base_merged.columns if c != 'merge_key' and c not in (sample_id_col or [])]

    def _agg_patient(df):
        rows = []
        for k, g in df.groupby('merge_key'):
            row = {'merge_key': k}
            for c in patient_level_cols:
                ser = g[c].dropna()
                if ser.empty:
                    row[c] = pd.NA
                else:
                    if pd.api.types.is_numeric_dtype(ser):
                        row[c] = float(ser.median())
                    else:
                        # keep first observed non-null (could be sample-specific; sample details preserved separately)
                        row[c] = str(ser.iloc[0])
            rows.append(row)
        return pd.DataFrame(rows)

    patient_df = _agg_patient(base_merged)

    # attach sample summary (concatenate sample ids and sites per patient)
    if samples_df is not None:
        def _join_unique(s):
            vals = sorted(set([str(x) for x in s.dropna().astype(str)]))
            return ';'.join(vals) if vals else pd.NA

        samples_grp = samples_df.groupby('merge_key').agg(
            sample_ids=(sample_id_col, _join_unique) if sample_id_col else (lambda s: pd.NA),
            sample_sites=(sample_site_col, _join_unique) if sample_site_col else (lambda s: pd.NA),
            sample_types=(sample_type_col, _join_unique) if sample_type_col else (lambda s: pd.NA),
        ).reset_index()
        patient_df = patient_df.merge(samples_grp, on='merge_key', how='left')

    # attach mutation aggregates
    if mut_agg is not None and not mut_agg.empty:
        patient_df = patient_df.merge(mut_agg, on='merge_key', how='left')

    patient_out = os.path.join(os.path.dirname(out), 'merged_impact_patient_level.csv')
    patient_df.to_csv(patient_out, index=False)
    try:
        patient_df.to_excel(patient_out.replace('.csv', '.xlsx'), index=False)
    except Exception:
        pass
    print('Wrote patient-level aggregated file:', patient_out)

    print('Wrote:', out)


if __name__ == '__main__':
    main()
