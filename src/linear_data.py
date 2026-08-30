import numpy as np
import pandas as pd
from pathlib import Path
import re



def read_linear_pep_data(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    
    if path.suffix.lower() != ".tsv":
            raise ValueError("tsv file needed")
    if not path.is_file():
        raise FileNotFoundError(f"no file")

    data_table = pd.read_csv(path, sep="\t")

    if data_table.empty:
        raise ValueError(".tsv file is empty")
    data_table.columns = data_table.columns.str.strip()

    needed_col = {"Spectrum", "Spectrum File", "Peptide", "Charge", "PeptideProphet Probability", "Protein"}

    no_col = needed_col - set(data_table.columns)

    if no_col:
        raise ValueError("no col: {no_col}")

    return data_table



def extract_linear_scan_id(spectrum_id: str) -> int:
    if not isinstance(spectrum_id, str):
        raise TypeError("spectrum_id is not string")

    spectrum_id = spectrum_id.strip()

    if not spectrum_id:
        raise ValueError("spectrum id is empty")

    parts = spectrum_id.rsplit(".", maxsplit=3)

    if len(parts) != 4:
        raise ValueError("spectra_id is incorrect: {spectrum_id}")

    first_scan_text = parts[-3]
    last_scan_text = parts[-2]

    try:
        first_scan = int(first_scan_text)
        last_scan = int(last_scan_text)

    except ValueError:
        raise ValueError(f"сannot extract scan number from: {spectrum_id}")


    if first_scan <= 0 or last_scan <= 0:
        raise ValueError(f"scan num is <0")

    if first_scan != last_scan:
        raise ValueError(" first and last scan num must match")

    return first_scan



def prep_linear_results(data_table: pd.DataFrame) -> pd.DataFrame:

    if not isinstance(data_table, pd.DataFrame):
        raise TypeError("data_table must be a pandas DataFrame")

    column_map = {
        "Spectrum": "spectrum_id",
        "Spectrum File": "spectrum_file",
        "Peptide": "peptide",
        "Modified Peptide": "modified_peptide",
        "Charge": "charge_reported",
        "Retention": "retention_time_seconds_reported",
        "Observed M/Z": "precursor_mz_reported",
        "PeptideProphet Probability": "psm_probability",
        "Expectation": "expectation",
        "Hyperscore": "hyperscore",
        "Protein": "protein",
        "Assigned Modifications": "assigned_modifications"
    }

    changed_data_table = data_table.rename(columns=column_map).copy()
    
    text_columns = ["spectrum_id", "spectrum_file", "peptide", "modified_peptide", "protein", "assigned_modifications"]
    num_columns = ["charge_reported", "retention_time_seconds_reported", "precursor_mz_reported", "psm_probability", "expectation", "hyperscore"]

    for column in text_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = changed_data_table[column].astype("string").str.strip()
            

    for column in num_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = pd.to_numeric(changed_data_table[column], errors="coerce")

    changed_data_table["scan_id"] = changed_data_table["spectrum_id"].apply(extract_linear_scan_id)
    spectrum_run_names = changed_data_table["spectrum_id"].str.rsplit(".", n=3).str[0]
    changed_data_table["run_id"] = spectrum_run_names.apply(lambda run_name: Path(run_name).stem)
    

    return changed_data_table



def filter_confident_linear_psms(data_table: pd.DataFrame, min_probability: float = 0.99,) -> pd.DataFrame:
    no_columns = {"run_id", "scan_id", "peptide", "psm_probability", "protein"} - set(data_table.columns)

    if no_columns:
        raise ValueError(f"no columns: {sorted(no_columns)}")

    if not 0 <= min_probability <= 1:
        raise ValueError("min_probability must be between [0, 1]")
    

    filtered_data = data_table.copy()

    peptide_values = filtered_data["peptide"].astype("string").str.strip()
    protein_values = filtered_data["protein"].astype("string").str.strip()
    
    probability_values = pd.to_numeric(filtered_data["psm_probability"], errors="coerce")

    filtered_data["peptide"] = peptide_values
    filtered_data["protein"] = protein_values
    filtered_data["psm_probability"] = probability_values

    peptide_present = peptide_values.notna() & peptide_values.str.len().gt(0)
    protein_present = protein_values.notna() & protein_values.str.len().gt(0)

    run_id_present = filtered_data["run_id"].notna() & filtered_data["run_id"].astype("string").str.len().gt(0)
    scan_id_valid = filtered_data["scan_id"].notna() & filtered_data["scan_id"].gt(0)

    probability_valid = probability_values.notna() & probability_values.ge(min_probability) & probability_values.le(1)
    is_decoy = protein_values.str.lower().str.contains(r"^(rev_|decoy_|reverse_)",regex=True,na=False)
    mask = peptide_present & protein_present & run_id_present & scan_id_valid & probability_valid & ~is_decoy
    

    result = filtered_data.loc[mask].copy().reset_index(drop=True)
    
    if result.empty:
        raise ValueError("no linear PSMs")

    return result



def select_best_psm_per_scan(data_table: pd.DataFrame) -> pd.DataFrame:
    no_columns = {"run_id", "scan_id", "psm_probability"} - set(data_table.columns)
    
    if no_columns:
        raise ValueError(f"no columns: {sorted(no_columns)}")

    sort_columns = ["run_id", "scan_id", "psm_probability"]
    ascending = [True, True, False]
    
    if "hyperscore" in data_table.columns:
        sort_columns.append("hyperscore")
        ascending.append(False)

    if "expectation" in data_table.columns:
        sort_columns.append("expectation")
        ascending.append(True)


    sorted_data = data_table.sort_values(by=sort_columns, ascending=ascending, na_position="last")

    best_psms = sorted_data.drop_duplicates(subset=["run_id", "scan_id"], keep="first")

    return best_psms.reset_index(drop=True)



def exclude_crosslinked_scans(linear_results: pd.DataFrame, positive_labels: pd.DataFrame) -> pd.DataFrame:
    missing_linear_columns = {"run_id",
        "scan_id"} - set(linear_results.columns)
    
    missing_positive_columns = {"run_id",
        "scan_id"} - set(positive_labels.columns)
    
    if missing_linear_columns and missing_positive_columns:
        raise ValueError("not all columns are present")



    positive_data = positive_labels.copy()

    if "label" in positive_data.columns:
        positive_data = positive_data.loc[positive_data["label"].eq(1)]

    positive_keys = positive_data[["run_id", "scan_id"]].drop_duplicates().copy()
    positive_keys["_is_crosslinked"] = True

    merged_data = linear_results.merge(positive_keys,on=["run_id", "scan_id"], how="left")

    safe_linear_results = merged_data.loc[merged_data["_is_crosslinked"].isna()].drop(columns="_is_crosslinked").reset_index(drop=True)

    if safe_linear_results.empty:
        raise ValueError("no linear PSMs remained")

    return safe_linear_results



def match_linear_spectra(linear_results: pd.DataFrame, spectra: list[dict]) -> list[dict]:
    no_columns = {"run_id", "scan_id", "peptide"} - set(linear_results.columns)

    if no_columns:
        raise ValueError(f"no columns: {sorted(no_columns)}")

    if not isinstance(spectra, list):
        raise TypeError("spectra must be a list")

    if not spectra:
        raise ValueError("spectra list is empty")


    spectrum_lookup = {}
    for spectrum in spectra:
        if not isinstance(spectrum, dict):
            raise TypeError("spectrum must be a dict")

        no_keys = {"run_id", "scan_id", "mz_array", "intensity_array"} - set(spectrum)
        if no_keys:
            raise KeyError( f"no keys: {sorted(no_keys)}")

        key = str(spectrum["run_id"]).strip(), int(spectrum["scan_id"])
        if key in spectrum_lookup:
            raise ValueError(f"duplicate spectrum key: {key}")

        spectrum_lookup[key] = spectrum

    metadata_columns = ["peptide", "modified_peptide", "protein","psm_probability", "hyperscore", "expectation", "charge_reported", "precursor_mz_reported", "retention_time_seconds_reported", "spectrum_file", "assigned_modifications"]
    matched_spectra = []
    no_spectrum_keys = []

    for i, row in linear_results.iterrows():
        key = str(row["run_id"]).strip(), int(row["scan_id"]),
    
        if key not in spectrum_lookup:
            no_spectrum_keys.append(key)
            continue

        matched_spectrum = spectrum_lookup[key].copy()

        for column in metadata_columns:
            if column in linear_results.columns:
                matched_spectrum[column] = row[column]

        matched_spectrum["source_type"] = "linear"
        matched_spectra.append(matched_spectrum)

    if no_spectrum_keys:
        raise ValueError(f"no spectrum keys")

    if not matched_spectra:
        raise ValueError("no linear PSMs")

    return matched_spectra



def load_linear_spectra(results_path: str | Path, spectra: list[dict], positive_labels: pd.DataFrame, min_probability: float = 0.99) -> list[dict]:
    raw_results = read_linear_pep_data(results_path)
    prepared_results = prep_linear_results(raw_results)

    confident_results = filter_confident_linear_psms(prepared_results, min_probability=min_probability)
    best_results = select_best_psm_per_scan(confident_results)
    safe_results = exclude_crosslinked_scans(best_results, positive_labels)
    linear_spectra = match_linear_spectra(safe_results, spectra)
    

    return linear_spectra