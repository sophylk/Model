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

    needed_col = { "Spectrum",   
            "Spectrum File",
            "Peptide",
            "Charge",
            "PeptideProphet Probability",
            "Protein"
        }

    no_col = needed_col - set(data_table.columns)

    if no_col:
        raise ValueError("no col: {no_col}")


    return data_table


def extract_linear_scan_id(spectrum_id: str) -> int:

    if not isinstance(spectrum_id, str):
        raise TypeError("spectrum_id is not string")

    spectrum_id = spectrum_id.strip()

    if not spectrum_id:
        raise ValueError("spectrum_id is empty")

    parts = spectrum_id.rsplit(".", maxsplit=3)

    if len(parts) != 4:
        raise ValueError("spectra_id is incorrect: {spectrum_id}")

    first_scan_text = parts[-3]
    last_scan_text = parts[-2]

    if (type(first_scan_text) == int and type(last_scan_text) == int):
        first_scan = int(first_scan_text)
        last_scan = int(last_scan_text)
    else:
        raise ValueError(f"cannot extract scan number from: {spectrum_id}")


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
    numeric_columns = ["charge_reported", "retention_time_seconds_reported", "precursor_mz_reported", "psm_probability", "expectation", "hyperscore"]

    for column in text_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = changed_data_table[column].astype("string").str.strip()
            

    for column in numeric_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = pd.to_numeric(changed_data_table[column], errors="coerce")

    changed_data_table["scan_id"] = changed_data_table["spectrum_id"].apply(extract_linear_scan_id)
    spectrum_run_names = changed_data_table["spectrum_id"].str.rsplit(".", n=3).str[0]
    changed_data_table["run_id"] = spectrum_run_names.apply(lambda run_name: Path(run_name).stem)
    

    return changed_data_table



def filter_confident_linear_psms(data_table: pd.DataFrame) -> pd.DataFrame:
    pass


def select_best_psm_per_scan(data_table: pd.DataFrame) -> pd.DataFrame:
    pass


def exclude_crosslinked_scans(linear_results: pd.DataFrame, positive_labels: pd.DataFrame) -> pd.DataFrame:
    pass


def match_linear_spectra(linear_results: pd.DataFrame, spectra: list[dict]) -> list[dict]:
    pass


def load_linear_spectra(results_path: str | Path, spectra: list[dict], positive_labels: pd.DataFrame) -> list[dict]:
    pass