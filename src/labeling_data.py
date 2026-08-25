from pyteomics import mzml
from pathlib import Path
import pandas as pd



def read_xlink_data(path: Path) -> pd.DataFrame:
    path = Path(path)

    if path.suffix.lower() != ".csv":
            raise ValueError("mzML file needed")
    if not path.exists():
        raise FileNotFoundError("no file")
    
    data_table = pd.read_csv(path)

    if data_table.empty:
        raise ValueError(f"table is empty")

    data_table.columns = data_table.columns.str.strip()

    needed_col = {
        "First Scan",
        "XlinkX Score",
        "Delta XlinkX Score",
        "Sequence A",
        "Sequence B",
        "Is Decoy",
        "Crosslinker",
        "Spectrum file path",
    }

    
    no_columns = needed_col - set(data_table.columns)
    
    if no_columns:
        raise ValueError("No columns: "f"{sorted(no_columns)}")

    return data_table


def adapt_run_id(file_path: str) -> str:

    if pd.isna(file_path):
        raise ValueError("empty file path")

    file_name = str(file_path).replace("\\", "/").rsplit("/", maxsplit=1)[-1]

    return Path(file_name).stem


def prepare_data(data_table: pd.DataFrame) -> pd.DataFrame:

    column_map = {
        "First Scan": "scan_id",
        "XlinkX Score": "xlinkx_score",
        "Delta XlinkX Score": "delta_score",
        "Sequence A": "peptide_a",
        "Sequence B": "peptide_b",
        "Is Decoy": "is_decoy",
        "Crosslinker": "crosslinker",
        "Crosslink Type": "crosslink_type",
        "Spectrum file path": "spectrum_file",
        "# Identified MS2 Scans": "identified_ms2_count",
        "Charge": "charge_reported",
        "m/z [Da]": "precursor_mz_reported",
        "RT [min]": "retention_time_reported",
    }

    changed_data_table = data_table.rename(columns=column_map)

    numeric_columns = ["scan_id", "xlinkx_score", "delta_score", "identified_ms2_count", "charge_reported", "precursor_mz_reported", "retention_time_reported"]
    text_columns = ["peptide_a", "peptide_b", "crosslinker", "crosslink_type", "spectrum_file"]


    for column in numeric_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = pd.to_numeric(changed_data_table[column], errors="coerce")


    for column in text_columns:
        if column in changed_data_table.columns:
            changed_data_table[column] = changed_data_table[column].astype("string").str.strip()

    changed_data_table["run_id"] = changed_data_table["spectrum_file"].apply(adapt_run_id)
    

    return changed_data_table


def filter_crosslinks(data_table: pd.DataFrame, min_score: float = 20.0, min_delta_score: float = 20.0) -> pd.DataFrame:
    decoy_values = data_table["is_decoy"].astype("string").str.strip().str.lower()
    is_not_decoy = decoy_values.isin({ "false","0","no","n"})
    is_dsso = data_table["crosslinker"].astype("string").str.contains("DSSO", case=False, na=False)

    peptide_a_present = data_table["peptide_a"].notna() & data_table["peptide_a"].str.len() > 0
    peptide_b_present = data_table["peptide_b"].notna() & data_table["peptide_b"].str.len() > 0
    

    check = data_table["scan_id"].notna() & data_table["xlinkx_score"].ge(min_score) & data_table["delta_score"].ge(min_delta_score) & is_not_decoy & is_dsso & peptide_a_present & peptide_b_present

    return data_table.loc[check].copy()


def crosslinked_labels(data_table: pd.DataFrame) -> pd.DataFrame:
    labels_table = data_table.copy()

    labels_table["scan_id"] = labels_table["scan_id"].astype(int)
    labels_table["label"] = 1

    labels_table = labels_table.sort_values(by=["xlinkx_score", "delta_score"], ascending=[False, False])
    labels_table = labels_table.drop_duplicates(subset=["run_id", "scan_id"], keep="first")

    columns = ["run_id", "scan_id", "label", "peptide_a", "peptide_b", "xlinkx_score", "delta_score", "crosslinker", "crosslink_type"]
    

    extra_columns = ["identified_ms2_count", "charge_reported", "precursor_mz_reported", "retention_time_reported"]

    for column in extra_columns:
        if column in labels_table.columns:
            columns.extend(column)

    return labels_table[columns].reset_index(drop=True)


def load_positive_labels(path: str | Path, min_score: float = 20.0, min_delta_score: float = 20.0) -> pd.DataFrame:
    raw_data = read_xlink_data(path)
    prep_data = prepare_data(raw_data)
    filtered_data = filter_crosslinks(prep_data, min_score=min_score, min_delta_score=min_delta_score)
    
    return crosslinked_labels(filtered_data)