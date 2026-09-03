from pathlib import Path
import pandas as pd

def read_chimera_results(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() != ".tsv":
            raise ValueError("tsv file needed")
    if not path.exists():
        raise FileNotFoundError("no file")
    
    data_table = pd.read_csv(path)
    if data_table.empty:
        raise ValueError(f"table is empty")

    data_table.columns = data_table.columns.str.strip()

    needed_col = {"Spectrum", "Peptide", "Probability", "Protein"}
    no_col = needed_col - set(data_table.columns)

    if no_col:
         raise ValueError("no columns: {no_col}")
    

    return data_table


def extract_chimera_scan_id(spectrum_id: str) -> int:
     
    pass