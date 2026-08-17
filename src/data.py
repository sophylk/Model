from pyteomics import mzml
from pathlib import Path
import re

def read_spectra(path: Path) -> list[dict]:

    if path.suffix.lower() != ".mzml":
        raise ValueError("mzML file needed")
    elif not path.exists():
        raise FileNotFoundError("no file")

    
    scans_info = []

    with mzml.read(path) as data:
        for scan in data:
            if scan['ms level'] == 2:

                selected_ion = (
                    scan["precursorList"]
                    ["precursor"][0]
                    ["selectedIonList"]
                    ["selectedIon"][0]
                )

                
                scans_info.append({
                    "scan_id": int(re.search(r'scan=(\d+)', scan['id']).group(1)),
                    "spectrum_id": scan['id'],
                    "ms_level": scan['ms level'],
                    "precursor_mz": float(selected_ion["selected ion m/z"]),
                    "charge": selected_ion.get("charge state"),
                    "retention_time": float(scan["scanList"]["scan"][0]["scan start time"]),

                    "mz_array": scan['m/z array'].tolist(),

                    "intensity_array": scan['intensity array'].tolist(),
                    "source_file": path.name,
                })



    return scans_info



def extract_scan_id(spectrum_id):
    
    pass


def parse_spectrum(spectrum, source_file):
    
    pass

