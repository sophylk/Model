from pyteomics import mzml
from pathlib import Path
import re

def read_spectra(path: Path) -> list[dict]:

    if path.suffix.lower() != ".mzml":
        raise ValueError("mzML file needed")
    if not path.exists():
        raise FileNotFoundError("no file")

  
    scans_info = []

    with mzml.read(str(path)) as data:
        for scan in data:
            if scan['ms level'] == 2:

                selected_ion = (
                    scan["precursorList"]
                    ["precursor"][0]
                    ["selectedIonList"]
                    ["selectedIon"][0]
                )


                scan_id = re.search(r"scan=(\d+)", scan["id"])
                if scan_id is None:
                    raise ValueError(f"No scan id: {scan['id']}")

                scan_id = int(scan_id.group(1))

                charge = selected_ion.get("charge state")

                if charge is not None:
                    charge = int(charge)

                
                scans_info.append({
                    "scan_id": scan_id,
                    "spectrum_id": scan['id'],
                    "ms_level": scan['ms level'],
                    "precursor_mz": float(selected_ion["selected ion m/z"]),
                    "charge": charge,
                    "retention_time": float(scan["scanList"]["scan"][0]["scan start time"]),
                    "mz_array": scan['m/z array'],
                    "intensity_array": scan['intensity array'],
                    "source_file": path.name,
                    "run_id": path.stem
                })



    return scans_info


