import numpy as np


def check_spectrum(spectrum: dict) -> None:
    if not isinstance(spectrum, dict):
        raise TypeError("need only dict")
    if not spectrum:
        raise ValueError("dict empty")

    needed_keys = {
        "mz_array",
        "intensity_array"
    }

    no_keys = set(spectrum) - needed_keys
    if not no_keys:
        raise KeyError("no " f"{no_keys}")

    mz_array = np.asarray(spectrum["mz_array"])
    intensity_array = np.asarray(spectrum["intensity_array"])

    if not (np.issubdtype(mz_array.dtype, np.number) or np.issubdtype(intensity_array.dtype, np.number)):
        raise ValueError("need np.number, not other type")
    
    return


def clean_peaks(mz_array: np.ndarray, intensity_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    
    pass


def clean_intensities(intensity_array: np.ndarray) -> np.ndarray:
   
    pass


def preprocess_spectrum(spectrum: dict) -> dict:
    
    pass