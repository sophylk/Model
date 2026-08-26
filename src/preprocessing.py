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

    no_keys =  needed_keys - set(spectrum)
    if no_keys:
        raise KeyError("no "f"{no_keys}")

    mz_array = np.asarray(spectrum["mz_array"])
    intensity_array = np.asarray(spectrum["intensity_array"])

    if not (np.issubdtype(mz_array.dtype, np.number) and np.issubdtype(intensity_array.dtype, np.number)):
        raise ValueError("need np.number, not other type")
    if not (mz_array.ndim == 1 and intensity_array.ndim == 1):
        raise ValueError("need 1 dim, not >1")
    if mz_array.size == 0 or intensity_array.size == 0:
        raise ValueError("arrays are empty")

    if mz_array.size != intensity_array.size:
        raise ValueError("len of two arrays must match")
  
    return None


def clean_peaks(mz_array: np.ndarray, intensity_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mz = np.asarray(mz_array, dtype=np.float64)
    intensities = np.asarray(intensity_array, dtype=np.float32)

    valid_peaks = np.isfinite(mz) & np.isfinite(intensities) & (mz > 0) & (intensities > 0)
    clear_mz = mz[valid_peaks]
    clean_intensities = intensities[valid_peaks]

    if clear_mz.size == 0:
        raise ValueError("clear_mz array is empty")

    clear_sort = np.argsort(clear_mz)
    sorted_mz = clear_mz[clear_sort]
    sorted_intensities = clean_intensities[clear_sort]


    return sorted_mz, sorted_intensities


def normalize_intensities(intensity_array: np.ndarray) -> np.ndarray:

    intensities = np.asarray(intensity_array, dtype=np.float32)

    if intensities.size == 0:
        raise ValueError("intensity_array is empty")

    max_intensity = np.max(intensities)
    if max_intensity <= 0:
        raise ValueError("max_intensity cannot be 0")

    normalized_intensities_array = intensities / max_intensity

    return normalized_intensities_array


def preprocess_spectrum(spectrum: dict) -> dict:

    check_spectrum(spectrum)
    cleaned_mz, cleaned_intensities = clean_peaks(spectrum["mz_array"], spectrum["intensity_array"])
    fixed_intensities = normalize_intensities(cleaned_intensities)

    new_spectrum = spectrum.copy()
    new_spectrum["mz_array"], new_spectrum["intensity_array"] = cleaned_mz, fixed_intensities
    new_spectrum["peak_count"] = len(cleaned_mz)



    return new_spectrum