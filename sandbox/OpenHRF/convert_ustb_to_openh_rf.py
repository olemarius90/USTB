"""Convert USTB UFF datasets to OpenH-RF (zea) HDF5 format.

Reads .uff files from Zenodo (or local cache) via pyuff-ustb and writes
.hdf5 files using the zea library in the OpenH-RF data format.

Usage:
    python convert_ustb_to_openh_rf.py --input data/ --output openh_rf/
    python convert_ustb_to_openh_rf.py --zenodo --output openh_rf/

Requirements:
    pip install pyuff-ustb zea numpy
"""

import argparse
import os
from pathlib import Path

import numpy as np

ZENODO_RECORD = "20261898"
ZENODO_BASE = f"https://zenodo.org/records/{ZENODO_RECORD}/files"

DATASETS = [
    # In-vivo human (cardiac, phased array P2-4)
    {"filename": "Verasonics_P2-4_parasternal_long_subject_1.uff", "probe": "P2-4", "tier": "in-vivo-human", "anatomy": "heart", "view": "parasternal long"},
    {"filename": "Verasonics_P2-4_parasternal_long_small.uff", "probe": "P2-4", "tier": "in-vivo-human", "anatomy": "heart", "view": "parasternal long"},
    {"filename": "Verasonics_P2-4_apical_four_chamber_subject_1.uff", "probe": "P2-4", "tier": "in-vivo-human", "anatomy": "heart", "view": "apical four-chamber"},
    # In-vivo human (carotid, linear array L7-4)
    {"filename": "L7_FI_carotid_cross_1.uff", "probe": "L7-4", "tier": "in-vivo-human", "anatomy": "carotid", "view": "cross-section"},
    {"filename": "L7_FI_carotid_cross_2.uff", "probe": "L7-4", "tier": "in-vivo-human", "anatomy": "carotid", "view": "cross-section"},
    {"filename": "L7_FI_carotid_cross_sub_2.uff", "probe": "L7-4", "tier": "in-vivo-human", "anatomy": "carotid", "view": "cross-section"},
    {"filename": "PICMUS_carotid_long.uff", "probe": "L7-4", "tier": "in-vivo-human", "anatomy": "carotid", "view": "longitudinal"},
    {"filename": "PICMUS_carotid_cross.uff", "probe": "L7-4", "tier": "in-vivo-human", "anatomy": "carotid", "view": "cross-section"},
    # Phantom (CIRS, tissue-mimicking)
    {"filename": "PICMUS_experiment_resolution_distortion.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "CIRS phantom", "view": "resolution"},
    {"filename": "PICMUS_experiment_contrast_speckle.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "CIRS phantom", "view": "contrast"},
    {"filename": "experimental_STAI_dynamic_range.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "dynamic range"},
    {"filename": "experimental_dynamic_range_phantom.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "dynamic range"},
    {"filename": "STAI_UFF_CIRS_phantom.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "CIRS phantom", "view": "STA"},
    {"filename": "L7_CPWC_193328.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "CPWC"},
    {"filename": "L7_FI_IUS2018.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "focused imaging"},
    {"filename": "L7_FI_Verasonics.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "focused imaging"},
    {"filename": "L7_FI_Verasonics_CIRS.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "CIRS phantom", "view": "focused imaging"},
    {"filename": "L7_FI_Verasonics_CIRS_points.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "CIRS phantom", "view": "point targets"},
    {"filename": "Alpinion_L3-8_FI_hypoechoic.uff", "probe": "L3-8", "tier": "phantom", "anatomy": "phantom", "view": "hypoechoic cyst"},
    {"filename": "Alpinion_L3-8_FI_hyperechoic_scatterers.uff", "probe": "L3-8", "tier": "phantom", "anatomy": "phantom", "view": "hyperechoic scatterers"},
    {"filename": "Alpinion_L3-8_CPWC_hypoechoic.uff", "probe": "L3-8", "tier": "phantom", "anatomy": "phantom", "view": "CPWC hypoechoic"},
    {"filename": "Alpinion_L3-8_CPWC_hyperechoic_scatterers.uff", "probe": "L3-8", "tier": "phantom", "anatomy": "phantom", "view": "CPWC hyperechoic"},
    {"filename": "ARFI_dataset.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "ARFI push"},
    {"filename": "SWE_L7_type_I.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "elastography phantom", "view": "shear wave"},
    {"filename": "SWE_L7_type_III.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "elastography phantom", "view": "shear wave"},
    {"filename": "SWE_L7_type_IV.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "elastography phantom", "view": "shear wave"},
    {"filename": "reference_RTB_data.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "RTB reference"},
    {"filename": "P4_FI_121444_45mm_focus.uff", "probe": "P4-1", "tier": "phantom", "anatomy": "phantom", "view": "focused phased array"},
    # Simulation (Field II)
    {"filename": "PICMUS_simulation_resolution_distortion.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "resolution"},
    {"filename": "PICMUS_simulation_contrast_speckle.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "contrast"},
    {"filename": "PICMUS_numerical_calib_v2.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "calibration"},
    {"filename": "FieldII_CPWC_simulation_v2.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "CPWC"},
    {"filename": "FieldII_CPWC_point_scatterers_res_v2.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "point scatterers"},
    {"filename": "FieldII_STAI_dynamic_range.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "dynamic range"},
    {"filename": "FieldII_STAI_simulated_dynamic_range.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "dynamic range"},
    {"filename": "FieldII_STAI_uniform_fov.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "uniform FOV"},
    {"filename": "FieldII_P4_point_scatterers.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "point scatterers"},
    {"filename": "FI_P4_point_scatterers.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "point scatterers"},
    {"filename": "FI_P4_cysts_center.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "cyst"},
    {"filename": "FieldII_speckle_simulation.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "speckle"},
    {"filename": "FieldII_speckle_DMASsimulation300000pts.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "DMAS speckle"},
    {"filename": "speckle_sim_FI_P4_probe_apod_1_speckle_long_many_angles.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "blocked array (full)"},
    {"filename": "speckle_sim_FI_P4_probe_apod_2_speckle_long_many_angles.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "blocked array (1/3)"},
    {"filename": "speckle_sim_FI_P4_probe_apod_3_speckle_long_many_angles.uff", "probe": "P4-1", "tier": "simulation", "anatomy": "simulated", "view": "blocked array (1/2)"},
    # Generalized Beamformer datasets
    {"filename": "L7_CPWC_TheGB.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "CPWC"},
    {"filename": "L7_FI_TheGB.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "focused imaging"},
    {"filename": "L7_DW_TheGB.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "diverging wave"},
    {"filename": "L7_STA_TheGB.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "synthetic transmit aperture"},
    # Large in-vivo
    {"filename": "invitro_20.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "tissue phantom", "view": "in-vitro"},
    {"filename": "insilico_20.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "in-silico"},
    {"filename": "insilico_side_100_M45.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "in-silico side"},
    {"filename": "beamformed_simulated_data.uff", "probe": "L7-4", "tier": "simulation", "anatomy": "simulated", "view": "DRA beamformed"},
    {"filename": "beamformed_experimental_data.uff", "probe": "L7-4", "tier": "phantom", "anatomy": "phantom", "view": "DRA beamformed"},
]


def download_if_missing(filename: str, data_dir: Path) -> Path:
    """Download .uff from Zenodo if not cached locally."""
    import urllib.request

    filepath = data_dir / filename
    if filepath.exists():
        return filepath
    data_dir.mkdir(parents=True, exist_ok=True)
    url = f"{ZENODO_BASE}/{filename}"
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, filepath)
    return filepath


def read_uff_channel_data(filepath: Path):
    """Read channel_data from a UFF file using pyuff-ustb."""
    from pyuff_ustb.objects.uff import Uff

    uff_file = Uff(str(filepath))
    try:
        return uff_file.read("channel_data")
    except Exception:
        return None


def uff_to_zea_scan(channel_data) -> dict:
    """Map pyuff-ustb ChannelData to a zea-compatible scan dictionary."""
    probe = channel_data.probe
    sequence = channel_data.sequence

    n_el = int(probe.N_elements)
    probe_geometry = np.zeros((n_el, 3), dtype=np.float32)
    probe_geometry[:, 0] = np.array(probe.x, dtype=np.float32).flatten()
    if hasattr(probe, 'y') and probe.y is not None:
        probe_geometry[:, 1] = np.array(probe.y, dtype=np.float32).flatten()
    if hasattr(probe, 'z') and probe.z is not None:
        probe_geometry[:, 2] = np.array(probe.z, dtype=np.float32).flatten()

    n_tx = len(sequence)
    sampling_frequency = float(channel_data.sampling_frequency)
    sound_speed = float(channel_data.sound_speed)

    polar_angles = np.zeros(n_tx, dtype=np.float32)
    azimuth_angles = np.zeros(n_tx, dtype=np.float32)
    focus_distances = np.full(n_tx, np.inf, dtype=np.float32)
    initial_times = np.zeros(n_tx, dtype=np.float32)

    for i, wave in enumerate(sequence):
        src = wave.source
        if hasattr(src, 'azimuth') and src.azimuth is not None:
            polar_angles[i] = float(src.azimuth)
        if hasattr(src, 'elevation') and src.elevation is not None:
            azimuth_angles[i] = float(src.elevation)
        if hasattr(src, 'distance') and src.distance is not None:
            focus_distances[i] = float(src.distance)
        if hasattr(wave, 'delay') and wave.delay is not None:
            initial_times[i] = float(wave.delay)

    center_frequency = 0.0
    if hasattr(channel_data, 'pulse') and channel_data.pulse is not None:
        if hasattr(channel_data.pulse, 'center_frequency'):
            center_frequency = float(channel_data.pulse.center_frequency or 0)

    scan = {
        "probe_geometry": probe_geometry,
        "sampling_frequency": sampling_frequency,
        "center_frequency": center_frequency,
        "initial_times": initial_times,
        "t0_delays": np.zeros((n_tx, n_el), dtype=np.float32),
        "tx_apodizations": np.ones((n_tx, n_el), dtype=np.float32),
        "focus_distances": focus_distances,
        "transmit_origins": np.zeros((n_tx, 3), dtype=np.float32),
        "polar_angles": polar_angles,
        "azimuth_angles": azimuth_angles,
        "sound_speed": sound_speed,
    }
    return scan


def convert_dataset(filepath: Path, output_dir: Path, meta: dict) -> bool:
    """Convert a single UFF file to OpenH-RF zea HDF5."""
    from zea import File

    channel_data = read_uff_channel_data(filepath)
    if channel_data is None:
        print(f"  SKIP (no channel_data): {filepath.name}")
        return False

    raw = np.array(channel_data.data, dtype=np.float32)
    # UFF shape: [samples, channels, waves, frames] -> [frames, tx, samples, elements, 1]
    if raw.ndim == 4:
        raw = raw.transpose(3, 2, 0, 1)  # [frames, waves, samples, channels]
    elif raw.ndim == 3:
        raw = raw.transpose(2, 0, 1)  # [waves, samples, channels]
        raw = raw[np.newaxis, ...]  # [1, waves, samples, channels]
    else:
        print(f"  SKIP (unexpected shape {raw.shape}): {filepath.name}")
        return False

    raw = raw[..., np.newaxis]  # add trailing dim -> [frames, tx, samples, el, 1]

    scan = uff_to_zea_scan(channel_data)
    output_path = output_dir / filepath.with_suffix(".hdf5").name

    metadata = {
        "annotations": {
            "anatomy": meta.get("anatomy", ""),
            "view": meta.get("view", ""),
        },
        "credit": "USTB - UltraSound ToolBox (University of Oslo)",
        "subject": {"type": meta.get("tier", "unknown")},
    }

    File.create(
        path=str(output_path),
        data={"raw_data": raw},
        scan=scan,
        metadata=metadata,
        probe_name=meta.get("probe", "unknown"),
        us_machine="Verasonics / Alpinion / Field II simulation",
        description=f"USTB dataset: {filepath.stem}. {meta.get('view', '')}. "
                    f"Tier: {meta.get('tier', '')}.",
    )
    print(f"  OK: {output_path.name} ({raw.shape})")
    return True


def main():
    parser = argparse.ArgumentParser(description="Convert USTB UFF to OpenH-RF")
    parser.add_argument("--input", type=str, default="data",
                        help="Directory with local .uff files")
    parser.add_argument("--output", type=str, default="openh_rf",
                        help="Output directory for .hdf5 files")
    parser.add_argument("--zenodo", action="store_true",
                        help="Download from Zenodo if file not in --input")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit number of datasets to convert (0=all)")
    args = parser.parse_args()

    data_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = DATASETS[:args.limit] if args.limit > 0 else DATASETS
    ok_count = 0

    for meta in datasets:
        fn = meta["filename"]
        filepath = data_dir / fn
        if not filepath.exists() and args.zenodo:
            filepath = download_if_missing(fn, data_dir)
        if not filepath.exists():
            print(f"  MISSING: {fn}")
            continue
        try:
            if convert_dataset(filepath, output_dir, meta):
                ok_count += 1
        except Exception as e:
            print(f"  ERROR ({fn}): {e}")

    print(f"\nConverted {ok_count}/{len(datasets)} datasets to {output_dir}/")


if __name__ == "__main__":
    main()
