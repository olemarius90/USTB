# OpenH-RF Contribution: USTB Datasets

This folder contains the proposal and conversion tools for contributing
USTB (UltraSound ToolBox) datasets to the
[OpenH-RF](https://github.com/open-h/OpenH-RF) initiative.

## Contents

| File | Description |
|------|-------------|
| `proposal.tex` | 2-page LaTeX proposal for the OpenH-RF RFP |
| `convert_ustb_to_openh_rf.py` | Python script: UFF → zea HDF5 conversion |
| `generate_labels.m` | MATLAB script: DAS beamforming labels for all datasets |

## Quick start

### Build the proposal PDF

```bash
cd sandbox/OpenHRF
pdflatex proposal.tex
```

### Convert datasets

```bash
pip install pyuff-ustb zea numpy
python convert_ustb_to_openh_rf.py --zenodo --output openh_rf/
```

### Generate labels (MATLAB)

```matlab
addpath(genpath(ustb_path()));
run('sandbox/OpenHRF/generate_labels.m');
```

## Dataset source

All 53 datasets are hosted on Zenodo:
- **Record 20261898**: https://zenodo.org/records/20261898
- **License**: MIT (code), CC BY 4.0 (data contribution to OpenH-RF)

## Links

- OpenH-RF RFP: https://github.com/open-h/OpenH-RF
- USTB: https://github.com/unioslo/USTB
- USTB website: https://www.ultrasoundtoolbox.com
