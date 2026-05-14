# Publishing USTB website assets

Pre-built MATLAB content is uploaded **separately** as three release assets (`examples-v1` on GitHub). That keeps runs **isolatable**: fix or refresh **examples**, **publications**, or **datasets** independently.

## Three publishers (`examples-v1` assets)

| Script | Artifact | Contents / deploy target |
|---|---|---|
| `./publish_examples.sh` | **`examples-html.tar.gz`** | `publish_all_examples.m` gallery → **`website/examples/`** |
| `./publish_publications.sh` | **`publications-html.tar.gz`** | TUSON (etc.) MATLAB `publish` output → merged into **`website/examples/`** (e.g. `TUSON/`) |
| `./publish_datasets.sh` | **`datasets-html.tar.gz`** | `export_dataset_previews_to_website` PNGs + `build_datasets_page.py` → **`website/`** root (`datasets.html`, `assets/images/datasets/`)

Shared MATLAB flags / Git Bash paths: **`publish_common.sh`** (sourced automatically).

Upload (each script):

```bash
./publish_examples.sh --upload unioslo/USTB
./publish_publications.sh --upload unioslo/USTB
./publish_datasets.sh --upload unioslo/USTB
```

## Examples only — Quick Start

```bash
./publish_examples.sh                 # Examples gallery only → examples-html.tar.gz
./publish_publications.sh            # Publication HTML only → publications-html.tar.gz
./publish_datasets.sh                # Dataset PNGs + datasets.html → datasets-html.tar.gz
```

## Prerequisites

- **MATLAB R2024b** or later with Signal Processing Toolbox
- **Field II v3.30** (optional, for simulation examples) — install to `/opt/field_ii`
- **Python 3** (for `generate_examples_index.py`)
- **Recompiled MEX file** — the MEX beamformer must be compiled against your system's TBB:

```matlab
cd('+mex');
mex('-R2018a', '-D_UNIX_', '-I/usr/include/tbb', ...
    'LDFLAGS="$LDFLAGS -Wl,-rpath,/usr/lib/x86_64-linux-gnu"', ...
    '-L/usr/lib/x86_64-linux-gnu', '-ltbb', 'source/das_c.cpp');
```

## How It Works

1. `./publish_examples.sh` runs MATLAB `publish_all_examples.m` (**examples only** — dataset previews are `./publish_datasets.sh`)
2. Outputs with **`Error using` / `Error in `** in the HTML body are stripped by `publish_examples.sh`
3. `generate_examples_index.py` writes `index.html`
4. Pack → **`examples-html.tar.gz`**, upload alongside the other two assets on **`examples-v1`**
5. **`deploy-website.yml`** downloads **three** `.tar.gz` files (primary repo, then **`olemarius90/USTB`** fallback): examples → `website/examples/`; publications → `website/examples/` (overlay); datasets → **`website/`** overlay.


## Publications page (`publications.html`)

Some iframes load HTML that is **not** under `examples/` in the repository (source lives in `sandbox/`). Those are published by **`publish_all_examples.m`** into the same `examples_html/` tree so the deployed path matches the site, e.g.:

| Website path | Source in repo |
|---|---|
| `examples/generalized_beamformer/CPWC_double_adaptive_redone.html` | `sandbox/The_Generalized_Beamformer/CPWC_double_adaptive_redone.m` |

After changing `publish_all_examples.m`, run `./publish_examples.sh` and upload **`examples-html.tar.gz`** to the **`examples-v1`** release on **`unioslo/USTB`** (and rely on CI `curl` fallback to your fork if needed).

```bash
./publish_examples.sh
./publish_examples.sh --upload unioslo/USTB
```

Then trigger **Deploy Website** on `unioslo/USTB` `master` (or merge the workflow paths fix) so GitHub Pages picks up the new tarball.

### Publication pages (`publications/`, TUSON, etc.)

Scripts under `publications/` are **not** part of `publish_all_examples` (they live outside `examples/`). Build and upload them with:

```bash
./publish_publications.sh
./publish_publications.sh --upload
# Or for the upstream repo:
./publish_publications.sh --upload unioslo/USTB
```

This produces `publications-html.tar.gz` (HTML + `publish()` figure PNGs) on **`examples-v1`**. **`deploy-website.yml`** overlays it onto **`website/examples/`** **after** the examples tarball, so paths like `website/examples/TUSON/.../Correction_of_simulated_blockage.html` match the **`publications.html`** iframes.

### Dataset page (`datasets.html` + PNGs)

Run **`./publish_datasets.sh`** (then **`--upload`**). That produces **`datasets-html.tar.gz`**, unpacked with **`-C website`** during deploy (**`datasets.html`** + **`assets/images/datasets/`**).

## Windows (Git Bash)

- **Git Bash paths:** the repo root appears as **`/c/...`** but Windows MATLAB **`cd()`** expects **`C:/...`** (or `\`). The **`publish_*.sh`** scripts source **`publish_common.sh`** which rewrites paths (`cygpath -m` or **`/c/x` → `C:/x`**) before **`-batch`**, otherwise **`publish_all_examples`** is “not found” while MATLAB says it exists under **`C:\...\ustb`**.
- MATLAB prints **"Unrecognized command line option: nodisplay"** if you use **Linux-only** `-nodisplay`. `publish_examples.sh` omits that flag when **Windows** is detected (`WINDIR` / `SYSTEMROOT`, or `OSTYPE` msys/cygwin/mingw) and only adds it on real Linux/macOS shells.
- **`slsc_mex.mexw64`** errors ("side-by-side configuration is incorrect") mean a **Visual C++ runtime** mismatch — reinstall the MSVC redist MATLAB ships with, or **rebuild** `+mex/slsc_mex` from source. Examples that depend on SLSC are **skipped** in `publish_all_examples.m` until the MEX loads.
- **`export_dataset_previews_to_website.m`** is **not** part of `./publish_examples.sh` — use **`./publish_datasets.sh`**.

## Manual Steps

```bash
# 1. Install Field II (optional)
mkdir -p /opt/field_ii
curl -L -A "Mozilla/5.0" -o /tmp/field_ii.tar.gz \
  "https://www.field-ii.dk/program_code/matlab_2021/Field_II_ver_3_30_linux.tar.gz"
tar -xzf /tmp/field_ii.tar.gz -C /opt/field_ii

# 2. Recompile MEX (if needed)
matlab -batch "cd('+mex'); mex('-R2018a','-D_UNIX_','-I/usr/include/tbb','LDFLAGS=\"\$LDFLAGS -Wl,-rpath,/usr/lib/x86_64-linux-gnu\"','-L/usr/lib/x86_64-linux-gnu','-ltbb','source/das_c.cpp')"

# 3. Publish examples (gallery only); upload if desired
./publish_examples.sh --upload unioslo/USTB

# 4. Publications (optional, separate tarball)
./publish_publications.sh --upload unioslo/USTB

# 5. Dataset preview PNGs + datasets.html (optional, separate tarball; long run)
./publish_datasets.sh --upload unioslo/USTB
```

## Skipped Examples

These examples are skipped by `publish_all_examples.m` (not attempted):

### External toolbox dependencies

| Area | Reason |
|---|---|
| `examples/FLUST/*` | Needs MUST toolbox |
| `examples/kWave/*` | Needs k-Wave toolbox |
| `examples/REFoCUS/*` | Causes segfault in headless MATLAB |
| `examples/field_II/*` | Field II **`field_init`** not available on typical dev machines |

### Hardware/data dependencies

| Example | Reason |
|---|---|
| `examples/verasonics/*` | Needs Verasonics hardware/data |
| `examples/alpinion/*` | Needs Alpinion hardware/data |
| `examples/acoustical_radiation_force_imaging/*` | Needs hardware data |

### Interactive or slow

| Example | Reason |
|---|---|
| `MATLAB_intro.m` | Uses `ginput()`, hangs in headless |
| `FI_L11_parfor_compared_to_fresnel.m` | Needs Parallel Computing Toolbox |
| `STAI_L11_speckle_parfor.m` | Needs Parallel Computing Toolbox |
| `FI_P4_cardiac_coherence.m` | Needs Parallel Computing Toolbox |
| `STAI_2D_array_cardiac.m` | 3D simulation, very long runtime |
| `CPWC_2D_array_cardiac.m` | 3D simulation, very long runtime |

*(Field II slowdowns live under **`examples/field_II/`**, which is skipped entirely — see External toolbox dependencies.)*

### Publishing batch (helpers, heavy downloads, SLSC MEX)

| Example | Reason |
|---|---|
| `dataset_preview_beamform.m` | Requires input arguments; not standalone |
| `export_dataset_previews_to_website.m` | Run **`publish_datasets.sh`** only — heavy multi-download previews |
| `export_png_like_b_data_plot.m` | Requires input arguments; not standalone |
| `dataset_smoke_test_all.m` | Many downloads; very long in batch |
| `FI_UFF_generalized_coherence_factor.m` | `slsc_mex` (often fails on Windows if VC++ runtime / MEX broken) |
| `FI_UFF_short_lag_spatial_coherence.m` | `slsc_mex` |
| `FI_UFF_multi_frame_processing.m` | Dataset URL may return HTTP 303 until download helper is updated |
| `resolve_channel_data_path.m`, `simple_process_dataset.m`, `website_slug_for_dataset.m` | Dataset-smoke helpers; **not** standalone **`publish`** targets |
| `CPWC_UFF_read.m`, `CPWC_UFF_write.m`, `FI_UFF_phased_array.m`, `FI_UFF_Verasonics_MLA.m` | 404/`questdlg`/MLA apod brittle in batch (Windows / recent MATLAB) |

### Fresnel / GPU / MLA (skipped by basename)

STA/RTB/MLA/matrix-array/low-PW/GPU/multiframe demos that reliably throw when **`publish()`** runs **`evalCode`** in batch (recent MATLAB): see the **`skip_files`** list inside **`publish_all_examples.m`** (`CPWC_linear_array_beamformer_speed`, `STA_linear_array*.m`, **`FI_phased_array_MLA`**, **`CPWC_matrix_array`**, …).

### Course exercises

| Example | Reason |
|---|---|
| `examples/UiO_course_IN4015_Ultrasound_Imaging/*` | Student exercises, some with unimplemented code |

### Runtime errors (`publish_examples.sh` cleanup)

Examples that **`publish`** with **`catchError`** may still emit error text inside HTML until `publish_examples.sh` strips those pages. Prefer skipping brittle scripts (above) so you do not churn release tarballs full of placeholders.

Historical examples that used to surface here are enumerated in **`skip_files`** inside **`publish_all_examples.m`**; many are now **skipped** instead of attempted.

## Published gallery (indicative)

Output varies by MATLAB version/toolboxes — check **`examples_html/index.html`** after **`./publish_examples.sh`**.

| Category | Examples *(non-exhaustive)* |
|---|---|
| Advanced Beamforming | `FI_UFF_delay_multiply_and_sum_resolution` |
| Field II | **Not published by default** — whole **`examples/field_II/`** directory is **`skip_dirs`** unless you customize it locally |
| Fresnel / Curvilinear | `DW_curvilinear_array`, `FI_curvilinear_array` |
| Fresnel / Linear | e.g. `CPWC_linear_array_tilt`, `DW_linear_array`, `FI_linear_array`, `RTB_linear_array_close_up` *(STA/standard RTB/multi-frame/low-PW/GPU snippets are in **`skip_files`**)*
| Fresnel / Phased | e.g. `DW_phased_array`, `FI_phased_array`, `FI_phased_array_RTB*` *(some MLA / multiframe demos skipped)* |
| PICMUS | All six *(experiment/simulation × resolution/contrast/invivo)* when downloads succeed |
| UFF | e.g. `CPWC_UFF_Alpinion`, `FI_UFF_Alpinion`, `FI_UFF_Verasonics_RTB`, `STAI_UFF_beamform_with_demodulation` |
