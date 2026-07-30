# AGENTS.md

Guide for AI agents working in this repository.

## Project overview

HIV-Quant is a bioimage analysis pipeline that segments nuclei in 3D from multi-channel confocal z-stacks (`.vsi` files) and quantifies per-channel intensity within them. It was originally built for HIV capsid/CPSF6/HA intensity measurement but is general-purpose — all key parameters are CLI options.

**Single-file architecture**: Everything lives in `quantify_nuclei_intensity.py` plus a Jupyter notebook (`quantify_nuclei_intensity.ipynb`). No package structure, no tests, no CI.

## Environment management (pixi)

This project uses [pixi](https://pixi.sh), not pip/conda directly:

```bash
pixi install              # create/update the environment
pixi run python <script>  # run a script inside the environment
pixi shell                 # drop into an activated shell
```

Dependencies are split across:
- `[dependencies]` in `pixi.toml` — conda-forge packages (Python 3.14, pixi-pycharm)
- `[pypi-dependencies]` in `pixi.toml` — PyPI packages (numpy, pandas, scipy, bioio, matplotlib, scikit-image, seaborn, jupyter)

**Adding dependencies**: Edit `pixi.toml`, then run `pixi install` to regenerate `pixi.lock`. Commit both files together. Note that `pixi.lock` is gitignored.

## Running the pipeline

```bash
pixi run python quantify_nuclei_intensity.py --data-dir ./data
```

Run with `--help` for full option list. Key CLI options:
- `--data-dir`: Directory containing `.vsi` files (default: `./data`)
- `--channel-names`: Ordered list of channel names matching acquisition order; must include `DAPI` (default: `DAPI HA CPSF6 Capsid`)
- `--nuclei-diameter-px`: Expected nucleus diameter in pixels for size filtering (default: 140)
- `--size-tolerance`: Fractional deviation tolerance from expected diameter (default: 0.3)
- `--condition-mapping`: JSON object mapping numeric file index to condition label

## Pipeline stages (in `quantify_nuclei_intensity.py`)

The script runs linearly through these stages:

1. **CLI parsing** (`parse_args`): All configuration comes from CLI args or module-level constants
2. **File discovery**: Globs `*.vsi` files under `--data-dir`
3. **Per-file processing** (`process_vsi_file`):
   - Loads `.vsi` via `BioImage.get_image_data("CZYX", T=0)` → 4D array `(channels, z, y, x)`
   - Extracts DAPI stack from configured channel index
   - Segments nuclei in 3D (`segment_nuclei_3d`)
   - Saves per-slice label images (`save_label_images`)
   - Extracts per-nucleus intensity metrics (`extract_intensity_metrics`)
   - Tags results with filename and condition (from `get_condition_from_filename`)
4. **Aggregation**: Concatenates all per-file DataFrames into `output/nuclei_measurements.csv`
5. **Summary statistics** (`summarize_by_condition`): Computes per-condition mean/std → `output/summary_statistics.csv`
6. **Visualization** (`plot_intensity_summary`): Per-channel boxplot+stripplot with log y-axis, normalized to DAPI → `output/intensity_summary.png`

### Key implementation details

**Segmentation** (`segment_nuclei_3d`):
- Gaussian smoothing with anisotropic sigma `[1.0, 2.0, 2.0]` (z vs xy)
- Triangle thresholding (not Otsu — Otsu under-segmented)
- Binary fill holes, erosion/dilation with ball(2)
- Connected component labeling in 3D
- Size filtering via `np.bincount` + lookup-table relabel (vectorized, not per-object scan)

**Intensity extraction** (`extract_intensity_metrics`):
- Uses `scipy.ndimage` functions (`mean`, `median`, `minimum`, `maximum`, `standard_deviation`, `sum_labels`) with labeled arrays
- Computes all stats for all nuclei at once per channel (vectorized, not per-nucleus loop)

**Condition mapping** (`get_condition_from_filename`):
- Parses leading numeric index from filenames like `10_Multichannel Z-Stack_20260622_67.vsi`
- Looks up index in `--condition-mapping` JSON object
- Returns "Unknown" if no match

**Plotting** (`plot_intensity_summary`):
- Normalizes each non-DAPI channel's mean intensity to that nucleus's `DAPI_mean`
- Log-scaled y-axis (intensities span wide range)
- Boxplot + stripplot (jittered points, not swarm — swarm drops points in dense conditions)
- DAPI excluded from panels (only used as normalization reference)

## Output structure

All outputs go to `./output/` (created automatically, not committed):
- `nuclei_measurements.csv` — per-nucleus metrics (one row per nucleus)
- `summary_statistics.csv` — per-condition summary (one row per condition)
- `intensity_summary.png` — visualization grid
- `label_images/<filename_stem>/z###.png` — per-slice overlays for segmentation quality checks

## Input data expectations

- Files must be `.vsi` format (Olympus)
- Filenames should start with a numeric index for condition mapping: `<index>_...`
- Reading `.vsi` goes through `bioio-bioformats`, which requires a JVM (downloads via `cjdk` on first run if needed)

## Important conventions

- **No hardcoded paths**: Use `DATA_DIR`, `OUTPUT_DIR`, `LABEL_IMAGE_DIR` module constants or CLI args
- **Channel order matters**: `--channel-names` position = channel index in the acquired image
- **DAPI is required**: Must be present in `--channel-names`; used for segmentation
- **Vectorized operations preferred**: The codebase avoids per-object/per-nucleus loops in favor of scipy.ndimage vectorized functions
- **Module-level defaults**: Constants near top of file provide CLI defaults; override via CLI rather than editing constants

## Gotchas

1. **`pixi.lock` is gitignored**: Always commit changes to `pixi.toml` — the lockfile won't be tracked
2. **JVM requirement**: First run may be slow due to JRE download for `bioio-bioformats`
3. **Filename parsing is fixed**: `get_condition_from_filename` expects `<index>_Multichannel Z-Stack_<date>_<n>.vsi` format. If your filenames differ, edit this function — there's no CLI option for the parsing pattern
4. **Size filtering formula**: The min/max voxel count bounds use an approximate spherical volume formula with arbitrary scaling factors (`/10` and `*10`). Adjust `--nuclei-diameter-px` and `--size-tolerance` if segmentation misses expected nuclei
5. **No error recovery**: If a file fails to process, it's skipped with a print message. Check console output for errors
6. **Memory usage**: Loading full 4D stacks into memory — large datasets may need adjustment

## Testing

No test suite exists. To verify changes:
1. Run the pipeline on sample data
2. Check console output for errors
3. Inspect `label_images/` for segmentation quality
4. Verify CSV outputs have expected columns and rows
5. Check plot renders correctly

## Common modifications

**Adding a new channel**: Add name to `--channel-names` in correct acquisition order position
**Changing nucleus size**: Adjust `--nuclei-diameter-px` and/or `--size-tolerance`
**Different experimental conditions**: Update `--condition-mapping` JSON
**Different filename format**: Edit `get_condition_from_filename` function
**Different segmentation parameters**: Modify `segment_nuclei_3d` (smoothing sigma, threshold method, morphology operations)
