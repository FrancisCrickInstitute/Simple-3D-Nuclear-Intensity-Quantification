# Simple 3D Nuclear Quantification

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.14](https://img.shields.io/badge/Python-3.14-green.svg)](https://www.python.org/downloads/)
[![Managed by Pixi](https://img.shields.io/badge/managed%20by-pixi-yellow.svg)](https://pixi.sh)
[![Platforms](https://img.shields.io/badge/platforms-linux--64%20%7C%20win--64-lightgrey.svg)](https://pixi.sh)
![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/HIV-Quant?style=plastic)

A general-purpose bioimage analysis pipeline for segmenting nuclei in 3D and quantifying per-channel intensity
within them, from multi-channel confocal z-stacks. Input is read via BioImage, so it accepts raw Olympus `.vsi`,
TIFF/OME-TIFF, CZI, and other common microscopy formats. Nuclei are segmented from a DAPI channel, then
per-nucleus intensity statistics are measured across the remaining channels, aggregated by experimental condition,
and written out as CSVs, a summary plot, and per-slice label image overlays.

This repo was originally built to quantify HIV capsid/CPSF6/HA intensity in infected-cell nuclei, but the
segmentation and measurement logic is not HIV-specific — see [Adapting to other experiments](#adapting-to-other-experiments)
below to reuse it for any experiment that needs per-nucleus intensity quantification across channels.

---

## Quick start for beginners

If you're new to Python, Jupyter, or command-line tools, follow these steps to get up and running on your own laptop (Windows, Mac, or Linux).

### Step 1: Install Git

You'll need Git to download this repository.

**Windows:**
1. Download Git from https://git-scm.com/download/win
2. Run the installer and accept all defaults

**Mac:**
1. Open Terminal (search "Terminal" in Spotlight)
2. Type `git --version` and press Enter
3. If not installed, it will prompt you to install Xcode Command Line Tools — follow the prompts

**Linux:**
```bash
sudo apt install git   # Ubuntu/Debian
# or
sudo yum install git   # CentOS/RHEL
```

### Step 2: Clone this repository

Open a terminal (Command Prompt on Windows, Terminal on Mac/Linux) and run:

```bash
git clone https://github.com/FrancisCrickInstitute/Simple-3D-Nuclear-Intensity-Quantification.git
cd Simple-3D-Nuclear-Intensity-Quantification
```

This downloads all the code into a folder called `Simple-3D-Nuclear-Intensity-Quantification`.

### Step 3: Install Pixi

Pixi manages the Python environment and dependencies for this project.

**All platforms:**
1. Visit https://pixi.sh/latest/get_started/
2. Follow the installation instructions for your operating system
3. Verify installation by opening a new terminal and typing:
   ```bash
   pixi --version
   ```

### Step 4: Set up the environment

From inside the repository folder (where you ran `cd` above), run:

```bash
pixi install
```

This creates a virtual environment with all required packages (numpy, pandas, scipy, bioio, etc.). It may take a few minutes on first run.

### Step 5: Find the path to your data

No need to copy files — you can leave your `.vsi` images wherever they are already stored. You just need to know the full folder path.

**Finding the path:**
- **Windows:** Open File Explorer, navigate to your data folder, click the address bar at the top, and copy the full path (e.g. `C:\Users\YourName\Documents\MyMicroscopeImages`)
- **Mac:** Right-click the folder, hold **Option**, and click "Copy [foldername] as Pathname"
- **Linux:** In a file browser, right-click and "Copy Location" or navigate there in a terminal and run `pwd`

You'll use this path in the next step. Each image's filename ends with a unique index used to look up its
condition, e.g. `p34_EXP1_11_Multichannel Z-Stack_20260805_140.vsi` ends in `140`. That trailing number is
mapped to an experimental condition via `CONDITION_MAPPING`.

### Step 6: Run the analysis

You have two options:

#### Option A: Use the Jupyter Notebook (recommended for beginners)

The notebook walks through the analysis step-by-step with visual feedback at each stage.

1. Start Jupyter:
   ```bash
   pixi run jupyter notebook
   ```
2. Your browser should open automatically. Click on `quantify_nuclei_intensity.ipynb`
3. In the **Configuration** cell (near the top), change `DATA_DIR` to point at your data folder — use the path you found in Step 5, with quotes around it:
   ```python
   DATA_DIR = "C:\\Users\\YourName\\Documents\\MyMicroscopeImages"
   ```
   (On Mac/Linux use forward slashes: `"/home/yourname/MyMicroscopeImages"`)
4. Read through the notebook cells — they explain what's happening at each step. Edit `CHANNEL_NAMES` and `CONDITION_MAPPING` as needed for your experiment.
5. To run the analysis:
   - Click on a cell (it will highlight with a blue border)
   - Press **Shift+Enter** to run that cell and move to the next one
   - Or click **Cell → Run All** from the menu to run everything at once
6. Results will appear in the notebook and be saved to the `output` folder

#### Option B: Run the script directly

For a quick, non-interactive run (point it at your data folder instead of `./data`):

```bash
pixi run python quantify_nuclei_intensity.py --data-dir "C:\Users\YourName\Documents\MyMicroscopeImages"
```

On Mac/Linux:
```bash
pixi run python quantify_nuclei_intensity.py --data-dir "/home/yourname/MyMicroscopeImages"
```

See [Usage](#usage) below for all customization options.

### Step 7: Check your results

After the analysis completes, look in the `output` folder:

- `nuclei_measurements.csv` — detailed measurements for every nucleus found
- `summary_statistics.csv` — summary table grouped by experimental condition
- `intensity_summary.png` — visualization showing intensity differences between conditions
- `label_images/` — folder with overlay images showing where nuclei were detected (useful for quality checking)

---

## Requirements

This project uses [pixi](https://pixi.sh) for environment management. All dependencies (numpy, pandas, scipy,
bioio, matplotlib, scikit-image, seaborn, etc.) are declared in `pixi.toml`.

```
pixi install
```

Reading `.vsi` (Olympus) files goes through `bioio-bioformats`, which relies on a JVM. On first run it will
download a JRE via `cjdk` if one isn't already available.

## Usage

```
pixi run python quantify_nuclei_intensity.py [OPTIONS]
```

Run with `--help` for the full option list. Available options (defaults match the original HIV experiment this
pipeline was built for):

| Option | Default | Description |
| --- | --- | --- |
| `--data-dir` | `./data` | Directory containing image files |
| `--channel-names` | `DAPI HA CPSF6 Capsid` | Ordered list of channel names, one per channel index in the acquired image. Must include `DAPI`, which is used for nucleus segmentation |
| `--min-nuclei-diameter-px` | `98` | Smallest nucleus diameter in pixels to keep, used to filter segmented objects by size |
| `--max-nuclei-diameter-px` | `182` | Largest nucleus diameter in pixels to keep, used to filter segmented objects by size |
| `--condition-mapping` | built-in HIV-Quant mapping | JSON object mapping the numeric file index parsed from each filename to an experimental condition label, e.g. `'{"1": "ConditionA", "2": "ConditionB"}'` |

Filenames are expected in the form `p<plate>_EXP<exp>_<group>_Multichannel Z-Stack_<date>_<id>` (any
extension), e.g. `p34_EXP1_11_Multichannel Z-Stack_20260805_140.vsi`. The trailing per-image index (`140`) is
looked up in `--condition-mapping` to assign each image to an experimental condition. Input files may be any
format BioImage can read (`.vsi`, `.tif`/`.ome.tiff`, `.czi`, `.lif`, `.nd2`, `.zarr`, `.oir`); the pipeline
discovers files by the `IMAGE_EXTENSIONS` list in `quantify_nuclei_intensity.py`, which you can extend if your
images use a different extension.

## Output

Results are written to `./output/`:

- `nuclei_measurements.csv` — per-nucleus intensity metrics (mean/median/min/max/std/total) for each channel
- `summary_statistics.csv` — per-condition mean/std of each metric
- `intensity_summary.png` — per-channel boxplot with individual nuclei overlaid as a swarm plot
- `label_images/<filename_stem>/z###.png` — one PNG per z-slice per input file, showing the DAPI signal with
  segmented nuclei overlaid, for visually checking segmentation quality

## Configuration

Channel names/order, expected nucleus size, and the condition mapping are all CLI options — see the table in
[Usage](#usage) — so no source changes are needed to point the pipeline at a different experiment. The
filename-parsing convention in `get_condition_from_filename` is still fixed in the script; see below if that needs
to change too.

## Adapting to other experiments

Nothing about the segmentation or measurement code is specific to HIV, capsid, CPSF6, or HA — those are just the
default `--channel-names`. To reuse the pipeline for a different multi-channel experiment:

- Pass `--channel-names` with your own stain names in acquisition order (one per channel index), e.g.
  `--channel-names DAPI GFP mCherry`. The name at each position becomes the column/plot label for that channel, and
  whichever position is named `DAPI` (required) is used for nucleus segmentation. Any number of channels is fine.
- Use `--min-nuclei-diameter-px`/`--max-nuclei-diameter-px` to match your expected nucleus size range, or adjust
  the thresholding/morphology steps in `segment_nuclei_3d` if your nuclear stain behaves differently.
- Use `--condition-mapping` to match your own experimental groups, or edit `get_condition_from_filename` if
  conditions aren't identified by a trailing numeric index in the filename.

Everything downstream (per-nucleus metrics, per-condition summary, plots, label image overlays) works off those
config values and needs no further changes.

## Troubleshooting

**"Command not found: pixi"**
- Make sure you installed Pixi correctly (Step 3 above)
- Try closing and reopening your terminal
- Run `pixi --version` to verify it's working

**"No image files found"**
- Check that your `--data-dir` (or `DATA_DIR` in the notebook) points to the correct folder containing your images
- Check that files use an extension in `IMAGE_EXTENSIONS` (`.vsi`, `.tif`, `.czi`, etc.). If your format differs, add its extension to `IMAGE_EXTENSIONS` in `quantify_nuclei_intensity.py`

**Processing is very slow**
- If your data is stored on a network drive or server, transfer speeds can be a bottleneck — reading large 3D stacks over a slow connection can take a long time
- It may be worth copying some of your files to your local hard drive first and pointing `--data-dir` at the local copy.

**Jupyter won't open**
- After running `pixi run jupyter notebook`, look for a URL in the terminal output (starts with `http://localhost:8888`)
- Copy and paste that URL into your browser manually

**First run is very slow**
- The first time you read a file that needs a JVM (e.g. `.vsi`, `.czi`, `.lif`), the system downloads a Java runtime (takes 1-2 minutes)
- Subsequent runs will be much faster

**Out of memory errors**
- Large 3D image stacks use significant RAM
- Close other applications while running the analysis
- Consider processing fewer images at once
