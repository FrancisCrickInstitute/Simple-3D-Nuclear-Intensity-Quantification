import argparse
import json
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bioio import BioImage
from scipy import ndimage
from scipy.ndimage import label, binary_fill_holes
from skimage import filters, morphology
from skimage.color import label2rgb

# Configuration defaults (overridable via CLI, see parse_args)
DATA_DIR = "./data"
OUTPUT_DIR = "./output"
CHANNEL_NAMES = ["DAPI", "HA", "CPSF6", "Capsid"]
NUCLEI_DIAMETER_PX = 140
SIZE_TOLERANCE = 0.3
LABEL_IMAGE_DIR = "./output/label_images"

# Default experimental conditions based on file indices
CONDITION_MAPPING = {
    7: "D37_RR-VLPs", 10: "D37_RR-VLPs", 11: "D37_RR-VLPs", 12: "D37_RR-VLPs", 13: "D37_RR-VLPs",
    14: "D102-VLPs", 15: "D102-VLPs", 16: "D102-VLPs", 22: "D102-VLPs",
    17: "Uninfected", 18: "Uninfected", 19: "Uninfected", 20: "Uninfected", 21: "Uninfected"
}

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)


def parse_condition_mapping(s):
    """Parse a JSON object string mapping file index -> condition label."""
    raw = json.loads(s)
    return {int(k): v for k, v in raw.items()}


def get_condition_from_filename(filename, condition_mapping):
    """Extract image index from filename and return condition."""
    # Filenames are formatted like "10_Multichannel Z-Stack_20260622_67.vsi",
    # where the leading number is the file index used in condition_mapping.
    numbers = re.findall(r"\d+", Path(filename).stem)
    if not numbers:
        return "Unknown"
    return condition_mapping.get(int(numbers[-1]), "Unknown")


def segment_nuclei_3d(dapi_stack, nuclei_diameter_px, size_tolerance):
    """
    Segment nuclei from 3D DAPI z-stack.

    Parameters:
    dapi_stack: 3D numpy array (z, y, x)
    nuclei_diameter_px: expected nucleus diameter in pixels
    size_tolerance: acceptable fractional deviation from nuclei_diameter_px

    Returns:
    labeled_3d: 3D labeled image with nucleus IDs
    """

    # Apply Gaussian smoothing to reduce noise
    dapi_smooth = ndimage.gaussian_filter(dapi_stack, sigma=[1.0, 2.0, 2.0])

    # Create initial binary mask using the triangle threshold
    threshold = filters.threshold_triangle(dapi_smooth)
    binary_mask = dapi_smooth > threshold

    # Fill holes in the binary mask
    binary_mask = binary_fill_holes(binary_mask)

    # Apply morphological operations to clean up
    binary_mask = morphology.binary_erosion(binary_mask, morphology.ball(2))
    binary_mask = morphology.binary_dilation(binary_mask, morphology.ball(2))

    # Label connected components in 3D
    labeled_3d, num_features = label(binary_mask)

    # Filter by size: keep nuclei within acceptable size range
    min_size = int(np.pi * (nuclei_diameter_px * (1 - size_tolerance) / 2) ** 2 / 10)
    max_size = int(np.pi * (nuclei_diameter_px * (1 + size_tolerance) / 2) ** 2 * 10)

    # Voxel count per label, then relabel sequentially in one pass rather than
    # rescanning the full volume once per candidate nucleus.
    voxel_counts = np.bincount(labeled_3d.ravel(), minlength=num_features + 1)
    keep = (voxel_counts >= min_size) & (voxel_counts <= max_size)
    keep[0] = False  # background is never a nucleus

    new_labels = np.zeros(num_features + 1, dtype=labeled_3d.dtype)
    new_labels[keep] = np.arange(1, keep.sum() + 1)

    return new_labels[labeled_3d]


def save_label_images(dapi_stack, labeled_nuclei, output_subdir):
    """
    Save one PNG per z-slice showing the DAPI signal with segmented nuclei overlaid.

    Parameters:
    dapi_stack: 3D numpy array (z, y, x) of raw DAPI intensities
    labeled_nuclei: 3D labeled image with nucleus IDs, same shape as dapi_stack
    output_subdir: directory to save the per-slice PNGs into (created if needed)
    """
    os.makedirs(output_subdir, exist_ok=True)

    dapi_min = np.percentile(dapi_stack, 1)
    dapi_max = np.percentile(dapi_stack, 99)
    dapi_normalized = np.clip((dapi_stack - dapi_min) / (dapi_max - dapi_min), 0, 1)

    for z in range(dapi_stack.shape[0]):
        overlay = label2rgb(labeled_nuclei[z], image=dapi_normalized[z], bg_label=0, alpha=0.4)
        plt.imsave(os.path.join(output_subdir, f"z{z:03d}.png"), overlay)


def extract_intensity_metrics(image_data, labeled_nuclei, nucleus_ids, channels):
    """
    Extract intensity metrics for every nucleus across all channels in one pass.

    Parameters:
    image_data: 4D array (channels, z, y, x)
    labeled_nuclei: 3D labeled image
    nucleus_ids: array-like of nucleus label IDs to measure
    channels: list of (channel_name, channel_index) pairs to measure

    Returns:
    pandas DataFrame with one row per nucleus and per-channel intensity metrics
    """
    metrics = {"nucleus_id": np.asarray(nucleus_ids)}

    # Compute each stat for all nuclei at once with scipy.ndimage, rather than
    # rescanning the full volume per nucleus per channel.
    for ch_name, ch_idx in channels:
        if ch_idx >= image_data.shape[0]:
            continue
        channel_data = image_data[ch_idx]
        metrics[f"{ch_name}_mean"] = ndimage.mean(channel_data, labeled_nuclei, nucleus_ids)
        metrics[f"{ch_name}_median"] = ndimage.median(channel_data, labeled_nuclei, nucleus_ids)
        metrics[f"{ch_name}_min"] = ndimage.minimum(channel_data, labeled_nuclei, nucleus_ids)
        metrics[f"{ch_name}_max"] = ndimage.maximum(channel_data, labeled_nuclei, nucleus_ids)
        metrics[f"{ch_name}_std"] = ndimage.standard_deviation(channel_data, labeled_nuclei, nucleus_ids)
        metrics[f"{ch_name}_total"] = ndimage.sum_labels(channel_data, labeled_nuclei, nucleus_ids)

    return pd.DataFrame(metrics)


def process_vsi_file(filepath, dapi_channel, channels, nuclei_diameter_px, size_tolerance, condition_mapping):
    """
    Process a single VSI file: segment nuclei and extract intensity metrics.

    Parameters:
    filepath: Path to VSI file
    dapi_channel: channel index of the DAPI (nuclear) stain, used for segmentation
    channels: list of (channel_name, channel_index) pairs to measure
    nuclei_diameter_px: expected nucleus diameter in pixels
    size_tolerance: acceptable fractional deviation from nuclei_diameter_px
    condition_mapping: maps the numeric file index parsed from each filename to a condition label

    Returns:
    pandas DataFrame with per-nucleus measurements
    """
    print(f"Processing: {filepath}")

    try:
        # Read image using bioio, collapsing the timepoint axis to get (C, Z, Y, X)
        bio_image = BioImage(filepath)
        image_data = bio_image.get_image_data("CZYX", T=0)
        dapi_stack = image_data[dapi_channel]

        # Segment nuclei in 3D
        labeled_nuclei = segment_nuclei_3d(dapi_stack, nuclei_diameter_px, size_tolerance)
        num_nuclei = int(labeled_nuclei.max())

        print(f"  Found {num_nuclei} nuclei")

        # Save per-slice label images for visual inspection
        filename = Path(filepath).name
        label_image_dir = os.path.join(LABEL_IMAGE_DIR, Path(filepath).stem)
        save_label_images(dapi_stack, labeled_nuclei, label_image_dir)

        # Extract metrics for all nuclei at once
        measurements = extract_intensity_metrics(
            image_data, labeled_nuclei, np.arange(1, num_nuclei + 1), channels
        )
        measurements["filename"] = filename
        measurements["condition"] = get_condition_from_filename(filename, condition_mapping)

        return measurements

    except Exception as e:
        print(f"  Error processing {filepath}: {e}")
        return None


def summarize_by_condition(results_df):
    """
    Compute per-condition summary statistics (mean/std of each intensity metric).

    Parameters:
    results_df: per-nucleus measurements DataFrame from process_vsi_file

    Returns:
    pandas DataFrame with one row per condition
    """
    summary_stats = []

    for condition in results_df["condition"].unique():
        condition_data = results_df[results_df["condition"] == condition]

        summary = {
            "condition": condition,
            "num_images": condition_data["filename"].nunique(),
            "num_nuclei": len(condition_data),
        }

        # Add mean values for each channel metric
        for col in condition_data.columns:
            if col not in ["nucleus_id", "filename", "condition"]:
                summary[f"{col}_mean"] = condition_data[col].mean()
                summary[f"{col}_std"] = condition_data[col].std()

        summary_stats.append(summary)

    return pd.DataFrame(summary_stats)


def plot_intensity_summary(results_df, plot_file, channels):
    """
    Save a per-channel boxplot with individual nuclei overlaid as a strip plot.

    Mean intensities are normalized to each nucleus's own DAPI mean intensity
    before plotting, and the y-axis is log-scaled, since raw intensities span
    a wide range across channels/conditions. DAPI itself is excluded from the
    plotted channels, since it's only used here as the normalization
    reference (its own panel would just be a ~1.0 flat line).

    Parameters:
    results_df: per-nucleus measurements DataFrame from process_vsi_file
    plot_file: path to save the PNG to
    channels: list of (channel_name, channel_index) pairs that were measured
    """
    plot_df = results_df.copy()
    plot_channels = [(name, idx) for name, idx in channels if name != "DAPI"]
    for channel, _ in plot_channels:
        plot_df[f"{channel}_mean_norm"] = plot_df[f"{channel}_mean"] / plot_df["DAPI_mean"]

    ncols = 2
    nrows = int(np.ceil(len(plot_channels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False)
    axes = axes.flatten()
    fig.suptitle("Nuclei Intensity Analysis Summary", fontsize=16)

    condition_order = sorted(plot_df["condition"].unique())

    for idx, (channel, _) in enumerate(plot_channels):
        ax = axes[idx]
        col = f"{channel}_mean_norm"

        sns.boxplot(
            data=plot_df, x="condition", y=col, order=condition_order,
            ax=ax, showfliers=False, color="lightgray",
        )
        sns.stripplot(
            data=plot_df, x="condition", y=col, order=condition_order,
            ax=ax, size=1.0, color="black", alpha=0.6, jitter=True,
        )
        ax.set_yscale("log")
        ax.set_title(channel)
        ax.set_xlabel("")
        ax.set_ylabel("Mean intensity per nucleus\n(normalized to DAPI)")
        ax.tick_params(axis="x", rotation=30)

    for ax in axes[len(plot_channels):]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(plot_file, dpi=150)
    plt.close(fig)


def main(data_dir, channel_names, nuclei_diameter_px, size_tolerance, condition_mapping):
    """Main analysis pipeline."""

    if "DAPI" not in channel_names:
        print("--channel-names must include 'DAPI' (used for nucleus segmentation)")
        return

    channels = [(name, idx) for idx, name in enumerate(channel_names)]
    dapi_channel = channel_names.index("DAPI")

    # Find all VSI files in data directory
    vsi_files = list(Path(data_dir).glob("*.vsi"))

    if not vsi_files:
        print(f"No VSI files found in {data_dir}")
        return

    print(f"Found {len(vsi_files)} VSI files")

    # Process each file and collect results
    all_measurements = []

    for vsi_file in sorted(vsi_files):
        df = process_vsi_file(
            str(vsi_file), dapi_channel, channels, nuclei_diameter_px, size_tolerance, condition_mapping
        )
        if df is not None and len(df) > 0:
            all_measurements.append(df)

    if not all_measurements:
        return

    results_df = pd.concat(all_measurements, ignore_index=True)

    # Save full results
    output_file = os.path.join(OUTPUT_DIR, "nuclei_measurements.csv")
    results_df.to_csv(output_file, index=False)
    print(f"\nSaved full measurements to: {output_file}")
    print(f"Total nuclei measured: {len(results_df)}")

    # Save summary statistics by condition
    summary_df = summarize_by_condition(results_df)
    summary_file = os.path.join(OUTPUT_DIR, "summary_statistics.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved summary statistics to: {summary_file}")

    print("\nSummary by condition:")
    print(summary_df[["condition", "num_images", "num_nuclei"]])

    # Save per-channel boxplot + swarm plot visualization
    plot_file = os.path.join(OUTPUT_DIR, "intensity_summary.png")
    plot_intensity_summary(results_df, plot_file, channels)
    print(f"Saved summary plot to: {plot_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Segment nuclei and quantify per-channel intensity from VSI z-stacks."
    )
    parser.add_argument(
        "--data-dir", default=DATA_DIR,
        help=f"Directory containing .vsi files (default: {DATA_DIR})",
    )
    parser.add_argument(
        "--channel-names", nargs="+", metavar="NAME", default=CHANNEL_NAMES,
        help="Ordered list of channel names, one per channel index in the acquired image, e.g. "
             "--channel-names DAPI HA CPSF6 Capsid. Must include 'DAPI', which is used for nucleus "
             f"segmentation (default: {' '.join(CHANNEL_NAMES)})",
    )
    parser.add_argument(
        "--nuclei-diameter-px", type=int, default=NUCLEI_DIAMETER_PX,
        help=f"Expected nucleus diameter in pixels, used to filter segmented objects by size "
             f"(default: {NUCLEI_DIAMETER_PX})",
    )
    parser.add_argument(
        "--size-tolerance", type=float, default=SIZE_TOLERANCE,
        help=f"Acceptable fractional deviation from --nuclei-diameter-px (default: {SIZE_TOLERANCE})",
    )
    parser.add_argument(
        "--condition-mapping", type=parse_condition_mapping, default=CONDITION_MAPPING,
        help="JSON object mapping the numeric file index parsed from each filename to an experimental "
             'condition label, e.g. \'{"1": "ConditionA", "2": "ConditionB"}\' '
             "(default: built-in HIV-Quant mapping)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.data_dir, args.channel_names, args.nuclei_diameter_px, args.size_tolerance, args.condition_mapping,
    )