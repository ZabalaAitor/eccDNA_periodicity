from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde


# Length range used for the analysis and plots
len_min = 0
len_max = 1000


# ============================================================
# KDE peak analysis
# ============================================================

def analyze_kde_peaks(
    lengths,
    bandwidth_method="silverman",
    x_resolution=1000,
):
    """
    Fit a KDE to a length distribution and detect its peaks.

    Parameters
    ----------
    lengths : array-like
        Circle length measurements.

    bandwidth_method : str or scalar
        Bandwidth method passed to scipy.stats.gaussian_kde.

    x_resolution : int
        Number of points used to evaluate the KDE.

    Returns
    -------
    peak_x : np.ndarray
        Positions of the detected KDE peaks.

    mean_peak_difference : float
        Mean distance between consecutive peaks.
        Returns NaN when fewer than two peaks are detected.
    """

    lengths = np.asarray(lengths)

    if lengths.size == 0:
        return np.array([]), np.nan

    # Fit the kernel density estimate
    kde = gaussian_kde(
        lengths,
        bw_method=bandwidth_method
    )

    # Define the x-axis over the observed length range
    x_values = np.linspace(
        lengths.min(),
        lengths.max(),
        x_resolution
    )

    # Evaluate the KDE
    y_values = kde(x_values)

    # Detect local maxima in the KDE
    peaks, _ = find_peaks(
        y_values
    )

    # Convert peak positions to integer base-pair values
    peak_x = x_values[peaks].round(0)

    # Calculate the mean spacing between consecutive peaks
    mean_peak_difference = (
        np.mean(np.diff(peak_x))
        if peak_x.size > 1
        else np.nan
    )

    return peak_x, mean_peak_difference


# ============================================================
# Length distribution plotting and KDE analysis
# ============================================================

def plot_length_distributions(
    dataset_group,
    output_file,
    peak_window=100,
    verbose=False,
    kde_output_file=None,
    bandwidth_method="silverman",
    x_resolution=1000,
):
    """
    Plot eccDNA length distributions and perform KDE peak analysis.

    Parameters
    ----------
    dataset_group : list
        List of (dataset_name, lengths) tuples.

    output_file : str or Path
        Path for the output PNG file. A PDF is also created.

    peak_window : int
        Width of the zoomed region centered on the main histogram peak.

    verbose : bool
        If True, each dataset occupies one column with the full
        distribution on the top row and the zoomed distribution
        on the bottom row.

        If False, each dataset occupies two panels:
        full distribution on the left and zoomed distribution
        on the right.

    kde_output_file : str or Path, optional
        Path for the KDE peak analysis CSV.

        If None, the CSV is saved next to output_file using:
        <output_stem>_kde_peaks.csv

    bandwidth_method : str or scalar
        Bandwidth method passed to scipy.stats.gaussian_kde.

    x_resolution : int
        Number of points used to evaluate the KDE.

    Returns
    -------
    pd.DataFrame
        KDE peak analysis results for each dataset.
    """

    dataset_group = list(
        dataset_group
    )

    num_datasets = len(
        dataset_group
    )

    if num_datasets == 0:
        return pd.DataFrame()

    # ========================================================
    # Figure layout
    # ========================================================

    if verbose:
        nrows = 2
        ncols = num_datasets
        figsize = (4 * ncols, 8)

    else:
        nrows = (num_datasets + 1) // 2
        ncols = 4
        figsize = (16, 3 * nrows)

    fig, axs = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        sharey=False
    )

    # Keep the axes array consistently two-dimensional
    axs = np.array(axs).reshape(
        nrows,
        ncols
    )

    # ========================================================
    # Histogram calculation
    # ========================================================

    def compute_smoothed_histogram(
        lengths,
        smooth=True
    ):
        """
        Calculate a histogram and optionally smooth it.

        The histogram range follows the original plotting
        implementation. Sigma=3 is used for the full distribution
        and sigma=1 for the zoomed distribution.
        """

        bins = np.arange(
            0,
            lengths.max() + 2
        )

        hist, bin_edges = np.histogram(
            lengths,
            bins=bins
        )

        bin_centers = (
            bin_edges[:-1] +
            bin_edges[1:]
        ) / 2

        sigma = 3 if smooth else 1

        smoothed_hist = gaussian_filter1d(
            hist,
            sigma=sigma
        )

        return (
            bin_centers,
            smoothed_hist
        )

    # Store KDE results for all datasets
    kde_results = []

    # ========================================================
    # Dataset plotting loop
    # ========================================================

    for i, (dataset_name, lengths) in enumerate(
        dataset_group
    ):

        lengths = np.asarray(
            lengths
        )

        if lengths.size == 0:
            continue

        # ====================================================
        # Determine subplot positions
        # ====================================================

        if verbose:
            # Full distribution on the top row
            # Zoomed distribution on the bottom row
            row = 0
            col_left = i
            col_right = i
            zoom_row = 1

        else:
            # Two panels per dataset:
            # left = full distribution
            # right = zoomed distribution
            row = i // 2
            col_left = (i % 2) * 2
            col_right = col_left + 1

        # ====================================================
        # Full length distribution
        # ====================================================

        bin_centers, smoothed_hist = (
            compute_smoothed_histogram(
                lengths,
                smooth=True
            )
        )

        # Number of circles within the selected length range
        full_count = (
            (lengths >= len_min) &
            (lengths <= len_max)
        ).sum()

        # Plot the full distribution in grey
        axs[
            row,
            col_left
        ].plot(
            bin_centers,
            smoothed_hist,
            color="#adadad",
            lw=2
        )

        axs[
            row,
            col_left
        ].set_xlim(
            len_min,
            len_max
        )

        axs[
            row,
            col_left
        ].tick_params(
            axis="both",
            labelsize=16
        )

        sns.despine(
            ax=axs[
                row,
                col_left
            ],
            top=True,
            right=True
        )

        # Display the dataset size
        if verbose:
            axs[
                row,
                col_left
            ].set_title(
                f"{dataset_name}\nn = {full_count}",
                fontsize=18
            )

        else:
            axs[
                row,
                col_left
            ].set_title(
                f"n = {full_count}",
                fontsize=14
            )

        # ====================================================
        # Main histogram peak
        # ====================================================

        peaks, _ = find_peaks(
            smoothed_hist
        )

        if len(peaks) > 0:

            # Select the highest detected histogram peak
            main_peak_idx = peaks[
                np.argmax(
                    smoothed_hist[peaks]
                )
            ]

            main_peak_pos = (
                bin_centers[
                    main_peak_idx
                ]
            )

            # Define the zoom window around the main peak
            range_min = max(
                main_peak_pos -
                peak_window // 2,
                0
            )

            range_max = (
                main_peak_pos +
                peak_window // 2
            )

        else:

            # Default zoom range when no peak is detected
            range_min = 130
            range_max = 230

        # ====================================================
        # Highlight the zoom region
        # ====================================================

        mask = (
            (bin_centers >= range_min) &
            (bin_centers <= range_max)
        )

        axs[
            row,
            col_left
        ].plot(
            bin_centers[mask],
            smoothed_hist[mask],
            color="black",
            lw=2
        )

        # Number of circles inside the zoom window
        zoom_count = (
            (lengths >= range_min) &
            (lengths <= range_max)
        ).sum()

        # ====================================================
        # Zoomed length distribution
        # ====================================================

        # Recalculate the histogram with minimal smoothing
        # to preserve the original zoomed-view appearance.
        zoom_bin_centers, zoom_hist = (
            compute_smoothed_histogram(
                lengths,
                smooth=False
            )
        )

        if verbose:

            zoom_ax = axs[
                zoom_row,
                col_right
            ]

            zoom_ax.plot(
                zoom_bin_centers,
                zoom_hist,
                color="black",
                lw=2
            )

            zoom_ax.set_xlim(
                range_min,
                range_max
            )

            zoom_ax.tick_params(
                axis="both",
                labelsize=16
            )

            sns.despine(
                ax=zoom_ax,
                top=True,
                right=True
            )

            zoom_ax.set_title(
                f"n = {zoom_count}",
                fontsize=18
            )

        else:

            zoom_ax = axs[
                row,
                col_right
            ]

            zoom_ax.plot(
                zoom_bin_centers,
                zoom_hist,
                color="black",
                lw=2
            )

            zoom_ax.set_xlim(
                range_min,
                range_max
            )

            zoom_ax.tick_params(
                axis="both",
                labelsize=16
            )

            sns.despine(
                ax=zoom_ax,
                top=True,
                right=True
            )

            zoom_ax.set_title(
                f"n = {zoom_count}",
                fontsize=14
            )

            zoom_ax.set_ylabel("")

        # ====================================================
        # Axis labels
        # ====================================================

        if verbose:

            axs[
                1,
                i
            ].set_xlabel(
                "Length (bp)",
                fontsize=18
            )

            if i == 0:

                axs[
                    0,
                    0
                ].set_ylabel(
                    "Frequency",
                    fontsize=18
                )

                axs[
                    1,
                    0
                ].set_ylabel(
                    "Frequency",
                    fontsize=18
                )

        else:

            # Add x-axis labels to the bottom row
            if i >= num_datasets - 2:

                axs[
                    row,
                    col_left
                ].set_xlabel(
                    "Length (bp)",
                    fontsize=18
                )

                axs[
                    row,
                    col_right
                ].set_xlabel(
                    "Length (bp)",
                    fontsize=18
                )

        # ====================================================
        # Dataset name between the two panels
        # ====================================================

        if not verbose:

            left_pos = axs[
                row,
                col_left
            ].get_position().x0

            right_pos = axs[
                row,
                col_right
            ].get_position().x1

            if i % 2 == 0:

                center_x = (
                    left_pos +
                    right_pos
                ) / 2

                top_y = (
                    0.95 -
                    i * 0.0378 +
                    0.03
                )

            else:

                center_x = (
                    left_pos +
                    right_pos
                ) / 2 + 0.05

                top_y = (
                    0.95 -
                    (i - 1) * 0.0378 +
                    0.03
                )

            fig.text(
                center_x,
                top_y,
                dataset_name,
                ha="center",
                va="bottom",
                fontsize=20
            )

        # ====================================================
        # KDE peak analysis
        # ====================================================

        peak_x, mean_peak_difference = (
            analyze_kde_peaks(
                lengths,
                bandwidth_method=bandwidth_method,
                x_resolution=x_resolution
            )
        )

        kde_results.append({
            "Dataset": dataset_name,
            "Peaks": peak_x.tolist(),
            "Main Peaks": peak_x[:2].tolist(),
            "Other Peaks": peak_x[2:].tolist(),
            "Mean Peak Difference": mean_peak_difference,
        })

    # ========================================================
    # Hide unused subplot axes
    # ========================================================

    if not verbose:

        total_dataset_slots = (
            nrows * 2
        )

        for j in range(
            num_datasets,
            total_dataset_slots
        ):

            row = j // 2

            col_left = (
                (j % 2) * 2
            )

            col_right = (
                col_left + 1
            )

            axs[
                row,
                col_left
            ].axis("off")

            axs[
                row,
                col_right
            ].axis("off")

    # ========================================================
    # Final figure layout
    # ========================================================

    plt.tight_layout(
        rect=[
            0,
            0,
            1,
            0.98
        ]
    )

    space = (
        0.3
        if verbose
        else 0.75
    )

    plt.subplots_adjust(
        hspace=space
    )

    # ========================================================
    # Save figure
    # ========================================================

    output_file = Path(
        output_file
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        output_file.with_suffix(".pdf"),
        bbox_inches="tight"
    )

    # ========================================================
    # Save KDE results
    # ========================================================

    kde_results = pd.DataFrame(
        kde_results
    )

    if kde_output_file is None:

        kde_output_file = (
            output_file.parent /
            f"{output_file.stem}_kde_peaks.csv"
        )

    else:

        kde_output_file = Path(
            kde_output_file
        )

    kde_output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    kde_results.to_csv(
        kde_output_file,
        index=False
    )

    # Display and close the figure
    plt.show()
    plt.close(fig)

    return kde_results