"""Short-time Fourier transform of eccDNA length distributions.

Builds length-density signals, removes the low-frequency envelope with a
Butterworth high-pass filter, and computes STFT spectrograms used to resolve
the nucleosome-scale (~180 bp) and helical-scale (~10.5 bp) periodicities.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import stft, butter, filtfilt
import matplotlib.ticker as ticker
import seaborn as sns
from pathlib import Path


# STFT and filtering parameters
CUTOFF_FREQ = 0.05
FS = 1
NPERSEG = 400
NOVERLAP = NPERSEG - 1

# Frequency range displayed in the STFT plots
FREQ_MIN = 0
FREQ_MAX = 0.2

# Length range analyzed
LEN_MIN = 0
LEN_MAX = 1000


def highpass_filter(data, cutoff_freq, fs=FS, order=4):
    """
    Apply a high-pass Butterworth filter to remove low-frequency components.

    Parameters
    ----------
    data : array-like
        Input signal.
    cutoff_freq : float
        High-pass cutoff frequency.
    fs : float
        Sampling frequency.
    order : int
        Butterworth filter order.

    Returns
    -------
    filtered_data : array-like
        High-pass filtered signal.
    """
    nyquist = 0.5 * fs
    normal_cutoff = cutoff_freq / nyquist

    b, a = butter(
        order,
        normal_cutoff,
        btype="high",
        analog=False
    )

    return filtfilt(b, a, data)


def zscore(x):
    """Apply z-score normalization to an array."""
    return (x - np.mean(x)) / np.std(x)


def process_dataset_for_stft(
    data_dict,
    dataset,
    ax,
    cutoff_freq=CUTOFF_FREQ,
    fs=FS,
    nperseg=NPERSEG,
    noverlap=NOVERLAP,
    freq_min=FREQ_MIN,
    freq_max=FREQ_MAX,
    normalize=True,
    legend=True,
    verbose=False
):
    """
    Calculate and plot the STFT magnitude of an eccDNA length distribution.

    The length distribution is first converted into a histogram, followed
    by high-pass filtering and optional z-score normalization. The STFT
    magnitude is then displayed as a function of length and frequency.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the dataset length arrays.

    dataset : str
        Dataset name. The corresponding key is expected to be
        ``<dataset>_lengths``.

    ax : matplotlib.axes.Axes
        Axes on which the STFT is plotted.

    cutoff_freq : float
        High-pass filter cutoff frequency.

    fs : float
        Sampling frequency used for the STFT.

    nperseg : int
        Number of points per STFT segment.

    noverlap : int
        Number of overlapping points between consecutive STFT segments.

    freq_min : float
        Minimum frequency displayed.

    freq_max : float
        Maximum frequency displayed.

    normalize : bool
        Whether to apply z-score normalization before the STFT.

    legend : bool
        Whether to display the colorbar label.

    verbose : bool
        Whether to use the extended x-axis tick configuration.
    """

    # Retrieve the length measurements for the dataset
    lengths = data_dict[f"{dataset}_lengths"]

    # Convert the length measurements into a frequency distribution
    length_bins = np.arange(LEN_MIN, LEN_MAX)
    length_distribution, _ = np.histogram(
        lengths,
        bins=length_bins
    )

    # Remove low-frequency components from the length distribution
    filtered_distribution = highpass_filter(
        length_distribution,
        cutoff_freq=cutoff_freq,
        fs=fs
    )

    # Normalize the filtered signal when requested
    if normalize:
        filtered_distribution = zscore(
            filtered_distribution
        )

    # Calculate the short-time Fourier transform
    frequencies, time, Zxx = stft(
        filtered_distribution,
        fs=fs,
        nperseg=nperseg,
        noverlap=noverlap
    )

    # Convert the STFT time axis back to the original length scale
    time = time + LEN_MIN

    # Plot STFT magnitude
    pcm = ax.pcolormesh(
        time,
        frequencies,
        np.abs(Zxx),
        shading="gouraud",
        rasterized=True
    )

    # Add the magnitude colorbar
    cbar = plt.colorbar(
        pcm,
        ax=ax
    )

    cbar.ax.tick_params(
        labelsize=16
    )

    if legend:
        cbar.set_label(
            "Magnitude",
            fontsize=18
        )
    else:
        cbar.set_label(
            "",
            fontsize=18
        )

    # Configure the x-axis
    ax.set_xlabel(
        "Length (bp)",
        fontsize=18
    )

    ax.set_xlim(
        0,
        time[-1]
    )

    # Limit the displayed frequency range
    ax.set_ylim(
        freq_min,
        freq_max
    )

    # Use fewer x-axis ticks in the standard layout
    if verbose:
        ax.set_xticks(
            [0, 200, 400, 600, 800, 1000]
        )
    else:
        ax.set_xticks(
            [0, 500, 1000]
        )

    ax.tick_params(
        axis="x",
        labelsize=16
    )

    # Frequencies corresponding to the periodicities of interest
    freq_ticks = [
        1 / 180,
        1 / 20,
        1 / 7,
        1 / 10.5,
        1 / 5
    ]

    ax.set_yticks(freq_ticks)

    ax.tick_params(
        axis="y",
        labelsize=16
    )

    # Display frequency values as periodicities in base pairs
    def periodicity_formatter(x, pos):
        """Format a frequency tick as the corresponding periodicity in base pairs."""
        if x == 0:
            return ""

        if abs(x - 1 / 10.5) < 1e-4:
            return "10.5"

        return f"{int(round(1 / x))}"

    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(
            periodicity_formatter
        )
    )

    ax.set_ylabel(
        "Periodicity (bp)",
        fontsize=18
    )


def plot_stft_distributions(
    data_dict,
    dataset_group,
    output_file,
    peak_params_left=None,
    peak_params_right=None,
    verbose=False
):
    """
    Plot STFT magnitude for multiple datasets.

    For each dataset, two STFT analyses are shown:

    - Left panel: broader/global periodicity analysis.
    - Right panel: local analysis emphasizing the ~10.5 bp periodicity.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the dataset length arrays.

    dataset_group : list
        Dataset names to plot.

    output_file : str or Path
        Path for the output PNG file. A PDF is saved alongside it.

    peak_params_left : dict, optional
        STFT parameters for the left/global analysis.

    peak_params_right : dict, optional
        STFT parameters for the right/local analysis.

    verbose : bool
        If True, datasets are arranged horizontally with the global
        analysis on the top row and the local analysis on the bottom row.
    """

    # Default filtering parameters for the two STFT analyses
    if peak_params_left is None:
        peak_params_left = dict(
            cutoff_freq=0.005
        )

    if peak_params_right is None:
        peak_params_right = dict(
            cutoff_freq=0.05
        )

    num_datasets = len(dataset_group)

    # Configure the figure layout
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

    # Ensure the axes array is always two-dimensional
    if nrows == 1:
        axs = np.expand_dims(
            axs,
            axis=0
        )

    if ncols == 2 and axs.ndim == 1:
        axs = np.expand_dims(
            axs,
            axis=0
        )

    for i, dataset_name in enumerate(dataset_group):

        if verbose:
            row_global = 0
            row_local = 1
            col_left = i
            col_right = i

            # Display the colorbar label only on the final dataset
            show_legend = (
                i == num_datasets - 1
            )

        else:
            row = i // 2
            col_left = 2 * (i % 2)
            col_right = col_left + 1

            # Display the colorbar label on the rightmost panel
            show_legend = (
                col_right == ncols - 1
            )

        if verbose:

            # -------- GLOBAL PERIODICITY --------
            ax_global = axs[
                row_global,
                col_left
            ]

            process_dataset_for_stft(
                data_dict,
                dataset_name,
                ax_global,
                **peak_params_left,
                legend=show_legend,
                verbose=verbose
            )

            ax_global.set_xlim(
                0,
                1000
            )

            ax_global.set_xticks(
                [0, 500, 1000]
            )

            ax_global.set_title(
                dataset_name,
                fontsize=18,
                pad=10
            )

            # -------- LOCAL PERIODICITY --------
            ax_local = axs[
                row_local,
                col_right
            ]

            process_dataset_for_stft(
                data_dict,
                dataset_name,
                ax_local,
                **peak_params_right,
                legend=show_legend,
                verbose=verbose
            )

            ax_local.set_xlim(
                0,
                1000
            )

            ax_local.set_xticks(
                [0, 500, 1000]
            )

            # Show y-axis labels only on the first dataset
            if i == 0:
                ax_global.set_ylabel(
                    "Periodicity (bp)",
                    fontsize=18
                )
                ax_local.set_ylabel(
                    "Periodicity (bp)",
                    fontsize=18
                )
            else:
                ax_global.set_ylabel("")
                ax_local.set_ylabel("")

            # Only the lower panel receives the x-axis label
            ax_global.set_xlabel("")
            ax_local.set_xlabel(
                "Length (bp)",
                fontsize=18
            )

        else:

            # Plot global periodicity in the left panel
            process_dataset_for_stft(
                data_dict,
                dataset_name,
                axs[row, col_left],
                **peak_params_left,
                legend=False,
                verbose=verbose
            )

            # Plot local ~10.5 bp periodicity in the right panel
            process_dataset_for_stft(
                data_dict,
                dataset_name,
                axs[row, col_right],
                **peak_params_right,
                legend=show_legend,
                verbose=verbose
            )

            # Show y-axis labels only in the first column
            for col in range(ncols):
                if col == 0:
                    axs[row, col].set_ylabel(
                        "Periodicity (bp)",
                        fontsize=18
                    )
                else:
                    axs[row, col].set_ylabel("")

            # Remove x-axis labels by default
            axs[row, col_left].set_xlabel("")
            axs[row, col_right].set_xlabel("")

            # Add x-axis labels to the bottom row
            if i == num_datasets - 1 or i == num_datasets - 2:
                axs[row, col_left].set_xlabel(
                    "Length (bp)",
                    fontsize=18
                )
                axs[row, col_right].set_xlabel(
                    "Length (bp)",
                    fontsize=18
                )

            # Position the dataset name above each pair of panels
            row = i // ncols
            col_left = 0
            col_right = 1

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
                    (left_pos + right_pos) / 2
                    - 0.045
                )
                top_y = (
                    0.95
                    - i * 0.0378
                    + 0.03
                )
            else:
                center_x = (
                    (left_pos + right_pos) / 2
                    + 0.5
                    - 0.045
                )
                top_y = (
                    0.95
                    - (i - 1) * 0.0378
                    + 0.03
                )

            fig.text(
                center_x,
                top_y,
                dataset_name,
                ha="center",
                va="bottom",
                fontsize=20
            )

    # Final figure layout and output
    plt.tight_layout(
        rect=[0, 0, 1, 0.98]
    )

    plt.subplots_adjust(
        hspace=0.6
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    pdf_file = str(output_file).replace(
        ".png",
        ".pdf"
    )

    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.show()


def extract_frequency_amplitude(
    data_dict,
    dataset,
    cutoff_freq=CUTOFF_FREQ,
    fs=FS,
    nperseg=NPERSEG,
    noverlap=NOVERLAP,
    len_min=LEN_MIN,
    len_max=LEN_MAX,
    target_freq=1 / 10.5,
    threshold=1 / 2,
    n_bootstraps=100,
    normalize=True
):
    """
    Extract the maximum STFT magnitude at a target frequency.

    The analysis uses the dataset's length distribution within the
    specified length range. The distribution is high-pass filtered,
    optionally normalized, and analyzed using the STFT.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the dataset length arrays.

    dataset : str
        Dataset name.

    cutoff_freq : float
        High-pass filter cutoff frequency.

    fs : float
        Sampling frequency.

    nperseg : int
        Number of points per STFT segment.

    noverlap : int
        Number of overlapping points between STFT segments.

    len_min : int
        Minimum circle length included.

    len_max : int
        Maximum circle length included.

    target_freq : float
        Frequency at which the STFT magnitude is extracted.
        The default corresponds to a 10.5 bp periodicity.

    threshold : float
        Minimum STFT magnitude required for the result to be retained.

    n_bootstraps : int
        Number of bootstrap iterations.
        Currently retained as a parameter for compatibility.

    normalize : bool
        Whether to apply z-score normalization before the STFT.

    Returns
    -------
    sample_size : int
        Number of circles included in the length range.

    max_amplitude : float or None
        Maximum STFT magnitude at the target frequency.
        Returns None when the magnitude is below the threshold.
    """

    # Retrieve the length measurements for the dataset
    lengths = data_dict[
        f"{dataset}_lengths"
    ]

    # Restrict the analysis to the selected length range
    filtered_lengths = [
        length
        for length in lengths
        if len_min <= length < len_max
    ]

    sample_size = len(
        filtered_lengths
    )

    if sample_size == 0:
        return 0, None

    # Convert the lengths into a histogram
    bins = np.arange(
        len_min,
        len_max
    )

    length_distribution, _ = np.histogram(
        filtered_lengths,
        bins=bins
    )

    # Remove low-frequency components
    filtered = highpass_filter(
        length_distribution,
        cutoff_freq=cutoff_freq,
        fs=fs
    )

    # Normalize the signal when requested
    if normalize:
        filtered = zscore(
            filtered
        )

    # Calculate the STFT
    frequencies, time, Zxx = stft(
        filtered,
        fs=fs,
        nperseg=nperseg,
        noverlap=noverlap
    )

    # Align the STFT time axis with the original length scale
    time = time + len_min

    # Find the frequency bin closest to the target frequency
    freq_idx = np.argmin(
        np.abs(
            frequencies - target_freq
        )
    )

    # Extract the STFT magnitude at the target frequency
    amplitude_at_freq = np.abs(
        Zxx[freq_idx, :]
    )

    max_amplitude = np.max(
        amplitude_at_freq
    )

    # Retain the amplitude only if it passes the threshold
    if max_amplitude >= threshold:
        return sample_size, max_amplitude

    return sample_size, None


def plot_amplitudes_vs_samples(
    data_dict,
    target_freq=1 / 10.5,
    threshold=0.01,
    normalize=True,
    len_min=LEN_MIN,
    len_max=LEN_MAX,
):
    """
    Plot maximum STFT magnitude at the target frequency
    against the number of circles in each dataset.

    All datasets are obtained directly from the *_lengths
    entries in data_dict.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the dataset length arrays.

    target_freq : float
        Frequency of interest.

    threshold : float
        Minimum STFT magnitude included in the plot.

    normalize : bool
        Whether to apply z-score normalization before the STFT.

    len_min : int
        Minimum circle length included.

    len_max : int
        Maximum circle length included.
    """

    # Identify all datasets represented by *_lengths entries
    datasets = [
        key.removesuffix("_lengths")
        for key in data_dict
        if key.endswith("_lengths")
    ]

    sample_sizes = []
    amplitudes = []

    for dataset in datasets:

        sample_size, amplitude = extract_frequency_amplitude(
            data_dict,
            dataset,
            target_freq=target_freq,
            threshold=threshold,
            normalize=normalize,
            len_min=len_min,
            len_max=len_max
        )

        # Include only datasets whose amplitude passes the threshold
        if amplitude is not None:
            sample_sizes.append(
                sample_size
            )
            amplitudes.append(
                amplitude
            )

    # Create the amplitude-versus-dataset-size plot
    fig, ax = plt.subplots(
        figsize=(6, 5)
    )

    ax.scatter(
        sample_sizes,
        amplitudes,
        color="#450558"
    )

    ax.set_xlabel(
        "Dataset Size",
        fontsize=18
    )

    ax.set_ylabel(
        "Magnitude at 10.5 bp",
        fontsize=18
    )

    # Use a logarithmic scale for dataset size
    ax.set_xscale("log")

    ax.grid(False)

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=16
    )

    sns.despine(
        ax=ax,
        top=True,
        right=True
    )

    plt.tight_layout()

    # Create the output directory if necessary
    output_dir = Path(
        "results/STFT"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # Save the figure in PNG and PDF formats
    plt.savefig(
        output_dir / "amplitude_vs_sample_size.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        output_dir / "amplitude_vs_sample_size.pdf",
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)