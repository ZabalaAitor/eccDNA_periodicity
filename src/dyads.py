"""Positioning of eccDNA breakpoints relative to nucleosome dyads.

Maps each breakpoint to its closest consensus dyad, builds observed-to-expected
enrichment profiles against randomised coordinate sets, detrends the profile to
isolate the minor-groove oscillation, and compares inward- with outward-facing
rotational positions.
"""

from cmath import phase
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.stats import (shapiro, ttest_ind, mannwhitneyu)

pd.options.mode.chained_assignment = None

PROFILE_COLOR = "#450558"
MINOR_GROOVE_OUT_COLOR = "#ecd40aff"
MINOR_GROOVE_IN_COLOR = "#450558"
LABEL_FONTSIZE = 18
TICK_FONTSIZE = 16

# --- CORE FUNCTIONS ---
def plot_circle_length_phase(
    circles_bed_path,
    output_dir=None,
    circle_size=200,
):
    """
    Calculate and plot circle-length phase relative to a 10.5 bp periodicity,
    using only circles with length <= circle_size.

    Parameters
    ----------
    circles_bed_path : str
        Path to the BED file containing circle coordinates.
    output_dir : str, optional
        Path to the output directory where plots will be saved.
    circle_size : int, optional
        Maximum circle length to include in the analysis.

    Returns
    -------
    dict
        Dictionary containing the calculated lengths, phases, counts, and bin edges.
    """

    # --- 1. Load BED ---
    bed = pd.read_csv(
        circles_bed_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"]
    )

    # --- 2. Calculate circle lengths ---
    lengths = (
        bed["end"].astype(float)
        - bed["start"].astype(float)
    ).to_numpy()

    print(f"Total circles before filtering: {len(lengths):,}")

    # --- 3. Filter circles ---
    lengths = lengths[lengths <= circle_size]

    print(
        f"Circles with length <= {circle_size} bp: "
        f"{len(lengths):,}"
    )

    print(
        f"Filtered length range: "
        f"{lengths.min():.1f}–{lengths.max():.1f} bp"
    )

    # --- 4. Define periodicity ---
    period = 10.5

    # --- 5. Calculate 10.5-bp phase ---
    phases = (
        (lengths % period)
        / period
        * 360.0
    )

    print(f"Periodicity: {period:.2f} bp")

    print(
        f"Phase range: "
        f"{phases.min():.2f}–{phases.max():.2f}°"
    )

    # --- 6. Phase distribution ---
    bins = np.arange(0, 361, 1)

    counts, edges = np.histogram(
        phases,
        bins=bins
    )

    print("\n10.5-bp phase distribution:")

    for i, count in enumerate(counts):
        print(
            f"{edges[i]:.1f}–{edges[i+1]:.1f}°: "
            f"{count}"
        )

    # --- 7. Plot ---
    fig, ax = plt.subplots(
        figsize=(7,5)
    )

    ax.bar(
        edges[:-1],
        counts,
        width=5,
        align="edge",
        color=PROFILE_COLOR
    )

    ax.set_xlabel(
        "Circle length phase relative to 10.5 bp (degrees)", 
        fontsize=LABEL_FONTSIZE
    )

    ax.set_ylabel(
        "Number of circles",
        fontsize=LABEL_FONTSIZE
    )

    # ax.set_title(
    #     f"Circle length phase relative to 10.5-bp periodicity "
    #     f"(length ≤ {max_length} bp)"
    # )

    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 30))

    ax.grid(
        axis="y",
        alpha=0.3
    )

    ax.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE,
    )

    sns.despine(
        ax=ax
    )

    plt.tight_layout()

    # --- 8. Save ---
    if output_dir is not None:
        plt.savefig(
            output_dir / "circle_length_phase.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.savefig(
            output_dir / "circle_length_phase.pdf",
            dpi=300,
            bbox_inches="tight"
        )

        print(
            f"\nPlot saved to: {output_dir}"
        )

    plt.show()

    return {
        "lengths": lengths,
        "phases": phases,
        "counts": counts,
        "bin_edges": edges,
        "period": period,
    }

def find_closest_consensus_dyad(circles, nucleosomes_dyads):
    """
    Find the closest consensus dyad to each circle breakpoint.

    Parameters
    ----------
    circles : pandas.DataFrame
        DataFrame containing circle coordinates.
    nucleosomes_dyads : pandas.DataFrame
        DataFrame containing nucleosome dyad coordinates.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the closest consensus dyad for each circle breakpoint.
    """

    unique_chromosomes = np.intersect1d(
        circles["chrom"].unique(),
        nucleosomes_dyads["chrom"].unique(),
    )

    chromosomes_to_keep = [f"chr{i}" for i in range(1, 23)]

    unique_chromosomes = np.intersect1d(
        unique_chromosomes,
        chromosomes_to_keep,
    )

    output = []

    for chrom in unique_chromosomes:

        circles_chrom = circles[
            circles["chrom"] == chrom
        ]

        dyad_positions = (
            nucleosomes_dyads.loc[
                nucleosomes_dyads["chrom"] == chrom,
                "start",
            ]
            .sort_values()
            .to_numpy()
        )

        starts = circles_chrom["start"].to_numpy()
        ends = circles_chrom["end"].to_numpy()

        def closest_dyads(breakpoints):

            """Return the closest consensus dyad to each breakpoint, with signed distance."""
            idx = np.searchsorted(
                dyad_positions,
                breakpoints,
            )

            left_idx = np.clip(
                idx - 1,
                0,
                len(dyad_positions) - 1,
            )

            right_idx = np.clip(
                idx,
                0,
                len(dyad_positions) - 1,
            )

            left_dist = np.abs(
                breakpoints -
                dyad_positions[left_idx]
            )

            right_dist = np.abs(
                breakpoints -
                dyad_positions[right_idx]
            )

            closest_idx = np.where(
                right_dist < left_dist,
                right_idx,
                left_idx,
            )

            return dyad_positions[closest_idx]

        # Start breakpoint
        output.append(
            pd.DataFrame({
                "circle_chrom":
                    circles_chrom["chrom"].to_numpy(),
                "circle_start": starts,
                "circle_end": ends,
                "breakpoint": starts,
                "dyad": closest_dyads(starts),
            })
        )

        # End breakpoint
        output.append(
            pd.DataFrame({
                "circle_chrom":
                    circles_chrom["chrom"].to_numpy(),
                "circle_start": starts,
                "circle_end": ends,
                "breakpoint": ends,
                "dyad": closest_dyads(ends),
            })
        )

    return pd.concat(
        output,
        ignore_index=True,
    )

def plot_non_centered_breakpoints(
    dyad_breakpoint_df,
    output_prefix,
    smooth_window_bp=5,
):
    """
    Plot breakpoint distance from closest dyad.

    Parameters
    ----------
    dyad_breakpoint_df : pandas.DataFrame
        DataFrame containing dyad and breakpoint information.
    output_prefix : str
        Prefix for the output file names.
    smooth_window_bp : int, optional
        Window size for smoothing the distribution.

    Returns
    -------
    None
    """

    df = dyad_breakpoint_df.copy()

    distance = (
        df["breakpoint"] -
        df["dyad"]
    )

    counts = (
        distance
        .value_counts()
        .sort_index()
    )

    # Remove distances larger than 1000 bp
    counts = counts[
        np.abs(counts.index.to_numpy()) <= 1000
    ]

    # Fold mirrored distances
    counts = (
        counts
        .groupby(np.abs(counts.index))
        .sum()
        .sort_index()
    )

    # Smooth
    counts = counts.rolling(
        smooth_window_bp,
        center=True,
        min_periods=1,
    ).mean()

    # Full plot
    plt.figure(figsize=(7, 5))

    plt.plot(
        counts.index,
        counts.values,
    )

    plt.xlabel("Distance from Dyad (bp)")
    plt.ylabel("Count")

    for x in [5, 15, 25, 35, 45, 55, 65, 73]:
        plt.axvline(
            x=x,
            linestyle="--",
            alpha=0.5,
        )

    plt.tight_layout()

    plt.savefig(
        f"{output_prefix}_non_centered.png",
        dpi=300,
    )
    plt.savefig(
        f"{output_prefix}_non_centered.pdf",
        dpi=300,
    )
    plt.show()
    plt.close()

    # Zoomed plot
    plt.figure(figsize=(7, 5))

    plt.plot(
        counts.index,
        counts.values,
    )

    plt.xlabel("Distance from Dyad (bp)")
    plt.ylabel("Count")

    for x in [5, 15, 25, 35, 45, 55, 65, 73]:
        plt.axvline(
            x=x,
            linestyle="--",
            alpha=0.5,
        )

    plt.xlim(-100, 100)

    plt.tight_layout()

    plt.savefig(
        f"{output_prefix}_non_centered_zoomed.png",
        dpi=300,
    )
    plt.savefig(
        f"{output_prefix}_non_centered_zoomed.pdf",
        dpi=300,
    )
    plt.show()
    plt.close()

def process_permutation(
    bed_file,
    consensus_dyads,
    output_dir,
):
    """
    Process one permutation BED file.
    
    Parameters
    ----------
    bed_file : pathlib.Path
        Path to the permutation BED file.
    consensus_dyads : pandas.DataFrame
        DataFrame containing consensus dyad coordinates.
    output_dir : pathlib.Path
        Directory for saving the output.

    Returns
    -------
    pathlib.Path
        Path to the output TSV file.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_tsv = (
        output_dir /
        f"{bed_file.stem}.tsv"
    )

    circles = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
        ],
    )

    out = find_closest_consensus_dyad(
        circles,
        consensus_dyads,
    )

    out.to_csv(
        output_tsv,
        sep="\t",
        index=False,
    )

    return output_tsv

# --- MAIN ANALYSIS ---

def analyze_circle_breakpoints(
    consensus_dyads_path,
    circles_bed_path,
    output_dir,
    permutations_dir=None,
    n_jobs=1,
    positioning_score_threshold=0.5,
    smooth_window_bp=5,
    run_permutations=True,
):
    """
    Run the circle breakpoint / consensus dyad analysis.

    Parameters
    ----------
    consensus_dyads_path : str or Path
        Consensus dyads TSV file.

    circles_bed_path : str or Path
        BED file containing circle coordinates.

    output_dir : str or Path
        Directory for analysis results.

    permutations_dir : str or Path, optional
        Directory containing permutation BED files.

    n_jobs : int
        Number of parallel processes for permutations.

    positioning_score_threshold : float
        Minimum positioning score for consensus dyads.

    smooth_window_bp : int
        Rolling window used for smoothing.

    run_permutations : bool
        Whether to analyse permutation BED files.

    Returns
    -------
    dict
        Dictionary containing the main analysis results.
    """

    consensus_dyads_path = Path(
        consensus_dyads_path
    )

    circles_bed_path = Path(
        circles_bed_path
    )

    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --- Load consensus dyads ---

    consensus_dyads = pd.read_csv(
        consensus_dyads_path,
        sep="\t",
    )

    consensus_dyads = consensus_dyads[
        consensus_dyads["positioning_score"]

    >= positioning_score_threshold

    ].copy()

    print(
        f"Consensus dyads after filtering: "
        f"{len(consensus_dyads):,}"
    )

    # --- Load circles ---

    circles = pd.read_csv(
        circles_bed_path,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
        ],
    )

    print(
        f"Circles: {len(circles):,}"
    )

    # --- Find closest dyads ---

    breakpoint_dyads = find_closest_consensus_dyad(
        circles,
        consensus_dyads,
    )

    print(
        f"Breakpoint-dyad pairs: "
        f"{len(breakpoint_dyads):,}"
    )

    # --- Save breakpoint results ---

    breakpoint_output = (
        output_dir /
        "closest_consensus_dyads_to_circle_breakpoints.tsv"
    )

    breakpoint_dyads.to_csv(
        breakpoint_output,
        sep="\t",
        index=False,
    )

    print(
        f"Saved: {breakpoint_output}"
    )

    # --- Plot ---

    plot_prefix = (
        output_dir /
        "distance_from_consensus_dyad_to_circle_breakpoints"
    )

    plot_non_centered_breakpoints(
        breakpoint_dyads,
        plot_prefix,
        smooth_window_bp=smooth_window_bp,
    )

    print(
        f"Plots saved with prefix: {plot_prefix}"
    )

    # --- Permutations ---

    permutation_results = []

    if run_permutations and permutations_dir is not None:

        permutations_dir = Path(
            permutations_dir
        )

        permutation_output_dir = (
            output_dir /
            "random_closest_consensus_dyads"
        )

        bed_files = sorted(
            permutations_dir.glob(
                "permutation_*.bed"
            )
        )

        print(
            f"Found {len(bed_files)} permutation files."
        )

        if n_jobs == 1:

            for i, bed_file in enumerate(
                bed_files,
                1,
            ):

                result = process_permutation(
                    bed_file,
                    consensus_dyads,
                    permutation_output_dir,
                )

                permutation_results.append(
                    result
                )

                print(
                    f"[{i}/{len(bed_files)}] "
                    f"Done: {result.name}"
                )

        else:

            with ProcessPoolExecutor(
                max_workers=n_jobs
            ) as executor:

                futures = {
                    executor.submit(
                        process_permutation,
                        bed_file,
                        consensus_dyads,
                        permutation_output_dir,
                    ): bed_file
                    for bed_file in bed_files
                }

                for i, future in enumerate(
                    as_completed(futures),
                    1,
                ):

                    result = future.result()

                    permutation_results.append(
                        result
                    )

                    print(
                        f"[{i}/{len(bed_files)}] "
                        f"Done: {result.name}"
                    )

    return {
        "consensus_dyads": consensus_dyads,
        "circles": circles,
        "breakpoint_dyads": breakpoint_dyads,
        "breakpoint_output": breakpoint_output,
        "permutation_results": permutation_results,
    }

def _load_breakpoint_distances(
    filepath,
    circle_size,
    max_distance=500,
):
    """
    Load breakpoint/dyad distances from one TSV file.
    
    Parameters
    ----------
    filepath : pathlib.Path
        Path to the TSV file containing breakpoint/dyad distances.
    circle_size : int
        Maximum size of circles to consider.
    max_distance : int, optional
        Maximum distance to consider.

    Returns
    -------
    pandas.Series
        Series containing the counts of distances.
    """

    df = pd.read_csv(
        filepath,
        sep="\t",
    )

    # Keep short circles
    df = df[
        (df["circle_end"] - df["circle_start"])
        <= circle_size
    ]

    distance = (
        df["dyad"] -
        df["breakpoint"]
    )

    distance = distance[
        distance.isin(
            range(
                -max_distance,
                max_distance + 1,
            )
        )
    ]

    return distance.value_counts().sort_index()

def _calculate_obs_exp(
    observed_distance,
    permutation_distances,
):
    """
    Calculate observed/expected ratio.
    
    Parameters
    ----------
    observed_distance : pandas.Series
        Series containing the observed distances.
    permutation_distances : dict
        Dictionary containing the expected distances for each permutation.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the observed/expected ratios.
    """

    obs_exp = {}

    for filename, expected_distance in (
        permutation_distances.items()
    ):
        ratio = (
            observed_distance /
            expected_distance
        )

        obs_exp[filename] = ratio

    return pd.DataFrame(obs_exp)

def _plot_colored_profile(
    x,
    y,
    output_file,
    xlabel,
    ylabel,
    xlim=(-60, 60),
    selected_points=None,
    width=2.5,
):
    """
    Plot profile with highlighted groove positions.
    
    Parameters
    ----------
    x : numpy.ndarray
        Array containing the x-coordinates.
    y : numpy.ndarray
        Array containing the y-coordinates.
    output_file : pathlib.Path
        Path to the output file.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    xlim : tuple, optional
        Limits for the x-axis.
    selected_points : numpy.ndarray, optional
        Array containing the selected points.
    width : float, optional
        Width of the highlighted segments.

    Returns
    -------
    None
    """

    if selected_points is None:
        selected_points = np.array(
            [5, 15, 25, 35, 45, 55],
            dtype=float,
        )

        selected_points = np.concatenate([
            selected_points,
            -selected_points,
        ])

    points = np.column_stack([
        x,
        y,
    ])

    segments = np.stack([
        points[:-1],
        points[1:],
    ], axis=1)

    segment_midpoints = (
        x[:-1] +
        x[1:]
    ) / 2

    highlight_segments = np.any(
        np.abs(
            segment_midpoints[:, None]
            - selected_points[None, :]
        ) <= width,
        axis=1,
    )

    colors = np.where(
        highlight_segments,
        "C3",
        "C2",
    )

    fig, ax = plt.subplots()

    collection = LineCollection(
        segments,
        colors=colors,
        linewidths=1.5,
        capstyle="projecting",
        antialiased=True,
    )

    ax.add_collection(collection)

    ax.autoscale_view(
        scalex=False,
        scaley=True,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color="C2",
            label="Minor groove inward oriented",
        ),
        Line2D(
            [0],
            [0],
            color="C3",
            label="Minor groove outward oriented",
        ),
    ]

    ax.legend(
        handles=legend_handles
    )

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        output_file.replace(".png", ".pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

# --- PLOT COLORED PROFILE ---

def _plot_colored_profile(
    x,
    y,
    output_file,
    xlabel,
    ylabel,
    xlim=(-60, 60),
    selected_points=None,
    width=2.5,
):
    """
    Plot profile with highlighted minor groove positions.
    
    Parameters
    ----------
    x : numpy.ndarray
        Array containing the x-coordinates.
    y : numpy.ndarray
        Array containing the y-coordinates.
    output_file : pathlib.Path
        Path to the output file.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    xlim : tuple, optional
        Limits for the x-axis.
    selected_points : numpy.ndarray, optional
        Array containing the selected points.
    width : float, optional
        Width of the highlighted segments.

    Returns
    -------
    None
    """

    if selected_points is None:
        selected_points = np.array(
            [5, 15, 25, 35, 45, 55],
            dtype=float,
        )

        selected_points = np.concatenate([
            selected_points,
            -selected_points,
        ])

    points = np.column_stack([
        x,
        y,
    ])

    segments = np.stack([
        points[:-1],
        points[1:],
    ], axis=1)

    segment_midpoints = (
        x[:-1] +
        x[1:]
    ) / 2

    highlight_segments = np.any(
        np.abs(
            segment_midpoints[:, None]
            - selected_points[None, :]
        ) <= width,
        axis=1,
    )

    segment_colors = np.where(
        highlight_segments,
        MINOR_GROOVE_OUT_COLOR,
        PROFILE_COLOR,
    )

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    collection = LineCollection(
        segments,
        colors=segment_colors,
        linewidths=1.8,
        capstyle="projecting",
        antialiased=True,
    )

    ax.add_collection(collection)

    ax.autoscale_view(
        scalex=False,
        scaley=True,
    )

    ax.set_xlabel(
        xlabel,
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_ylabel(
        ylabel,
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_xlim(*xlim)

    ax.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE,
    )

    sns.despine(
        ax=ax
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        output_file.replace(".png", ".pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    plt.close(fig)

# --- MAIN ENRICHMENT ANALYSIS ---

def analyze_circle_breakpoint_enrichment(
    observed_file,
    permutations_dir,
    output_dir,
    circle_size=200,
    max_distance=500,
    smooth_window=21,
    circles_bed_path=None
):
    """
    Calculate observed/expected enrichment of consensus dyads
    around circle breakpoints.

    Parameters
    ----------
    observed_file : str or Path
        Path to the TSV file containing observed breakpoint/dyad distances.
    permutations_dir : str or Path
        Path to the directory containing permutation files.
    output_dir : str or Path
        Path to the output directory.
    circle_size : int, optional
        Size of the circles to analyze.
    max_distance : int, optional
        Maximum distance to consider.
    smooth_window : int, optional
        Window size for smoothing.
    circles_bed_path : str or Path, optional
        Path to the BED file containing circle coordinates.
    """

    observed_file = Path(observed_file)
    permutations_dir = Path(permutations_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --- 1. CIRCLE LENGTH NUCLEOSOMAL PHASE ---

    print(
        "Plotting circle length nucleosomal phase..."
    )

    circle_phase = plot_circle_length_phase(
        circles_bed_path,
        output_dir,
        circle_size=circle_size
    )
    # --- OBSERVED ---

    print("Loading observed data...")

    observed_distance = _load_breakpoint_distances(
        observed_file,
        circle_size=circle_size,
        max_distance=max_distance,
    )

    print(
        f"Observed breakpoint pairs: "
        f"{observed_distance.sum():,}"
    )

    # --- PERMUTATIONS ---

    permutation_files = sorted(
        permutations_dir.glob("*.tsv")
    )

    print(
        f"Found {len(permutation_files)} "
        "permutation files."
    )

    expected_distances = {}

    for permutation_file in permutation_files:

        print(
            f"Loading: {permutation_file.name}"
        )

        expected_distances[
            permutation_file.name
        ] = _load_breakpoint_distances(
            permutation_file,
            circle_size=circle_size,
            max_distance=max_distance,
        )

    # --- OBSERVED / EXPECTED ---

    obs_exp_df = _calculate_obs_exp(
        observed_distance,
        expected_distances,
    )

    mean_obs_exp = (
        obs_exp_df
        .mean(axis=1)
        .sort_index()
    )

    # --- SMOOTHING ---

    # Smoothing for the full profile
    smoothed_obs_exp = (
        mean_obs_exp
        .rolling(
            window=21,
            center=True,
            min_periods=1,
        )
        .mean()
    )

    # Smoothing for the zoomed profile
    smoothed_obs_exp_zoom = (
        mean_obs_exp
        .rolling(
            window=5,
            center=True,
            min_periods=1,
        )
        .mean()
    )

    # --- SAVE O/E TABLE ---

    obs_exp_output = (
        output_dir /
        "observed_expected_ratio.tsv"
    )

    obs_exp_df.to_csv(
        obs_exp_output,
        sep="\t",
    )

    mean_obs_exp_output = (
        output_dir /
        "mean_observed_expected_ratio.tsv"
    )

    mean_obs_exp.to_csv(
        mean_obs_exp_output,
        sep="\t",
        header=["observed_expected_ratio"],
    )

    # --- COMMON POSITIONS ---

    selected_points = np.array(
        [5, 15, 25, 35, 45, 55],
        dtype=float,
    )

    selected_points = np.concatenate([
        selected_points,
        -selected_points,
    ])

    width = 2.5

    # --- UNFOLDED PROFILE WITH CONFIDENCE INTERVAL ---

    x = mean_obs_exp.index.to_numpy(
        dtype=float
    )

    # Mean observed/expected ratio
    y = smoothed_obs_exp.to_numpy(
        dtype=float
    )

    # --- Confidence interval from permutation distributions ---

    # Calculate 95% confidence interval across permutations
    ci_lower = (
        obs_exp_df
        .quantile(
            0.025,
            axis=1,
        )
        .sort_index()
    )

    ci_upper = (
        obs_exp_df
        .quantile(
            0.975,
            axis=1,
        )
        .sort_index()
    )

    # Smooth the confidence limits with the same
    # window used for the main profile
    ci_lower_smoothed = (
        ci_lower
        .rolling(
            window=smooth_window,
            center=True,
            min_periods=1,
        )
        .mean()
    )

    ci_upper_smoothed = (
        ci_upper
        .rolling(
            window=smooth_window,
            center=True,
            min_periods=1,
        )
        .mean()
    )

    ci_low = ci_lower_smoothed.to_numpy(
        dtype=float
    )

    ci_high = ci_upper_smoothed.to_numpy(
        dtype=float
    )

    # --- PLOT ---

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    # --- Nucleosome region ---

    ax.axvspan(
        -83.5,
        83.5,
        color="#56B4E9",
        alpha=0.20,
        zorder=0,
    )

    # --- 95% confidence interval ---

    ax.fill_between(
        x,
        ci_low,
        ci_high,
        color="#BDBDBD",
        alpha=0.35,
        linewidth=0,
        label="95% CI",
        zorder=1,
    )

    # --- Observed / expected profile ---

    ax.plot(
        x,
        y,
        color=PROFILE_COLOR,
        linewidth=1.8,
        label="Observed/Expected",
        zorder=2,
    )

    # --- Labels ---

    ax.set_xlabel(
        "Breakpoint distance from dyad (bp)",
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_ylabel(
        "Observed/Expected Ratio",
        fontsize=LABEL_FONTSIZE,
    )

    ax.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE,
    )

    # --- Legend ---

    ax.legend(
        handles=[
            Patch(
                facecolor="#56B4E9",
                alpha=0.20,
                edgecolor="none",
                label="Nucleosome + linker region (±83.5 bp)",
            ),
            Patch(
                facecolor="#BDBDBD",
                alpha=0.35,
                edgecolor="none",
                label="95% CI",
            ),
            Line2D(
                [0],
                [0],
                color=PROFILE_COLOR,
                linewidth=1.8,
                label="Observed/Expected",
            ),
        ],
        loc="upper right",
        frameon=False,
        fontsize=TICK_FONTSIZE,
    )

    # --- Style ---

    sns.despine(
        ax=ax
    )

    fig.tight_layout()

    # --- SAVE ---

    unfolded_output = (
        output_dir /
        "obs_exp_ratio.png"
    )

    fig.savefig(
        unfolded_output,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        unfolded_output.with_suffix(".pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # ZOOMED PROFILE
    # SMOOTHING = 5
    # ========================================================

    y_zoom = smoothed_obs_exp_zoom.to_numpy(
        dtype=float
    )

    points = np.column_stack([
        x,
        y_zoom,
    ])

    segments = np.stack([
        points[:-1],
        points[1:],
    ], axis=1)

    segment_midpoints = (
        x[:-1] +
        x[1:]
    ) / 2

    highlight_segments = np.any(
        np.abs(
            segment_midpoints[:, None]
            - selected_points[None, :]
        ) <= width,
        axis=1,
    )

    segment_colors = np.where(
        highlight_segments,
        MINOR_GROOVE_IN_COLOR,
        MINOR_GROOVE_OUT_COLOR,
    )

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    collection = LineCollection(
        segments,
        colors=segment_colors,
        linewidths=1.8,
        capstyle="projecting",
        antialiased=True,
    )

    ax.add_collection(
        collection
    )

    ax.autoscale_view(
        scalex=False,
        scaley=True,
    )

    ax.set_xlabel(
        "Breakpoint distance from dyad (bp)",
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_ylabel(
        "Observed/Expected Ratio",
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_xlim(
        -60,
        60,
    )

    ax.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=MINOR_GROOVE_IN_COLOR,
            linewidth=2,
            label="Minor groove inward oriented",
        ),
        Line2D(
            [0],
            [0],
            color=MINOR_GROOVE_OUT_COLOR,
            linewidth=2,
            label="Minor groove outward oriented",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        frameon=False,
        fontsize=TICK_FONTSIZE,
    )

    sns.despine(
        ax=ax
    )

    fig.tight_layout()

    zoom_output = (
        output_dir /
        "obs_exp_ratio_zoom.png"
    )

    fig.savefig(
        zoom_output,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        zoom_output.with_suffix(".pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close(fig)

    # --- DETRENDING ---

    y_raw = mean_obs_exp.to_numpy(
        dtype=float
    )

    baseline = savgol_filter(
        y_raw,
        window_length=21,
        polyorder=2,
    )

    y_smoothed = savgol_filter(
        y_raw,
        window_length=5,
        polyorder=2,
    )

    residual = (
        y_smoothed -
        baseline
    )

    residual_blurred = gaussian_filter1d(
        residual,
        sigma=1.5,
    )

    # --- DETRENDED PROFILE ---

    points = np.column_stack([
        x,
        residual_blurred,
    ])

    segments = np.stack([
        points[:-1],
        points[1:],
    ], axis=1)

    segment_midpoints = (
        x[:-1] +
        x[1:]
    ) / 2

    highlight_segments = np.any(
        np.abs(
            segment_midpoints[:, None]
            - selected_points[None, :]
        ) <= width,
        axis=1,
    )

    segment_colors = np.where(
        highlight_segments,
        MINOR_GROOVE_IN_COLOR,
        MINOR_GROOVE_OUT_COLOR,
    )

    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    collection = LineCollection(
        segments,
        colors=segment_colors,
        linewidths=1.8,
        capstyle="projecting",
        antialiased=True,
    )

    ax.add_collection(
        collection
    )

    ax.autoscale_view(
        scalex=False,
        scaley=True,
    )

    ax.set_xlabel(
        "Distance (bp)",
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_ylabel(
        "Deviation from baseline",
        fontsize=LABEL_FONTSIZE,
    )

    ax.set_xlim(
        -60,
        60,
    )

    ax.axhline(
        0,
        color="black",
        linewidth=1,
        linestyle="--",
    )

    ax.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=MINOR_GROOVE_IN_COLOR,
            linewidth=2,
            label="Minor groove inward oriented",
        ),
        Line2D(
            [0],
            [0],
            color=MINOR_GROOVE_OUT_COLOR,
            linewidth=2,
            label="Minor groove outward oriented",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        frameon=False,
        fontsize=TICK_FONTSIZE,
    )

    sns.despine(
        ax=ax
    )

    fig.tight_layout()

    detrended_output = (
        output_dir /
        "detrended.png"
    )

    fig.savefig(
        detrended_output,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        detrended_output.with_suffix(".pdf"),
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
    plt.close(fig)

    # --- MINOR GROOVE IN / OUT ---

    minor_groove_in = [
        p + delta
        for p in selected_points
        for delta in (-2, -1, 0, 1, 2)
    ]

    minor_groove_in = np.asarray(
        minor_groove_in,
        dtype=int,
    )

    minor_groove_out = [
        0,
        10,
        20,
        30,
        40,
        50,
    ]

    minor_groove_out += [
        -p
        for p in minor_groove_out
    ]

    minor_groove_out = [
        p + delta
        for p in minor_groove_out
        for delta in (-2, -1, 0, 1, 2)
    ]

    minor_groove_out = np.asarray(
        minor_groove_out,
        dtype=int,
    )

    x_int = x.astype(int)

    range_mask = (
        (x_int >= -60) &
        (x_int <= 60)
    )

    in_mask = (
        np.isin(
            x_int,
            minor_groove_in,
        )
        & range_mask
    )

    out_mask = (
        np.isin(
            x_int,
            minor_groove_out,
        )
        & range_mask
    )

    minor_groove_in_values = (
        residual_blurred[in_mask]
    )

    minor_groove_out_values = (
        residual_blurred[out_mask]
    )

    # --- VIOLIN PLOT ---

    data = [
        minor_groove_in_values,
        minor_groove_out_values,
    ]

    fig, ax = plt.subplots(
        figsize=(6, 5)
    )

    parts = ax.violinplot(
        data,
        positions=[1, 2],
        widths=0.7,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    # --- VIOLIN COLORS ---

    violin_colors = [
        MINOR_GROOVE_IN_COLOR,
        MINOR_GROOVE_OUT_COLOR,
    ]

    for body, color in zip(
        parts["bodies"],
        violin_colors,
    ):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.25)
        body.set_linewidth(1.2)

    # --- MEDIAN COLORS ---

    if "cmedians" in parts:
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(2.0)

    # --- JITTERED DATA POINTS ---

    rng = np.random.default_rng(42)

    jitter_in = rng.uniform(
        -0.08,
        0.08,
        size=len(minor_groove_in_values),
    )

    jitter_out = rng.uniform(
        -0.08,
        0.08,
        size=len(minor_groove_out_values),
    )

    # Minor groove in
    ax.scatter(
        np.full(
            len(minor_groove_in_values),
            1,
        ) + jitter_in,
        minor_groove_in_values,
        color=MINOR_GROOVE_IN_COLOR,
        s=20,
        edgecolors="none",
    )

    # Minor groove out
    ax.scatter(
        np.full(
            len(minor_groove_out_values),
            2,
        ) + jitter_out,
        minor_groove_out_values,
        color=MINOR_GROOVE_OUT_COLOR,
        s=20,
        edgecolors="none",
    )

    # --- AXIS ---

    ax.set_xticks([1, 2])

    ax.set_xticklabels(
        [
            "Minor groove \ninward oriented",
            "Minor groove \noutward oriented",
        ],
        fontsize=TICK_FONTSIZE,
    )

    ax.set_ylabel(
        "Deviation from baseline",
        fontsize=LABEL_FONTSIZE,
    )

    ax.tick_params(
        axis="y",
        labelsize=TICK_FONTSIZE,
    )

    sns.despine(
        ax=ax
    )

    fig.tight_layout()

    # --- SAVE PNG + PDF ---

    violin_output = (
        output_dir /
        "obs_exp_ratio_folded_detrended_boxplot.png"
    )

    fig.savefig(
        violin_output,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        violin_output.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    plt.show()
    plt.close(fig)

        # --- STATISTICAL TESTS ---

    print("\n" + "=" * 60)
    print("STATISTICAL COMPARISON: MINOR GROOVE ORIENTATION")
    print("=" * 60)

    in_values = np.asarray(
        minor_groove_in_values,
        dtype=float,
    )

    out_values = np.asarray(
        minor_groove_out_values,
        dtype=float,
    )

    # Remove NaNs / inf
    in_values = in_values[
        np.isfinite(in_values)
    ]

    out_values = out_values[
        np.isfinite(out_values)
    ]

    print(
        f"Inward:  n = {len(in_values)}"
    )

    print(
        f"Outward: n = {len(out_values)}"
    )

    print(
        f"Inward median:  {np.median(in_values):.5f}"
    )

    print(
        f"Outward median: {np.median(out_values):.5f}"
    )

    print(
        f"Inward mean:    {np.mean(in_values):.5f}"
    )

    print(
        f"Outward mean:   {np.mean(out_values):.5f}"
    )

    # --- NORMALITY ---

    # Shapiro-Wilk is not reliable for very large samples.
    # If n > 5000, use a random subset of 5000 observations.

    rng_stats = np.random.default_rng(42)

    if len(in_values) > 5000:
        in_normality = rng_stats.choice(
            in_values,
            size=5000,
            replace=False,
        )
    else:
        in_normality = in_values

    if len(out_values) > 5000:
        out_normality = rng_stats.choice(
            out_values,
            size=5000,
            replace=False,
        )
    else:
        out_normality = out_values

    shapiro_in = shapiro(
        in_normality
    )

    shapiro_out = shapiro(
        out_normality
    )

    print("\n--- Normality ---")

    print(
        f"Inward:  W = {shapiro_in.statistic:.4f}, "
        f"p = {shapiro_in.pvalue:.4e}"
    )

    print(
        f"Outward: W = {shapiro_out.statistic:.4f}, "
        f"p = {shapiro_out.pvalue:.4e}"
    )

    # --- GROUP COMPARISON ---

    if (
        shapiro_in.pvalue > 0.05
        and shapiro_out.pvalue > 0.05
    ):

        # --- Welch's t-test ---

        test = ttest_ind(
            in_values,
            out_values,
            equal_var=False,
        )

        test_name = "Welch's t-test"

        statistic = test.statistic
        pvalue = test.pvalue

        # Cohen's d
        mean_in = np.mean(in_values)
        mean_out = np.mean(out_values)

        var_in = np.var(
            in_values,
            ddof=1,
        )

        var_out = np.var(
            out_values,
            ddof=1,
        )

        pooled_sd = np.sqrt(
            (
                var_in +
                var_out
            ) / 2
        )

        # Same group order as the t-test above (in_values, out_values), so that
        # the sign of Cohen's d matches the sign of the t statistic.
        cohens_d = (
            mean_in -
            mean_out
        ) / pooled_sd

        effect_name = "Cohen's d"
        effect_size = cohens_d

    else:

        # --- Mann–Whitney U ---

        test = mannwhitneyu(
            in_values,
            out_values,
            alternative="two-sided",
        )

        test_name = "Mann–Whitney U"

        statistic = test.statistic
        pvalue = test.pvalue

        # Rank-biserial correlation
        n1 = len(in_values)
        n2 = len(out_values)

        rank_biserial = (
            2 * statistic /
            (n1 * n2)
            - 1
        )

        effect_name = "Rank-biserial r"
        effect_size = rank_biserial

    # --- STATISTICAL REPORT ---

    print("\n--- Group comparison ---")

    print(
        f"Test: {test_name}"
    )

    print(
        f"Statistic = {statistic:.4f}"
    )

    print(
        f"p-value = {pvalue:.4e}"
    )

    print(
        f"{effect_name} = {effect_size:.4f}"
    )

    print("=" * 60)

    # --- SAVE STATISTICS ---

    statistics_output = (
        output_dir /
        "minor_groove_statistics.tsv"
    )

    statistics_df = pd.DataFrame({
        "comparison": [
            "Minor groove inward vs outward"
        ],
        "n_inward": [
            len(in_values)
        ],
        "n_outward": [
            len(out_values)
        ],
        "median_inward": [
            np.median(in_values)
        ],
        "median_outward": [
            np.median(out_values)
        ],
        "mean_inward": [
            np.mean(in_values)
        ],
        "mean_outward": [
            np.mean(out_values)
        ],
        "shapiro_W_inward": [
            shapiro_in.statistic
        ],
        "shapiro_p_inward": [
            shapiro_in.pvalue
        ],
        "shapiro_W_outward": [
            shapiro_out.statistic
        ],
        "shapiro_p_outward": [
            shapiro_out.pvalue
        ],
        "test": [
            test_name
        ],
        "statistic": [
            statistic
        ],
        "p_value": [
            pvalue
        ],
        "effect_size_name": [
            effect_name
        ],
        "effect_size": [
            effect_size
        ],
    })

    statistics_df.to_csv(
        statistics_output,
        sep="\t",
        index=False,
    )

    print(
        f"Statistics saved to: "
        f"{statistics_output}"
    )

    return {
        "observed_distance": observed_distance,
        "expected_distances": expected_distances,
        "obs_exp": obs_exp_df,
        "mean_obs_exp": mean_obs_exp,
        "smoothed_obs_exp": smoothed_obs_exp,
        "residual": residual,
        "residual_blurred": residual_blurred,
        "minor_groove_in": minor_groove_in_values,
        "minor_groove_out": minor_groove_out_values,
        "statistics": statistics_df,
        "output_dir": output_dir,
    }