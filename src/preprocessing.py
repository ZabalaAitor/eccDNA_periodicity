from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def merge_bed_files(file_list, output_file):
    """"
    Merge multiple BED files into a single file.

    Parameters
    ----------
    file_list : list of str or Path
        List of BED file paths to merge.

    output_file : str or Path
        Path to the output merged BED file.
    """
    with open(output_file, "w") as outfile:
        for f in file_list:
            with open(f) as infile:
                for line in infile:
                    outfile.write(line)

def prepare_and_plot_dataset_sizes(
    base_dir="data",
    len_min=0,
    len_max=1000,
    output_dir="results/dataset_size",
    figsize=(8, 8),
    palette=None,
):
    """
    Process BED files, calculate dataset sizes, assign protocols,
    sort datasets by protocol and size, and generate the plot.

    Parameters
    ----------
    base_dir : str or Path
        Base directory containing the dataset folders.

    len_min : int or float
        Minimum eccDNA length to retain.

    len_max : int or float
        Maximum eccDNA length to retain.

    output_dir : str or Path
        Directory where PNG and PDF plots are saved.

    figsize : tuple
        Figure size.

    palette : dict or None
        Seaborn palette for protocols.

    Returns
    -------
    df : pandas.DataFrame
        Dataset names, sizes, and protocols in plotting order.

    data_dict : dict
        Filtered data and length arrays for all datasets.
    """

    # ========================================================
    # Configuration
    # ========================================================

    base_dir = Path(base_dir)
    output_dir = Path(output_dir)

    folders = [
        "arabidopsis",
        "chicken",
        "human",
        "mice",
        "mice1L3KO",
    ]

    aggregated_folders = [
        "whole",
        "microDNA",
        "Circle-Seq",
        "ATAC-seq",
    ]

    # Aggregated datasets excluded from the main plot
    excluded_datasets = {
        "whole",
        "microDNA",
        "Circle-Seq",
        "ATAC-seq",
    }

    # Protocol order in the plot
    protocol_order = [
        "MD",
        "ATAC-seq",
        "Circle-Seq",
        "Other",
    ]

    # Default palette
    if palette is None:
        palette = {
            "Circle-Seq": "#4C72B0",
            "ATAC-seq": "#DD8452",
            "MD": "#55A868",
            "Other": "#999999",
        }

    # ========================================================
    # 1. Read and filter datasets
    # ========================================================

    data_dict = {}

    for folder in folders + aggregated_folders:

        folder_path = base_dir / folder

        if not folder_path.exists():
            continue

        for bed_file in folder_path.glob("*.bed"):

            dataset = bed_file.stem

            data = pd.read_csv(
                bed_file,
                sep="\t",
                header=None,
                low_memory=False,
            )

            # Calculate eccDNA length
            data["length"] = data[2] - data[1]

            # Filter by length
            filtered = data.loc[
                data["length"].between(
                    len_min,
                    len_max,
                )
            ].copy()

            # Store filtered data
            data_dict[f"{dataset}_filtered_data"] = filtered

            # Store lengths
            data_dict[f"{dataset}_lengths"] = (
                filtered["length"].to_numpy()
            )

    # ========================================================
    # 2. Create dataset-size dataframe
    # ========================================================

    records = [
        {
            "Dataset": key.removesuffix("_lengths"),
            "Size": len(values),
        }
        for key, values in data_dict.items()
        if key.endswith("_lengths")
        and key.removesuffix("_lengths") not in excluded_datasets
    ]

    df = pd.DataFrame(records)

    # ========================================================
    # 3. Assign protocol
    # ========================================================

    def get_protocol(dataset):

        if dataset.endswith("-MD"):
            return "MD"

        if dataset.endswith("-AS"):
            return "ATAC-seq"

        if dataset.endswith("-CS"):
            return "Circle-Seq"

        return "Other"

    df["Protocol"] = df["Dataset"].map(get_protocol)

    # ========================================================
    # 4. Set protocol order
    # ========================================================

    df["Protocol"] = pd.Categorical(
        df["Protocol"],
        categories=protocol_order,
        ordered=True,
    )

    # ========================================================
    # 5. Sort by protocol and size
    # ========================================================

    df = (
        df.sort_values(
            ["Protocol", "Size"],
            ascending=[True, False],
            kind="stable",
        )
        .reset_index(drop=True)
    )

    # ========================================================
    # 6. Plot
    # ========================================================

    fig, ax = plt.subplots(
        figsize=figsize
    )

    sns.barplot(
        data=df,
        y="Dataset",
        x="Size",
        hue="Protocol",
        dodge=False,
        palette=palette,
        order=df["Dataset"],
        ax=ax,
    )

    ax.set_xscale("log")

    ax.set_xlabel(
        "#eccDNAs (log scale)",
        fontsize=18,
    )

    ax.set_ylabel("")

    ax.tick_params(
        axis="x",
        labelsize=16,
    )

    ax.tick_params(
        axis="y",
        labelsize=16,
    )

    sns.despine(ax=ax)

    ax.get_legend().remove()

    fig.tight_layout()

    # ========================================================
    # 7. Save
    # ========================================================

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output_dir / "dataset_sizes_per_sample_ordered.png",
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        output_dir / "dataset_sizes_per_sample_ordered.pdf",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ========================================================
    # 8. Return results
    # ========================================================

    return df, data_dict