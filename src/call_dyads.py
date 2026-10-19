"""Nucleosome dyad calling from MNase-seq midpoint signal.

Scores every genomic position by the fraction of fragment midpoints falling in a
central 31-bp window relative to a 201-bp background window, and reports the
local maxima of that positioning score as candidate dyads.
"""

import bisect

import numpy as np
import pyBigWig


# --- Settings ---

BIGWIG_FILE = "GSE36979_mnase_mids_combined_147.bw"
OUTPUT_FILE = "gaffney_consensus_dyads_hg18.tsv"

CHUNK_SIZE = 5000000
HALO = 1_000


REGION_THRESHOLD = 0.4
MINIMUM_PEAK_SCORE = 0.4
STRONG_PEAK_SCORE = 0.5

MINIMUM_FAR_COUNT = 10
MINIMUM_SEPARATION = 147


def window_sum(values, radius):
    """For every position, sum values within ±radius."""

    padded = np.pad(values, radius)
    cumulative = np.concatenate([[0], np.cumsum(padded)])
    width = 2 * radius + 1

    return ( cumulative[width : width + len(values)] - cumulative[: len(values)])


def find_regions(mask):
    """Find contiguous True regions as [start, end) intervals."""
    padded = np.concatenate([[False], mask, [False]])
    changes = np.flatnonzero(padded[1:] != padded[:-1])

    return zip(changes[::2], changes[1::2])


def choose_peak(score, start, end):
    """Choose the highest-scoring position in a region."""
    region_scores = score[start:end]
    maximum = region_scores.max()

    tied_positions = (np.flatnonzero(region_scores == maximum) + start)

    # If several positions tie, use the one closest to the region centre.
    region_centre = (start + end - 1) / 2

    return int(tied_positions[np.argmin(np.abs(tied_positions - region_centre))])


def remove_overlapping_calls(calls):
    """
    Retain the strongest dyad when calls are separated by less than
    MINIMUM_SEPARATION.
    """
    calls_by_strength = sorted(
        calls,
        key=lambda x: (-x["score"], -x["far_count"]),
    )

    accepted = []
    accepted_positions = []

    for call in calls_by_strength:
        position = call["dyad"]

        index = bisect.bisect_left(
            accepted_positions,
            position,
        )

        too_close_left = (
            index > 0
            and position - accepted_positions[index - 1]
            < MINIMUM_SEPARATION
        )

        too_close_right = (
            index < len(accepted_positions)
            and accepted_positions[index] - position
            < MINIMUM_SEPARATION
        )

        if too_close_left or too_close_right:
            continue

        accepted_positions.insert(index, position)
        accepted.insert(index, call)

    return accepted


def call_dyads_on_chromosome(bigwig,chrom,chrom_length):
    """Call candidate dyads on a single chromosome, processed in chunks."""
    calls = []
    for chunk_start in range(0, chrom_length, CHUNK_SIZE):
        # define the end
        chunk_end = min(chunk_start + CHUNK_SIZE,chrom_length)

        fetch_start = max(0, chunk_start - HALO)
        fetch_end = min(chrom_length, chunk_end + HALO)

        values = bigwig.values(chrom,fetch_start,fetch_end,numpy=True)
        values = np.nan_to_num(np.asarray(values, dtype=float),nan=0.0)

        # get the near and far counts for each position
        near_count = window_sum(values, radius=15)
        far_count = window_sum(values, radius=100)
        # score nucleosome positioning as the ratio of near to far counts, with a minimum of 0.0
        score = np.divide(near_count,far_count,out=np.zeros_like(near_count),where=far_count > 0)
        candidate_regions = score > REGION_THRESHOLD

        for region_start, region_end in find_regions(candidate_regions):
            peak_index = choose_peak(score,region_start,region_end)

            dyad = fetch_start + peak_index

            # Prevent duplicate calls from overlapping chunks.
            if not chunk_start <= dyad < chunk_end:
                continue

            peak_score = score[peak_index]
            peak_far_count = far_count[peak_index]

            if peak_score < MINIMUM_PEAK_SCORE:
                continue

            if peak_far_count < MINIMUM_FAR_COUNT:
                continue

            calls.append({"dyad": int(dyad),"score": float(peak_score),
                    "near_count": int(round(near_count[peak_index])),
                    "far_count": int(round(peak_far_count)),
                }
            )

    calls.sort(key=lambda x: x["dyad"])

    return remove_overlapping_calls(calls)


#open file
with pyBigWig.open(BIGWIG_FILE) as bigwig, open(
    OUTPUT_FILE,
    "w",
) as output:

    output.write(
        "chrom\tstart\tend\tpositioning_score\t"
        "near_count\tfar_count\tstrong\n"
    )
    
    total_calls = 0

    for chrom, chrom_length in bigwig.chroms().items():
        calls = call_dyads_on_chromosome(
            bigwig,
            chrom,
            chrom_length,
        )

        for call in calls:
            dyad = call["dyad"]
            score = call["score"]
            strong = int(score >= STRONG_PEAK_SCORE)

            output.write(
                f"{chrom}\t"
                f"{dyad}\t"
                f"{dyad + 1}\t"
                f"{score:.6f}\t"
                f"{call['near_count']}\t"
                f"{call['far_count']}\t"
                f"{strong}\n"
            )

        total_calls += len(calls)

        print(
            f"{chrom}: {len(calls):,} consensus dyads"
        )

print(f"\nTotal: {total_calls:,} consensus dyads")
print(f"Saved to: {OUTPUT_FILE}")
