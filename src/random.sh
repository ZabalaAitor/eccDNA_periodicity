#!/usr/bin/env bash
set -euo pipefail

# Generate randomised eccDNA coordinate sets for the observed/expected profile.
# Usage: bash src/random.sh [dyads_dir] [reference_dir] [n_permutations] [n_jobs]
#
# Each permutation is seeded with its own index, so the whole set is reproducible.
# Circles stay on their original chromosome (-chrom) and ENCODE blacklist regions
# are excluded. Note that the random stream of `bedtools shuffle` is not guaranteed
# to be stable across bedtools versions.

DYADS_DIR="${1:-results/dyads}"
REFERENCE_DIR="${2:-data/reference}"
N_PERMUTATIONS="${3:-1000}"
N_JOBS="${4:-8}"

CIRCLES_BED="${DYADS_DIR}/non_genic_circles.bed"
GENOME_FILE="${REFERENCE_DIR}/hg38.genome"
BLACKLIST_BED="${REFERENCE_DIR}/ENCFF356LFX.bed"
RANDOM_DIR="${DYADS_DIR}/random"

mkdir -p "$RANDOM_DIR"
seq 1 "$N_PERMUTATIONS" | xargs -P "$N_JOBS" -I {} sh -c \
    "bedtools shuffle \
        -i '${CIRCLES_BED}' \
        -g '${GENOME_FILE}' \
        -excl '${BLACKLIST_BED}' \
        -chrom \
        -seed {} \
        > '${RANDOM_DIR}/permutation_{}.bed'"

ls "$RANDOM_DIR" | wc -l
