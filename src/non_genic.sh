#!/usr/bin/env bash
set -euo pipefail

# Keep only eccDNA that do not overlap an annotated gene.
# Usage: bash src/non_genic.sh [dyads_dir]
#
# Any overlap of at least one base pair with a GENCODE gene-level feature removes
# the circle, so the retained set is intergenic under a conservative definition.

DYADS_DIR="${1:-results/dyads}"

CIRCLES_BED="${DYADS_DIR}/merged_hg38.bed"
GENIC_BED="${DYADS_DIR}/hg38_genic_regions.bed"
OUTPUT_BED="${DYADS_DIR}/non_genic_circles.bed"

bedtools intersect -v \
    -a "$CIRCLES_BED" \
    -b "$GENIC_BED" \
    > "$OUTPUT_BED"

wc -l "$OUTPUT_BED"
